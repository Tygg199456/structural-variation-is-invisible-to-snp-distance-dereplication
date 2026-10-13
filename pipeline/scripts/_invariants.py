#!/usr/bin/env python3
"""
Pipeline invariants, asserted at the point of writing rather than patched
filter by filter.

WHY THIS MODULE EXISTS
----------------------
Anchor strains have now been lost, or nearly lost, at five separate
reduction points:

    1. GCA/GCF pairing            KIM10+ deleted as a "mirror"
    2. metadata QC thresholds     (bypass added)
    3. dnadiff gate               (bypass added)
    4. tier_b_coverage            (bypass omitted from the contract)
    5. called_frac threshold      Nepal516 flagged for removal

Adding a bypass at the Nth filter guarantees an Nth+1 filter will miss it.
The eight whitelisted genomes are the phylogenetic anchors and the basis of
comparability with Morelli 2010 / Cui 2013; losing one silently invalidates
lineage assignment downstream.

So instead: every rule that writes a Tier A or Tier B subset calls
`assert_anchors_present()` before it finishes. Where the filtering logic
lives stops mattering.

    from _invariants import assert_anchors_present
    ...
    assert_anchors_present(kept_accessions, config["whitelist"], "tier_b_coverage")
"""
import sys


def _acc_set(whitelist):
    return {v["accession"] if isinstance(v, dict) else str(v)
            for v in whitelist.values()}


def assert_anchors_present(kept, whitelist, step_name, hard=True):
    """
    Raise (or warn) if any whitelisted anchor is missing from `kept`.

    kept       iterable of accessions surviving this step
    whitelist  config["whitelist"] mapping name -> {"accession": ...}
    step_name  the reduction point, named in the error
    hard       False downgrades to a warning, for steps where an anchor
               legitimately may not qualify -- dating taxa selection, where
               an anchor without a collection_date cannot be a tip and
               forcing it in would fabricate a date.
    """
    kept = set(kept)
    anchors = _acc_set(whitelist)
    missing = sorted(anchors - kept)
    if not missing:
        return

    names = {v["accession"]: k for k, v in whitelist.items()
             if isinstance(v, dict)}
    detail = ", ".join(f"{names.get(a, '?')} ({a})" for a in missing)
    msg = (f"ANCHOR INVARIANT VIOLATED at '{step_name}': "
           f"{len(missing)}/{len(anchors)} anchors dropped -> {detail}\n"
           f"  Anchors are the phylogenetic reference points and the basis "
           f"of comparability with the published population scheme. Every "
           f"reduction step must exempt them.\n"
           f"  Fix the filter, do not relax the invariant.")
    if hard:
        raise SystemExit(msg)
    print(f"WARNING: {msg}", file=sys.stderr)


def assert_no_reference_metric_filter(step_name, metric_name, used_as_filter):
    """
    Guard the principle that reference-based metrics are covariates, not
    filters.

    Any metric computed against a single reference -- dnadiff AlignedBases,
    called_frac, mapped fraction -- measures divergence from that reference
    as much as it measures quality. In a study of divergence and genome
    decay it is collinear with the signal. Three such metrics have each
    flagged the same Russian / Central Asian lineage.

    Genuine quality outliers are identified by RESIDUAL from the divergence
    trend (scripts/reference_bias_report.py), never by an absolute cut.
    """
    if used_as_filter:
        raise SystemExit(
            f"DESIGN VIOLATION at '{step_name}': '{metric_name}' is computed "
            f"against a single reference and is being used as an exclusion "
            f"filter.\n"
            f"  Such metrics are divergence measures in disguise and are "
            f"collinear with the signal this study measures.\n"
            f"  Record it as a covariate, and identify true outliers by "
            f"residual from the divergence trend "
            f"(scripts/reference_bias_report.py)."
        )


def summarise_reduction(step_name, before, after, log_path=None):
    """One uniform line per reduction point, so the audit trail is complete."""
    line = (f"[reduction] {step_name}: {before} -> {after} "
            f"({before - after} removed, {100 * after / max(before, 1):.1f}% kept)")
    print(line, file=sys.stderr)
    if log_path:
        with open(log_path, "a") as fh:
            fh.write(line + "\n")
    return line
