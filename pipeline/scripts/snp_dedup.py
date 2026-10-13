#!/usr/bin/env python3
"""
Dereplication by core-SNP distance. v2 fixes two defects in v1.

DEFECT 1 -- SINGLE LINKAGE CHAINS ONCE THE THRESHOLD EXCEEDS ZERO
v1 used union-find, i.e. single linkage. At t=0 that is safe because
identity is transitive. At t=25 it is not: A--B=25, B--C=25, C--D=25 chains
into one cluster of diameter 75, silently. Since measurement showed that
independent assemblies of one strain differ by 3-21 SNPs, t=0 is untenable
and single linkage goes with it.

v2 uses COMPLETE linkage: every pair inside a cluster is within the
threshold, so cluster diameter is bounded. "The same isolate sequenced
twice" should mean mutually near-identical, not connected by a chain.

DEFECT 2 -- THE HARD CAP DISCARDED THE MOST DISTINCTIVE GENOMES
v1 sorted by descending cluster size and kept the top N, so when the cap
bit it preferentially kept representatives of large outbreak clusters and
dropped singletons. Singletons are usually the most divergent genomes;
large clusters are the over-sampling the cap exists to control. The
behaviour was exactly backwards for a phylogenetic study. Tier A (251) never
triggers it, so it would first have bitten on Tier B (3,441).

v2 uses greedy farthest-point (maximin) selection seeded with the lineage
anchors, which keeps the spread of diversity rather than the bulk of the
sampling.

ANI IS STILL NOT USED: all Y. pestis pairs exceed 99.99% ANI, so ANI has no
resolving power here.

Thresholds apply to the RAW distance matrix (observed differences among
jointly called sites), not the coverage-scaled one.
"""
import csv
import os
import re
import sys
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

snakemake = globals().get("snakemake")  # noqa: F821

SENTINEL = 10 ** 6


def pick_representative(members, accs, meta, priority, median_len):
    def key(i):
        acc = accs[i]
        m = meta.get(acc, {})
        score = []
        for rule in priority:
            if rule == "refseq_gcf":
                score.append(0 if acc.startswith("GCF") else 1)
            elif rule == "has_collection_date":
                score.append(0 if m.get("collection_date") else 1)
            elif rule == "has_geo_loc_name":
                score.append(0 if m.get("geo_loc_name") else 1)
            elif rule == "closest_to_median_length":
                try:
                    score.append(abs(int(m.get("total_bp", 0)) - median_len))
                except ValueError:
                    score.append(10 ** 9)
        score.append(acc)
        return tuple(score)

    return sorted(members, key=key)[0]


def farthest_point(D, candidates, seeds, k):
    """Greedy maximin: repeatedly take the candidate farthest from the set."""
    chosen = list(seeds)
    if not chosen:
        chosen = [candidates[0]]
    pool = [c for c in candidates if c not in chosen]
    finite = np.where(D >= SENTINEL, np.nan, D).astype(float)
    mind = np.nanmin(finite[np.ix_(pool, chosen)], axis=1) if pool else None
    while pool and len(chosen) < k:
        j = int(np.nanargmax(mind))
        pick = pool.pop(j)
        chosen.append(pick)
        mind = np.delete(mind, j)
        if pool:
            mind = np.fmin(mind, finite[np.ix_(pool, [pick])].ravel())
    return chosen


def main():
    raw_path = Path(str(snakemake.input.dist))
    D_full = np.load(raw_path).astype(float)
    order = [ln.strip() for ln in
             raw_path.with_suffix(".order.txt").read_text().splitlines()
             if ln.strip()]

    cand = list(csv.DictReader(open(snakemake.input.cand), delimiter="\t"))
    cand_acc = {r["accession"] for r in cand}
    meta = {r["accession"]: r for r in
            csv.DictReader(open(snakemake.input.meta), delimiter="\t")}

    sel = [i for i, a in enumerate(order) if a in cand_acc]
    accs = [order[i] for i in sel]
    D = D_full[np.ix_(sel, sel)]
    D = np.minimum(D, SENTINEL)
    np.fill_diagonal(D, 0)
    D = (D + D.T) / 2.0

    # ----------------------------------------------------------------
    # DELETION-PROFILE CONSTRAINT
    #
    # SNP distance is computed only over JOINTLY CALLED sites, so a large
    # deletion contributes exactly zero to it: the deleted positions are
    # uncalled in the deleted genome and never enter the count. A pgm-
    # derivative and its pgm+ parent therefore look identical, merge at t=20,
    # and the representative rule (refseq > dated > ...) silently decides
    # whether that lineage is recorded as pgm+ or pgm- in the PAV matrix --
    # which is the signal Fig.2 and Fig.3 measure.
    #
    # Measured 2026-10-10 on Tier A: 25 of 36 multi-member clusters had mixed
    # deletion profiles; cluster 108 (size 12) held 7 genomes with no large
    # deletion and 4 carrying the 105 kb pgm deletion.
    #
    # Genomes with different deletion profiles are therefore forced apart.
    profile_path = getattr(snakemake.input, "profile", None)
    n_del_blocked = 0
    n_unprofiled = -1
    blocked_pairs = []          # (acc_a, acc_b, raw_d) for the whitelist subset
    # SY_DEDUP_PROFILE_MODE: full (default) | off | threestate
    #   off        no profile constraint at all  (the 175 / 36 / 25 baseline)
    #   threestate three-state-aware constraint: chromosome-block differences
    #              always block; a plasmid difference blocks ONLY when BOTH
    #              genomes have that plasmid resolved (neither UNRESOLVED),
    #              because non-submission is not evidence of loss.
    mode = os.environ.get("SY_DEDUP_PROFILE_MODE", "full")
    if profile_path and mode != "off":
        prof = {}
        with open(profile_path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                prof[r["accession"]] = (
                    r.get("replicons_absent", "") or "",
                    r.get("big_blocks_5kb", "") or "",
                )

        def _norm(a):
            return re.sub(r"\.\d+$", "", str(a).strip())

        def _chrom(a):
            rec = prof.get(a)
            if rec is None:
                return None
            return frozenset(b for b in rec[1].split(";") if b.startswith("NC_003143:"))

        tstate = None
        if mode in ("threestate", "conservative"):
            tstate = {}
            ts_path = os.environ.get(
                "SY_DEDUP_THREESTATE",
                "results/qc/replicon_three_state_seqreport.tsv")
            for r in csv.DictReader(open(ts_path), delimiter="\t"):
                tstate[_norm(r["accession"])] = (r["pCD1_state"], r["pPCP1_state"],
                                                 r["pMT1_state"])

        sig = [prof.get(a) for a in accs]
        n_unprofiled = sum(1 for s in sig if s is None)
        for i in range(len(accs)):
            if sig[i] is None:
                continue
            for j in range(i + 1, len(accs)):
                if sig[j] is None:
                    continue
                if mode == "threestate":
                    block = (_chrom(accs[i]) != _chrom(accs[j]))
                    if not block:
                        ta = tstate.get(_norm(accs[i]))
                        tb = tstate.get(_norm(accs[j]))
                        if ta and tb:
                            block = any(ta[k] != "UNRESOLVED" and tb[k] != "UNRESOLVED"
                                        and ta[k] != tb[k] for k in range(3))
                elif mode == "conservative":
                    # missing evidence blocks, it does not merge: an UNRESOLVED
                    # plasmid is NOT evidence of sameness (same asymmetry as
                    # the polarity rule). Two UNRESOLVED stay equal and do not
                    # block; any other difference blocks.
                    ta = tstate.get(_norm(accs[i]))
                    tb = tstate.get(_norm(accs[j]))
                    block = (_chrom(accs[i]) != _chrom(accs[j])) or (
                        ta is not None and tb is not None and ta != tb)
                else:
                    block = sig[i] != sig[j]
                if block:
                    # capture the SNP distance BEFORE masking it, so the
                    # control can assert "distance alone said merge"
                    d_raw = float(D[i, j])
                    D[i, j] = D[j, i] = SENTINEL
                    n_del_blocked += 1
                    blocked_pairs.append((accs[i], accs[j], d_raw))

    whitelist = {v["accession"] for v in snakemake.params.whitelist.values()}

    # Persist the whitelist subset of the blocked pairs. The full list is
    # ~5.2M pairs on Tier B and useless to store or eyeball, but the anchors
    # are 8 accessions, so this is a few dozen rows and it contains exactly
    # the control cases. scripts/sanity_duplicate_control.py consumes it as
    # its fourth and strongest condition:
    #   distance says merge -> it did not merge -> this pair was blocked by
    #   the deletion profile, named.
    wl_pairs = []
    if profile_path:
        wl_pairs = [x for x in blocked_pairs
                    if x[0] in whitelist or x[1] in whitelist]
        wl_path = Path(snakemake.output.clusters).with_name(
            Path(snakemake.output.clusters).name.replace(
                "_clusters.tsv", "_blocked_whitelist_pairs.tsv"))
        with open(wl_path, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["accession_a", "accession_b", "raw_d",
                        "whitelist_member"])
            for a_, b_, d_ in sorted(wl_pairs, key=lambda x: (x[0], x[1])):
                w.writerow([a_, b_, f"{d_:.0f}",
                            a_ if a_ in whitelist else b_])

    t = float(snakemake.params.max_dist)
    method = getattr(snakemake.params, "linkage", "complete")

    if t <= 0:
        raise SystemExit(
            "FATAL: max_snp_distance <= 0. Measurement showed independent "
            "assemblies of one strain differ by 3-21 SNPs, so a zero "
            "threshold leaves duplicate records in the set. Run "
            "scripts/dedup_threshold_sweep.py and set the threshold it "
            "selects."
        )

    Z = linkage(squareform(D, checks=False), method=method)
    labels = fcluster(Z, t=t, criterion="distance")

    clusters = {}
    for i, lab in enumerate(labels):
        clusters.setdefault(int(lab), []).append(i)

    lengths = [int(meta.get(a, {}).get("total_bp", 0) or 0) for a in accs]
    median_len = int(np.median(lengths)) if lengths else 0

    keep_idx, cluster_rows, max_diam = [], [], 0.0
    for cid, members in sorted(clusters.items()):
        anchors = [i for i in members if accs[i] in whitelist]
        rep = (sorted(anchors)[0] if anchors
               else pick_representative(members, accs, meta,
                                        snakemake.params.priority, median_len))
        keep_idx.append(rep)
        for a in anchors:
            if a != rep and a not in keep_idx:
                keep_idx.append(a)
        if len(members) > 1:
            sub = D[np.ix_(members, members)]
            max_diam = max(max_diam, float(sub[sub < SENTINEL].max()))
        cluster_rows.append({
            "cluster_id": cid, "size": len(members),
            "representative": accs[rep],
            "diameter": int(D[np.ix_(members, members)].max())
            if len(members) > 1 else 0,
            "members": ",".join(sorted(accs[i] for i in members))})

    if max_diam > 2 * t:
        raise SystemExit(
            f"FATAL: max cluster diameter {max_diam:.0f} exceeds 2x the "
            f"threshold ({t}). Complete linkage bounds diameter at the "
            f"threshold, so this indicates a broken distance matrix, not a "
            f"borderline case. Do not proceed."
        )

    cap = getattr(snakemake.params, "hard_cap", None)
    capped = False
    if cap and len(keep_idx) > int(cap):
        cap = int(cap)
        seeds = [i for i in keep_idx if accs[i] in whitelist]
        keep_idx = farthest_point(D, keep_idx, seeds, cap)
        capped = True

    keep = [accs[i] for i in keep_idx]

    # ANCHOR INVARIANT: patches at five separate reduction points all failed to
    # keep the anchors in (GCA/GCF pairing dropped KIM10+; the residue-based
    # called_frac rule would have dropped Nepal516 and KIM10+). A sixth bypass
    # would just move the problem, so the check is an invariant instead.
    from _invariants import assert_anchors_present
    assert_anchors_present(keep, snakemake.params.whitelist, "snp_dedup")

    with open(snakemake.output.keep, "w") as fh:
        fh.write("accession\n")
        for a in keep:
            fh.write(a + "\n")
    with open(snakemake.output.clusters, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["cluster_id", "size",
                                           "representative", "diameter",
                                           "members"], delimiter="\t")
        w.writeheader()
        w.writerows(cluster_rows)

    msg = "\n".join([
        f"[snp_dedup] input={len(accs)} threshold={t} linkage={method}",
        f"[snp_dedup] clusters={len(clusters)} kept={len(keep)} "
        f"collapsed={len(accs) - len(keep)}",
        f"[snp_dedup] max cluster diameter={max_diam:.0f} "
        f"(bound {2 * t:.0f})",
        f"[snp_dedup] profile mode={mode} blocked pairs={n_del_blocked} "
        f"unprofiled={n_unprofiled}",
        f"[snp_dedup] blocked pairs involving whitelist anchors="
        f"{len(wl_pairs)} -> ..._blocked_whitelist_pairs.tsv",
        f"[snp_dedup] largest cluster="
        f"{max(len(m) for m in clusters.values())}",
        f"[snp_dedup] anchors retained="
        f"{sum(1 for a in keep if a in whitelist)}/{len(whitelist)}",
        f"[snp_dedup] hard cap applied={capped}"
        + (" (greedy farthest-point, anchor-seeded)" if capped else ""),
    ])
    print(msg, file=sys.stderr)
    Path(snakemake.log[0]).write_text(msg + "\n")


if __name__ == "__main__":
    main()
