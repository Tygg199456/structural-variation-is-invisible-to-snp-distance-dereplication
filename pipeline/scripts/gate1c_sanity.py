#!/usr/bin/env python3
"""
Positive control for the Tier A dereplication chain. v2 recalibrates three
expectations that v1 got wrong.

WHAT v1 ASSUMED, AND WHAT MEASUREMENT SHOWED
--------------------------------------------
  called_frac ~ 1.0 on complete genomes
      WRONG. The CO92 reference compared against ITSELF calls 95.5% of
      chromosome positions. Shredding an assembly into 250 bp pseudo-reads
      cannot map IS-element copies uniquely, and Y. pestis carries
      IS100/IS285/IS1541/IS1661 across roughly 3-4% of the chromosome.
      0.955 is this pipeline's ceiling, not 1.0.
      v2 measures each genome against that empirical ceiling.

  duplicate assemblies at distance 0
      WRONG. Independent assemblies of one strain differ by 3-21 SNPs
      (median 7) -- different sequencing eras, different assemblers,
      base-level disagreement in uniquely mappable regions.
      v2 requires distance <= the dereplication threshold, and co-clustering.

  retained 190-225
      WRONG, because it followed from the t=0 assumption.
      v2 takes the expected count from the threshold sweep instead of from
      a guess.

v2 also adds a chaining check, because complete linkage is only trustworthy
if its diameter bound actually holds.

Usage:
    python3 scripts/gate1c_sanity.py \
        results/dedup/snp_distance_raw_tierA.npy \
        results/dedup/genotypes_tierA.npz \
        results/dedup/tierA_clusters.tsv \
        results/dedup/tierA_keep.tsv \
        --threshold 25 --expected-retained 212 --tolerance 15
"""
import argparse
import csv
import sys
from pathlib import Path

import numpy as np

SENTINEL = 10 ** 6

DUPLICATE_ASSEMBLIES = {
    "CO92":       ["GCF_000009065.1", "GCF_001293415.1"],
    "Antiqua":    ["GCF_000013825.1", "GCF_000834825.1"],
    "PestoidesF": ["GCF_000016445.1", "GCF_000834315.1"],
    "A1122":      ["GCF_000222975.1", "GCF_000834755.1"],
    "Angola":     ["GCF_000018805.1", "GCF_000834845.1"],
}

# Empirical ceiling: reference vs itself through this pipeline.
# Override with --self-ceiling if re-measured.
SELF_COMPARISON_CEILING = 0.955
# A genome must reach this fraction OF THE CEILING, not of 1.0.
MIN_FRAC_OF_CEILING = 0.90


def norm(acc):
    return acc.split("_", 1)[1] if "_" in acc else acc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dist")
    ap.add_argument("npz")
    ap.add_argument("clusters")
    ap.add_argument("keep")
    ap.add_argument("--threshold", type=float, required=True,
                    help="dereplication threshold chosen by the sweep")
    ap.add_argument("--expected-retained", type=int, required=True,
                    help="retained count the sweep reported for that threshold")
    ap.add_argument("--tolerance", type=int, default=10)
    ap.add_argument("--self-ceiling", type=float,
                    default=SELF_COMPARISON_CEILING)
    ap.add_argument("--deletion-profile", default=None,
                    help="per-genome deletion profile TSV (optional); genomes "
                         "with real large deletions are reported, never failed")
    a = ap.parse_args()

    D = np.load(a.dist).astype(float)
    order = [ln.strip() for ln in
             Path(a.dist).with_suffix(".order.txt").read_text().splitlines()
             if ln.strip()]
    idx = {norm(x): i for i, x in enumerate(order)}

    z = np.load(a.npz, allow_pickle=True)
    called = z["called_frac"] if "called_frac" in z else None
    samples = [s.decode() if isinstance(s, bytes) else str(s)
               for s in z["samples"]]

    clusters = list(csv.DictReader(open(a.clusters), delimiter="\t"))
    keep = [ln.strip() for ln in
            Path(a.keep).read_text().splitlines()[1:] if ln.strip()]

    failures, notes, warnings = [], [], []
    floor = a.self_ceiling * MIN_FRAC_OF_CEILING

    # ---- 1. called fraction ------------------------------------------------
    # REWRITTEN 2026-10-10. The previous version FAILed any genome below 90% of
    # the self-comparison ceiling (0.859). Measured against reality that flag
    # selected a lineage, not a defect: all 12 genomes it caught were Russian
    # or Nepalese, two of them lineage anchors, and uncalled_region_profile
    # showed 47 of the 50 lowest-callability genomes carry REAL deletions
    # (105 kb pgm locus, whole pMT1 / pCD1), i.e. the phenomenon this study
    # measures. A reference-based callability metric is collinear with the
    # signal, so it cannot be an absolute filter.
    #
    # Final form: report the distribution; warn on the lowest decile; FAIL only
    # below min_joint_called_frac, where a genome cannot support its own
    # placement in the tree.
    HARD_FLOOR = 0.70
    if called is not None:
        dec = np.percentile(called, [0, 10, 25, 50, 75, 90, 100])
        print(f"[1] called_frac  n={len(called)}  "
              f"min={called.min():.4f} median={np.median(called):.4f} "
              f"max={called.max():.4f}")
        print(f"    deciles 0/10/25/50/75/90/100 = "
              + "/".join(f"{v:.3f}" for v in dec))
        print(f"    pipeline ceiling (ref vs itself) = {a.self_ceiling:.3f}")

        hard = [(samples[i], float(called[i]))
                for i in range(len(called)) if called[i] < HARD_FLOOR]
        d10 = float(dec[1])
        low_decile = [(samples[i], float(called[i]))
                      for i in range(len(called)) if called[i] < d10]

        if hard:
            failures.append(
                f"{len(hard)} genomes below the hard floor {HARD_FLOOR}: "
                f"{hard[:5]}\n"
                f"    Below {HARD_FLOOR} a genome cannot support its own "
                f"placement in the tree.")
        else:
            print(f"    PASS -- every genome above the hard floor "
                  f"{HARD_FLOOR} (min observed {called.min():.4f})")
        if low_decile:
            warnings.append(
                f"lowest decile (called_frac < {d10:.3f}): {len(low_decile)} "
                f"genomes, e.g. {low_decile[:5]} -- reported, not failed")
    else:
        notes.append("called_frac absent from npz; check 1 skipped")

    # genomes whose uncalled regions are real deletions are FINDINGS, not faults.
    if getattr(a, "deletion_profile", None):
        try:
            import csv as _csv
            n_bio = 0
            with open(a.deletion_profile) as fh:
                for r in _csv.DictReader(fh, delimiter="\t"):
                    if int(r.get("n_big") or 0) > 0 or (r.get("replicons_absent") or ""):
                        n_bio += 1
            print(f"[1b] genomes with real large deletions in the profile: "
                  f"{n_bio} -- reported as findings, never failed")
        except Exception as exc:   # profile is optional
            notes.append(f"deletion profile unreadable: {exc}")

    # ---- 2. blocked pairs -------------------------------------------------
    off = ~np.eye(D.shape[0], dtype=bool)
    n_blocked = int(((D >= SENTINEL) & off).sum() // 2)
    print(f"[2] blocked pairs = {n_blocked}")
    if n_blocked:
        failures.append(
            f"{n_blocked} Tier A pairs blocked by the joint-called guard. "
            f"On complete genomes with ~0.95 called fraction, joint overlap "
            f"is ~0.90 and nothing should be blocked at 0.70.")
    else:
        print("    PASS")

    # ---- 3. duplicate assemblies ------------------------------------------
    print(f"[3] same-strain control (threshold {a.threshold:.0f})")
    member_of = {}
    for row in clusters:
        for m in row["members"].split(","):
            member_of[norm(m.strip())] = row["cluster_id"]

    for strain, accs in DUPLICATE_ASSEMBLIES.items():
        present = [x for x in accs if norm(x) in idx]
        if len(present) < 2:
            notes.append(f"{strain}: {len(present)}/2 assemblies present "
                         f"(the dnadiff gate may have removed one)")
            continue
        i, j = idx[norm(present[0])], idx[norm(present[1])]
        d = float(D[i, j])
        same = (member_of.get(norm(present[0])) ==
                member_of.get(norm(present[1])))
        ok = d <= a.threshold and same
        print(f"    {strain:<12} d={d:<7.0f} same_cluster={str(same):<6} "
              f"[{'PASS' if ok else 'FAIL'}]")
        if not ok:
            failures.append(
                f"{strain}: d={d:.0f}, merged={same}. Expected d <= "
                f"{a.threshold:.0f} and merged. If d is far above the "
                f"threshold the sweep picked the wrong value; if d is small "
                f"but merged is False, clustering is broken.")

    # ---- 4. chaining ------------------------------------------------------
    diams = [int(r.get("diameter", 0) or 0) for r in clusters]
    max_diam = max(diams) if diams else 0
    print(f"[4] max cluster diameter = {max_diam} (bound {2 * a.threshold:.0f})")
    if max_diam > 2 * a.threshold:
        failures.append(
            f"diameter {max_diam} exceeds twice the threshold. Complete "
            f"linkage bounds this; exceeding it means single linkage is "
            f"still in use somewhere, or the matrix is asymmetric.")
    else:
        print("    PASS -- no chaining")

    # ---- 5. retained count vs the sweep -----------------------------------
    lo = a.expected_retained - a.tolerance
    hi = a.expected_retained + a.tolerance
    print(f"[5] retained = {len(keep)}  (sweep predicted "
          f"{a.expected_retained} +/- {a.tolerance})")
    if not (lo <= len(keep) <= hi):
        failures.append(
            f"retained {len(keep)}, outside {lo}-{hi}. The sweep and the "
            f"dedup run disagree, so they are not using the same matrix, "
            f"threshold or linkage.")
    else:
        print("    PASS")

    print()
    for n in notes:
        print(f"NOTE  {n}")
    for w in warnings:
        print(f"WARN  {w}")
    if failures:
        print()
        for f in failures:
            print(f"FAIL  {f}")
        print("\nDo not proceed to the next stage until these are resolved.")
        sys.exit(1)
    print("\nAll hard checks passed.")
    if warnings:
        print(f"{len(warnings)} warning(s) above are reported, not failures: "
              f"reference-based callability is a covariate, not a filter.")


if __name__ == "__main__":
    main()
