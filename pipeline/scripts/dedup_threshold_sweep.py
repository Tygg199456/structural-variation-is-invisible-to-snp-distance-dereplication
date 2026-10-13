#!/usr/bin/env python3
"""
Choose the dereplication threshold from data instead of from a guess.

WHY THIS EXISTS
---------------
`max_snp_distance: 0` was a number I wrote without evidence, on the
assumption that two assemblies of one strain are identical. Measurement says
otherwise: independent assemblies of CO92, Antiqua, Pestoides F, A1122 and
Angola differ by 3-21 SNPs, and the reference compared against ITSELF calls
only 95.5% of sites because IS elements (IS100/IS285/IS1541/IS1661) defeat
unique mapping of shredded pseudo-reads.

So the pipeline has a technical noise floor, and the dereplication threshold
has to sit above it. This script measures where, by sweeping thresholds and
reporting what each one does.

SELECTION RULE -- fixed here in advance, not chosen after seeing the output:

    Take the SMALLEST threshold t such that
      (1) all 5 known same-strain assembly pairs merge,
      (2) max_cluster_diameter <= 2 * t,
      (3) max_cluster_size <= 25.

  (1) the noise floor must be crossed, or dereplication fails its purpose
  (2) complete linkage bounds diameter at t; 2*t allows for scaling slack,
      and anything beyond that means the implementation is wrong
  (3) a cluster above 25 is no longer "the same isolate twice" but genuine
      epidemiological diversity being thinned, which needs its own argument

COMPLETE LINKAGE, NOT SINGLE LINKAGE
------------------------------------
At t=0 single linkage is safe: identity is transitive. At t=25 it is not.
A--B=25, B--C=25, C--D=25 chains into one cluster whose diameter is 75, and
nothing reports it. "Same strain" should mean every pair in the cluster is
within the noise floor, which is complete linkage.

Usage:
    python3 scripts/dedup_threshold_sweep.py \
        results/dedup/snp_distance_raw_tierA.npy \
        results/dedup/tierA_refined.tsv
"""
import csv
import sys
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

THRESHOLDS = [0, 3, 5, 10, 15, 20, 25, 30, 50, 75, 100]
SENTINEL = 10 ** 6

DUPLICATE_ASSEMBLIES = {
    "CO92":       ["GCF_000009065.1", "GCF_001293415.1"],
    "Antiqua":    ["GCF_000013825.1", "GCF_000834825.1"],
    "PestoidesF": ["GCF_000016445.1", "GCF_000834315.1"],
    "A1122":      ["GCF_000222975.1", "GCF_000834755.1"],
    "Angola":     ["GCF_000018805.1", "GCF_000834845.1"],
}

MAX_DIAMETER_FACTOR = 2
MAX_CLUSTER_SIZE = 25


def norm(acc):
    return acc.split("_", 1)[1] if "_" in acc else acc


def main():
    dist_path = sys.argv[1]
    cand_path = sys.argv[2] if len(sys.argv) > 2 else None

    D = np.load(dist_path).astype(float)
    order = [ln.strip() for ln in
             Path(dist_path).with_suffix(".order.txt").read_text().splitlines()
             if ln.strip()]

    if cand_path:
        keep = {norm(r["accession"]) for r in
                csv.DictReader(open(cand_path), delimiter="\t")}
        sel = [i for i, a in enumerate(order) if norm(a) in keep]
        D = D[np.ix_(sel, sel)]
        order = [order[i] for i in sel]

    n = len(order)
    idx = {norm(a): i for i, a in enumerate(order)}
    D = np.minimum(D, SENTINEL)
    np.fill_diagonal(D, 0)
    D = (D + D.T) / 2.0                      # enforce symmetry for squareform

    Z = linkage(squareform(D, checks=False), method="complete")

    dup_idx = {}
    for strain, accs in DUPLICATE_ASSEMBLIES.items():
        present = [norm(a) for a in accs if norm(a) in idx]
        if len(present) == 2:
            dup_idx[strain] = (idx[present[0]], idx[present[1]])

    print(f"genomes = {n}   known same-strain pairs available = {len(dup_idx)}")
    print(f"observed same-strain distances: "
          f"{ {s: int(D[i, j]) for s, (i, j) in dup_idx.items()} }")
    print()
    hdr = (f"{'t':>5} {'clusters':>9} {'retained':>9} {'singl':>6} "
           f"{'maxDiam':>8} {'maxSize':>8} {'dupMerged':>10} {'rule':>6}")
    print(hdr)
    print("-" * len(hdr))

    rows, passing = [], []
    for t in THRESHOLDS:
        labels = fcluster(Z, t=max(t, 1e-9), criterion="distance")
        clusters = {}
        for i, lab in enumerate(labels):
            clusters.setdefault(lab, []).append(i)

        diam = 0.0
        max_size = 0
        for members in clusters.values():
            max_size = max(max_size, len(members))
            if len(members) > 1:
                sub = D[np.ix_(members, members)]
                diam = max(diam, float(sub[sub < SENTINEL].max()))

        merged = sum(1 for (i, j) in dup_idx.values()
                     if labels[i] == labels[j])
        singles = sum(1 for m in clusters.values() if len(m) == 1)

        ok = (merged == len(dup_idx)
              and diam <= MAX_DIAMETER_FACTOR * max(t, 1)
              and max_size <= MAX_CLUSTER_SIZE)
        if ok:
            passing.append(t)

        rows.append((t, len(clusters), len(clusters), singles, diam,
                     max_size, merged, ok))
        print(f"{t:>5} {len(clusters):>9} {len(clusters):>9} {singles:>6} "
              f"{diam:>8.0f} {max_size:>8} "
              f"{merged:>6}/{len(dup_idx):<3} {'PASS' if ok else '':>6}")

    print()
    if passing:
        chosen = min(passing)
        n_ret = next(r[2] for r in rows if r[0] == chosen)
        print("=" * 60)
        print(f"SELECTED THRESHOLD: {chosen}")
        print(f"  retained after dereplication: {n_ret}")
        print(f"  rule: smallest t with all same-strain pairs merged, "
              f"diameter <= {MAX_DIAMETER_FACTOR}t, max cluster <= "
              f"{MAX_CLUSTER_SIZE}")
        print("=" * 60)
        print()
        print("Write into config.yaml:")
        print(f"  dedup:")
        print(f"    max_snp_distance: {chosen}      "
              f"# measured {Path(dist_path).name}")
        print(f"    linkage: complete")
        print(f"  tiers:")
        print(f"    A:")
        print(f"      measured_2026_10_09:")
        print(f"        after_snp_dedup: {n_ret}")
    else:
        print("NO THRESHOLD SATISFIES THE RULE.")
        print("  If dupMerged never reaches the full count, the distance")
        print("  matrix or the accession list is wrong -- do not widen the")
        print("  rule to make it pass.")
        print("  If maxSize blows past 25 before all pairs merge, this")
        print("  dataset has no clean separation between technical noise and")
        print("  real diversity, and dereplication must be reconsidered")
        print("  (e.g. cgMLST allele distance instead of raw SNPs).")
        sys.exit(1)


if __name__ == "__main__":
    main()
