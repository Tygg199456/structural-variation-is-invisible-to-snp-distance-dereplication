#!/usr/bin/env python3
"""
Two-directional positive control for dereplication.

Replaces sanity check [3], whose old form asserted that all five
same-strain assembly pairs must merge. That premise is false: measurement
showed GCF_001293415.1, labelled CO92, is pgm-minus / pPCP1-minus and is
therefore a laboratory derivative, not a second submission of the reference
strain. Requiring it to merge would assert that two different organisms are
the same genome.

So the control now tests BOTH directions:

    four true duplicate pairs  -> must merge
    the CO92 pair              -> must NOT merge, and the deletion profile
                                  must be what prevented it

THE THIRD CONDITION IS WHAT MAKES THIS A TEST OF THE MECHANISM
Asserting only "different cluster" is too weak: it passes trivially if the
pair happens to be distant in SNP space, in which case the deletion
constraint was never exercised at all. The assertion is a conjunction:

    1. raw SNP distance <= threshold    -> distance alone says MERGE
    2. different cluster_id             -> yet they did not merge
    3. deletion profiles differ         -> and only the profile can explain it

Measured: the CO92 pair sits at raw d = 13 against a threshold of 20, so
condition 1 holds and the control is live.

4. optional, strongest: the pair appears in snp_dedup's blocked-pairs list.
   Requires snp_dedup to write blocked pairs involving whitelisted
   accessions (the full Tier B list is ~5.2M pairs; whitelist-only is a few
   dozen and covers exactly the control cases).

A NOTE ON STALE INPUTS
Run this against the CURRENT cluster file. An earlier report of
"CO92 same_cluster=True" came from a sanity run that finished about an hour
BEFORE the deletion-constrained dereplication, so it described the
pre-constraint clustering. Comparing a stale log against current cluster
files is its own failure mode.

Usage:
    python3 scripts/sanity_duplicate_control.py \
        results/dedup/snp_distance_raw_tierA.npy \
        results/dedup/tierA_clusters.tsv \
        results/qc/deletion_profile_all.tsv \
        --threshold 20 \
        --blocked-pairs results/dedup/tierA_blocked_whitelist_pairs.tsv
"""
import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

SENTINEL = 10 ** 6
BIN = 5000

MUST_MERGE = {
    "Antiqua":    ["GCF_000013825.1", "GCF_000834825.1"],
    "PestoidesF": ["GCF_000016445.1", "GCF_000834315.1"],
    "A1122":      ["GCF_000222975.1", "GCF_000834755.1"],
    "Angola":     ["GCF_000018805.1", "GCF_000834845.1"],
}

MUST_SEPARATE = {
    # GCF_001293415.1 is labelled CO92 but is pgm-minus / pPCP1-minus.
    # Verify its true strain identity from BioSample/BioProject and name it
    # correctly in the manuscript -- CO92 is THE reference strain and a
    # mislabelled derivative in the archive is worth a sentence.
    "CO92_vs_derivative": ["GCF_000009065.1", "GCF_001293415.1"],
}


def norm(a):
    return a.split("_", 1)[1] if "_" in a else a


def load_profiles(path):
    """accession -> set of binned large-block keys plus absent replicons."""
    prof = defaultdict(set)
    if not Path(path).exists():
        return prof
    for r in csv.DictReader(open(path), delimiter="\t"):
        acc = r["accession"]
        absent = r.get("replicons_absent", "") or ""
        for rep in (x for x in absent.split(";") if x):
            prof[acc].add(f"ABSENT:{rep}")
        blob = r.get("big_blocks_5kb") or r.get("pgm_sized_coords") or ""
        for tok in (x.strip() for x in blob.replace(",", ";").split(";") if x.strip()):
            prof[acc].add(tok)
        # tolerate a schema that only carries counts
        if not prof[acc]:
            n = r.get("n_blocks_ge10kb")
            if n and int(n) > 0:
                prof[acc].add(f"NBLOCKS:{n}")
    return prof


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dist")
    ap.add_argument("clusters")
    ap.add_argument("profiles")
    ap.add_argument("--threshold", type=float, required=True)
    ap.add_argument("--blocked-pairs",
                    help="optional: snp_dedup's whitelist blocked-pairs list")
    a = ap.parse_args()

    D = np.load(a.dist).astype(float)
    order = [ln.strip() for ln in
             Path(a.dist).with_suffix(".order.txt").read_text().splitlines()
             if ln.strip()]
    idx = {norm(x): i for i, x in enumerate(order)}

    cluster_of = {}
    for r in csv.DictReader(open(a.clusters), delimiter="\t"):
        for m in r["members"].split(","):
            cluster_of[norm(m.strip())] = r["cluster_id"]

    prof = load_profiles(a.profiles)

    blocked = set()
    if a.blocked_pairs and Path(a.blocked_pairs).exists():
        for r in csv.DictReader(open(a.blocked_pairs), delimiter="\t"):
            k = tuple(sorted((norm(r["accession_a"]), norm(r["accession_b"]))))
            blocked.add(k)

    failures, notes = [], []

    print("=" * 80)
    print(f"Dereplication two-directional control   threshold = {a.threshold:g}")
    print("=" * 80)

    # ---- direction 1: true duplicates must merge -------------------------
    print("\n[A] true duplicate pairs -- MUST MERGE")
    for name, accs in MUST_MERGE.items():
        present = [x for x in accs if norm(x) in idx]
        if len(present) < 2:
            notes.append(f"{name}: {len(present)}/2 assemblies present")
            continue
        i, j = idx[norm(present[0])], idx[norm(present[1])]
        d = float(D[i, j])
        same = (cluster_of.get(norm(present[0])) ==
                cluster_of.get(norm(present[1])))
        ok = d <= a.threshold and same
        print(f"    {name:<14} d={d:<7.0f} same_cluster={str(same):<6} "
              f"[{'PASS' if ok else 'FAIL'}]")
        if not ok:
            failures.append(
                f"{name}: d={d:.0f}, merged={same}. Expected d <= "
                f"{a.threshold:g} AND merged. If d is small but merged is "
                f"False, the deletion constraint is over-splitting true "
                f"duplicates -- check their profiles differ for a real reason.")

    # ---- direction 2: the derivative must separate, for the right reason -
    print("\n[B] CO92 vs its derivative -- MUST SEPARATE, BY DELETION PROFILE")
    for name, accs in MUST_SEPARATE.items():
        present = [x for x in accs if norm(x) in idx]
        if len(present) < 2:
            notes.append(f"{name}: {len(present)}/2 assemblies present -- "
                         f"control cannot run")
            continue
        a0, a1 = present
        i, j = idx[norm(a0)], idx[norm(a1)]
        d = float(D[i, j])
        c0, c1 = cluster_of.get(norm(a0)), cluster_of.get(norm(a1))
        p0, p1 = prof.get(a0, set()), prof.get(a1, set())
        exclusive = (p0 ^ p1)

        cond1 = d <= a.threshold                       # distance says merge
        cond2 = (c0 is not None and c1 is not None and c0 != c1)
        cond3 = len(exclusive) >= 1
        cond4 = tuple(sorted((norm(a0), norm(a1)))) in blocked if blocked \
            else None

        print(f"    {a0}  cluster={c0}")
        print(f"    {a1}  cluster={c1}")
        print(f"    1. raw d = {d:.0f} <= {a.threshold:g}          "
              f"[{'PASS' if cond1 else 'FAIL'}]  distance alone says MERGE")
        print(f"    2. different cluster                 "
              f"[{'PASS' if cond2 else 'FAIL'}]  yet they did not merge")
        print(f"    3. deletion profiles differ ({len(exclusive)})     "
              f"[{'PASS' if cond3 else 'FAIL'}]  the only mechanism that fits")
        if cond4 is not None:
            print(f"    4. pair in blocked-pairs list        "
                  f"[{'PASS' if cond4 else 'FAIL'}]  mechanism observed "
                  f"directly")
        if exclusive:
            print(f"       exclusive features: "
                  f"{sorted(exclusive)[:6]}")

        if not cond1:
            failures.append(
                f"{name}: raw d = {d:.0f} exceeds the threshold, so SNP "
                f"distance alone would have separated them and the deletion "
                f"constraint was never exercised. This control is not live; "
                f"pick a pair that is close in SNP space.")
        if not cond2:
            failures.append(
                f"{name}: still in the same cluster ({c0}). The deletion "
                f"constraint did not fire. Check that snp_dedup received "
                f"deletion_profile_all.tsv and that both accessions appear "
                f"in it, and confirm you are reading the CURRENT cluster "
                f"file rather than a pre-constraint run.")
        if not cond3:
            failures.append(
                f"{name}: separated, but the deletion profiles are identical, "
                f"so something else split them. The control does not "
                f"demonstrate the mechanism it claims to test.")
        if cond4 is False:
            failures.append(
                f"{name}: profiles differ and clusters differ, but the pair "
                f"is absent from the blocked-pairs list. Either the list is "
                f"incomplete or the separation came from elsewhere.")

    print()
    for n in notes:
        print(f"NOTE  {n}")
    if failures:
        print()
        for f in failures:
            print(f"FAIL  {f}")
        sys.exit(1)
    print("\nBoth directions pass: duplicates merge, the derivative separates,")
    print("and the deletion profile is demonstrably the reason.")


if __name__ == "__main__":
    main()
