#!/usr/bin/env python3
"""
Find genuine quality outliers by comparing each genome to ITS OWN lineage,
not to a regression on a divergence proxy.

WHY THE REGRESSION WAS REPLACED
-------------------------------
`called_frac ~ distance_to_CO92` is circular. Distance is counted only over
callable sites, so a divergent genome loses callability and simultaneously
loses the SNPs that would have revealed its divergence. Both axes are
suppressed by one mechanism, and a truly divergent strain lands at
"short distance, low callability" -- exactly where the model calls it an
outlier. Measured result: R^2 = 0.017 and 44/255 genomes flagged, two of
them anchors.

The scaled matrix does not rescue it either: scaled = raw / (cf_i * cf_j)
divides cf out explicitly, so regressing cf on it is cf ~ 1/cf and inflates
R^2 instead. Raw gives a lower bound, scaled an upper bound; neither is the
answer.

No continuous divergence measure escapes this -- cgMLST, k-mer and
multi-reference distances all fail in precisely the regions that are most
divergent.

So the question changes. Not "is callability explained by divergence" but
"after controlling for lineage, is anything left". Lineage becomes a FACTOR,
not a continuous proxy:

    1. coarse hierarchical clustering of the raw distance matrix
       (~8-15 groups, approximating populations)
    2. robust z of called_frac WITHIN each group
    3. below 3 MADs of its OWN lineage -> genuine defect
       low along with its whole lineage -> lineage effect, keep

Residual caveat, to state in Methods: the clustering uses a distance matrix
that is itself deflated by callability. Distances span 25-1617 while the
deflation is ~20%, so coarse group membership is not changed by it, but the
caveat belongs in the text.

Run `uncalled_region_profile.py` FIRST. Genomes whose uncalled positions sit
in large contiguous blocks carry real deletions and are biological
observations, not QC cases; they should be excluded from this test rather
than flagged by it.

Usage:
    python3 scripts/lineage_stratified_qc.py \
        results/dedup/snp_distance_raw_tierA.npy \
        results/dedup/genotypes_tierA.npz \
        results/qc/tier_report.tsv \
        --n-groups 12 \
        --deletion-profile results/qc/uncalled_profile.tsv \
        --out results/qc/lineage_stratified_qc.tsv
"""
import argparse
import csv
import sys
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

SENTINEL = 10 ** 6
Z_CUT = -3.0


def norm(a):
    return a.split("_", 1)[1] if "_" in a else a


def robust_z(v):
    med = np.median(v)
    mad = np.median(np.abs(v - med))
    scale = 1.4826 * mad
    if scale <= 0:
        scale = v.std() or 1.0
    return (v - med) / scale


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dist")
    ap.add_argument("npz")
    ap.add_argument("meta")
    ap.add_argument("--n-groups", type=int, default=12)
    ap.add_argument("--deletion-profile",
                    help="uncalled_profile.tsv; BIOLOGICAL rows are excused")
    ap.add_argument("--out", default="results/qc/lineage_stratified_qc.tsv")
    a = ap.parse_args()

    D = np.load(a.dist).astype(float)
    order = [ln.strip() for ln in
             Path(a.dist).with_suffix(".order.txt").read_text().splitlines()
             if ln.strip()]
    D = np.minimum(D, SENTINEL)
    np.fill_diagonal(D, 0)
    D = (D + D.T) / 2.0

    z = np.load(a.npz, allow_pickle=True)
    samples = [s.decode() if isinstance(s, bytes) else str(s)
               for s in z["samples"]]
    cf = dict(zip(samples, map(float, z["called_frac"])))
    meta = {r["accession"]: r for r in
            csv.DictReader(open(a.meta), delimiter="\t")}

    excused = {}
    if a.deletion_profile and Path(a.deletion_profile).exists():
        for r in csv.DictReader(open(a.deletion_profile), delimiter="\t"):
            if r.get("verdict", "").startswith("BIOLOGICAL"):
                excused[r["accession"]] = r["verdict"]

    Z = linkage(squareform(D, checks=False), method="average")
    labels = fcluster(Z, t=a.n_groups, criterion="maxclust")

    groups = {}
    for i, lab in enumerate(labels):
        groups.setdefault(int(lab), []).append(i)

    rows = []
    for gid, members in sorted(groups.items()):
        vals, idxs = [], []
        for i in members:
            v = cf.get(order[i])
            if v is not None:
                vals.append(v)
                idxs.append(i)
        if len(vals) < 4:
            for i, v in zip(idxs, vals):
                rows.append(_row(order[i], gid, len(members), v, np.nan,
                                 np.nan, "group_too_small", meta, excused))
            continue
        arr = np.array(vals)
        zs = robust_z(arr)
        med = float(np.median(arr))
        for i, v, zz in zip(idxs, arr, zs):
            acc = order[i]
            if acc in excused:
                flag = "biological_deletion"
            elif zz < Z_CUT:
                flag = "QC_OUTLIER"
            else:
                flag = "on_lineage_trend"
            rows.append(_row(acc, gid, len(members), float(v), med,
                             float(zz), flag, meta, excused))

    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: (r["group"],
                                                r["called_frac"])))

    print("=" * 84)
    print(f"Lineage-stratified callability QC   groups = {a.n_groups}")
    print("=" * 84)
    print(f"{'group':>6}{'n':>5}{'cf_median':>11}{'cf_min':>9}"
          f"{'outliers':>10}  representative geo")
    for gid, members in sorted(groups.items()):
        sub = [r for r in rows if r["group"] == gid]
        if not sub:
            continue
        cfs = [r["called_frac"] for r in sub]
        geos = [r["geo_loc_name"] for r in sub if r["geo_loc_name"]]
        top = max(set(geos), key=geos.count) if geos else "-"
        n_out = sum(1 for r in sub if r["flag"] == "QC_OUTLIER")
        print(f"{gid:>6}{len(sub):>5}{np.median(cfs):>11.3f}"
              f"{min(cfs):>9.3f}{n_out:>10}  {top}")

    outliers = [r for r in rows if r["flag"] == "QC_OUTLIER"]
    bio = [r for r in rows if r["flag"] == "biological_deletion"]
    print()
    print(f"genuine QC outliers (within-lineage z < {Z_CUT}) : {len(outliers)}")
    for r in outliers:
        print(f"    {r['accession']}  cf={r['called_frac']:.3f} "
              f"lineage_median={r['lineage_median_cf']:.3f} "
              f"z={r['within_lineage_z']:.1f}  "
              f"{r['strain']} / {r['geo_loc_name']}")
    print(f"excused as real deletions                     : {len(bio)}")
    print()
    print("Genomes on their own lineage's trend are a lineage effect and "
          "must be")
    print("kept, with called_frac recorded as a covariate.")
    print(f"written: {a.out}")
    print("=" * 84)


def _row(acc, gid, gsize, v, med, zz, flag, meta, excused):
    m = meta.get(acc, {})
    return {"accession": acc, "group": gid, "group_size": gsize,
            "called_frac": round(v, 4),
            "lineage_median_cf": round(med, 4) if med == med else "",
            "within_lineage_z": round(zz, 2) if zz == zz else "",
            "flag": flag,
            "deletion_verdict": excused.get(acc, ""),
            "strain": m.get("strain", ""),
            "geo_loc_name": m.get("geo_loc_name", "")}


if __name__ == "__main__":
    main()
