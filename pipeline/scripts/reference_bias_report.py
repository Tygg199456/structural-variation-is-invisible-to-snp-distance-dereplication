#!/usr/bin/env python3
"""
Separate reference bias from genuine quality problems.

THE PROBLEM THIS SOLVES
-----------------------
Three QC metrics have now flagged the same set of Russian / Central Asian
genomes:

    CheckM completeness   -> penalises GENE LOSS
    dnadiff AlignedBases  -> penalises DIVERGENCE FROM CO92
    called_frac           -> penalises DIVERGENCE FROM CO92

Two mechanisms, one direction: all three measure "how far is this genome
from the modern Orientalis norm". That is the quantity this study exists to
measure. Using any of them as a filter deletes the signal and calls it
quality control, and changing the number does not help -- the tool is wrong,
not the setting.

The right statistic is the RESIDUAL. Regress the metric on divergence from
the reference, then:

    on the trend      -> lineage effect. KEEP, record as a covariate.
    off the trend     -> genuine defect (contamination, chimera, wrong
                         species). THIS is what a QC filter should catch.

Distance to CO92 is already in the raw distance matrix, so the regression
needs no new computation.

Outputs feed the Fig.S5 panel and Table 1's covariate column.

Usage:
    python3 scripts/reference_bias_report.py \
        results/dedup/snp_distance_raw_tierA.npy \
        results/dedup/genotypes_tierA.npz \
        results/qc/tier_report.tsv \
        --reference GCF_000009065.1 \
        --out results/qc/reference_bias.tsv
"""
import argparse
import csv
import sys
from pathlib import Path

import numpy as np

SENTINEL = 10 ** 6
RESIDUAL_Z = 3.0          # flag beyond this many robust SDs below the trend


def norm(a):
    return a.split("_", 1)[1] if "_" in a else a


def robust_z(resid):
    med = np.median(resid)
    mad = np.median(np.abs(resid - med))
    scale = 1.4826 * mad
    if scale <= 0:
        scale = resid.std() or 1.0
    return (resid - med) / scale


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dist")
    ap.add_argument("npz")
    ap.add_argument("meta")
    ap.add_argument("--reference", default="GCF_000009065.1")
    ap.add_argument("--out", default="results/qc/reference_bias.tsv")
    a = ap.parse_args()

    D = np.load(a.dist).astype(float)
    order = [ln.strip() for ln in
             Path(a.dist).with_suffix(".order.txt").read_text().splitlines()
             if ln.strip()]
    idx = {norm(x): i for i, x in enumerate(order)}

    z = np.load(a.npz, allow_pickle=True)
    samples = [s.decode() if isinstance(s, bytes) else str(s)
               for s in z["samples"]]
    called = dict(zip(samples, map(float, z["called_frac"])))

    meta = {r["accession"]: r for r in
            csv.DictReader(open(a.meta), delimiter="\t")}

    ri = idx.get(norm(a.reference))
    if ri is None:
        sys.exit(f"FATAL: reference {a.reference} not in the distance matrix")

    rows = []
    for acc in order:
        i = idx[norm(acc)]
        d = float(D[i, ri])
        if d >= SENTINEL or acc == a.reference:
            continue
        cf = called.get(acc)
        if cf is None:
            continue
        m = meta.get(acc, {})
        rows.append({
            "accession": acc,
            "dist_to_reference": int(d),
            "called_frac": round(cf, 4),
            "dnadiff_aligned_pct": m.get("dnadiff_aligned_pct", ""),
            "checkm_completeness": m.get(
                "checkm_completeness_COVARIATE_ONLY", ""),
            "strain": m.get("strain", ""),
            "geo_loc_name": m.get("geo_loc_name", ""),
        })

    if len(rows) < 10:
        sys.exit("FATAL: too few genomes with both distance and called_frac")

    x = np.array([r["dist_to_reference"] for r in rows], dtype=float)
    y = np.array([r["called_frac"] for r in rows], dtype=float)

    slope, intercept = np.polyfit(x, y, 1)
    pred = slope * x + intercept
    resid = y - pred
    zs = robust_z(resid)

    ss_res = float(((y - pred) ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1 - ss_res / ss_tot if ss_tot else float("nan")

    for r, p, rr, zz in zip(rows, pred, resid, zs):
        r["predicted_called_frac"] = round(float(p), 4)
        r["residual"] = round(float(rr), 4)
        r["residual_z"] = round(float(zz), 2)
        r["flag"] = "QC_OUTLIER" if zz < -RESIDUAL_Z else "on_trend"

    rows.sort(key=lambda r: r["residual_z"])
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    outliers = [r for r in rows if r["flag"] == "QC_OUTLIER"]
    lowest = rows[:12]

    print("=" * 74)
    print(f"Reference-bias regression   reference = {a.reference}")
    print("=" * 74)
    print(f"n genomes            : {len(rows)}")
    print(f"called_frac ~ distance_to_reference")
    print(f"  slope              : {slope:.3e} per SNP")
    print(f"  intercept          : {intercept:.4f}")
    print(f"  R^2                : {r2:.3f}")
    print()
    if r2 >= 0.20:
        print("INTERPRETATION: callability tracks divergence from the")
        print("reference. An absolute called_frac threshold is therefore a")
        print("divergence filter in disguise and must not be used to exclude")
        print("genomes. Record called_frac as a covariate instead.")
    else:
        print("INTERPRETATION: weak relationship. Low callability is NOT")
        print("mainly explained by divergence here, so the low-callability")
        print("genomes deserve individual inspection.")
    print()
    print(f"QC outliers (residual z < -{RESIDUAL_Z}): {len(outliers)}")
    for r in outliers:
        print(f"    {r['accession']}  cf={r['called_frac']} "
              f"pred={r['predicted_called_frac']} z={r['residual_z']}  "
              f"{r['strain']} / {r['geo_loc_name']}")
    if not outliers:
        print("    none -- every low-callability genome sits on the")
        print("    divergence trend, i.e. lineage effect, not a defect.")
    print()
    print("Lowest called_frac, with trend context:")
    print(f"  {'accession':<20}{'cf':>8}{'pred':>8}{'z':>7}{'dist':>8}  "
          f"strain / geo")
    for r in lowest:
        print(f"  {r['accession']:<20}{r['called_frac']:>8.3f}"
              f"{r['predicted_called_frac']:>8.3f}{r['residual_z']:>7.1f}"
              f"{r['dist_to_reference']:>8}  "
              f"{r['strain']} / {r['geo_loc_name']}")
    print()
    print(f"written: {a.out}")
    print("Feeds the Fig.S5 panel and the Table 1 covariate column.")
    print("=" * 74)


if __name__ == "__main__":
    main()
