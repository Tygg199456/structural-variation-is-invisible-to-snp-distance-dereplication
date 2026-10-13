#!/usr/bin/env python3
"""
All-pairs SNP distance. v3 emits BOTH a raw and a scaled matrix.

WHY BOTH
--------
v2 emitted only `diff * (n_sites / C)`, which extrapolates observed
differences across positions that were never compared. That is defensible as
a descriptive distance but wrong as the basis for a MERGE decision: merging
two genomes should require observed evidence of near-identity, not
extrapolated evidence.

So:
    snp_distance_raw_<set>.npy   differing positions among JOINTLY CALLED
                                 sites. Interpretable as a SNP count.
                                 -> dereplication thresholds apply to THIS
    snp_distance_<set>.npy       coverage-scaled. Diagnostic only.

Low-overlap pairs are handled by the joint-called guard, not by scaling.

THE NOISE FLOOR THIS MEASURES AGAINST
-------------------------------------
Shredding an assembly into pseudo-reads and remapping cannot call IS-element
copies uniquely. Measured on this dataset: the CO92 reference compared
against itself calls 95.5% of chromosome positions, and independent
assemblies of one strain differ by 3-21 SNPs. Distances below roughly 25 are
within that technical floor, which is why `max_snp_distance: 0` was wrong.
`dedup_threshold_sweep.py` turns that observation into a chosen threshold.
"""
import sys
from pathlib import Path

import numpy as np

snakemake = globals().get("snakemake")  # noqa: F821

BASES = (1, 2, 3, 4)          # 0 = missing
SENTINEL = np.int32(10 ** 6)


def main():
    min_joint = float(getattr(snakemake.params, "min_joint_called_frac", 0.70))

    z = np.load(snakemake.input.npz, allow_pickle=True)
    G = z["genotypes"]
    samples = [s.decode() if isinstance(s, bytes) else str(s)
               for s in z["samples"]]
    n, n_sites = G.shape
    print(f"[snp_distance] matrix {G.shape}, {G.nbytes / 1e6:.1f} MB",
          file=sys.stderr)

    called = (G != 0).astype(np.float32)
    C = called @ called.T

    agree = np.zeros((n, n), dtype=np.float32)
    for b in BASES:
        M = (G == b).astype(np.float32)
        agree += M @ M.T

    diff = C - agree                               # raw, observed

    with np.errstate(divide="ignore", invalid="ignore"):
        scaled = np.where(C > 0, diff * (n_sites / np.maximum(C, 1)), np.inf)

    D_raw = np.rint(diff).astype(np.int32)
    D_scl = np.where(np.isfinite(scaled), np.rint(scaled),
                     SENTINEL).astype(np.int32)

    joint_frac = C / float(n_sites)
    blocked = joint_frac < min_joint
    np.fill_diagonal(blocked, False)
    for M in (D_raw, D_scl):
        M[blocked] = SENTINEL
        np.fill_diagonal(M, 0)

    # Prefer the declared output (see rules/02_dedup.smk); fall back to the
    # derived name so the script still works when run standalone.
    raw_path = Path(getattr(snakemake.output, "raw", None) or str(
        snakemake.output.dist).replace("snp_distance_", "snp_distance_raw_"))
    np.save(raw_path, D_raw)
    np.save(snakemake.output.dist, D_scl)
    for p in (raw_path, Path(snakemake.output.dist)):
        p.with_suffix(".order.txt").write_text("\n".join(samples) + "\n")

    off = ~np.eye(n, dtype=bool)
    usable = D_raw[off & (D_raw != SENTINEL)]
    pcts = [1, 5, 25, 50, 75, 95, 99]
    qs = {p: int(np.percentile(usable, p)) for p in pcts} if usable.size else {}

    # DRIFT WARNING (decided 2026-10-10): min_joint_called_frac stays at 0.70.
    # It has never fired (Tier A and Tier B both blocked = 0; lowest callable
    # fraction seen is 0.902), because the upstream QC already guarantees it.
    # It is redundant insurance, not dead code -- but if a future, more
    # permissive QC lets fragmented drafts through, this is what catches them.
    # Warn while there is still margin, so drift is visible before it bites.
    drift = ""
    if "called_frac" in z:
        lowest = float(z["called_frac"].min())
        if lowest < min_joint + 0.10:
            drift = (f"WARNING: lowest called_frac {lowest:.4f} is within 0.10 "
                     f"of the joint-called guard {min_joint}; the guard is "
                     f"about to start blocking pairs.")

    lines = [
        f"[snp_distance] n={n} sites={n_sites:,} "
        f"min_joint_called_frac={min_joint}",
        f"[snp_distance] called_frac min={z['called_frac'].min():.4f} "
        f"median={np.median(z['called_frac']):.4f} "
        f"max={z['called_frac'].max():.4f}"
        if "called_frac" in z else "[snp_distance] called_frac not stored",
        f"[snp_distance] RAW distance percentiles: {qs}",
        f"[snp_distance] pairs at raw distance 0  : "
        f"{int(((D_raw == 0) & off).sum() // 2)}",
        f"[snp_distance] pairs at raw distance <=25: "
        f"{int(((D_raw <= 25) & off).sum() // 2)}",
        f"[snp_distance] blocked pairs            : "
        f"{int((blocked).sum() // 2)}",
        f"[snp_distance] joint-called guard {min_joint}: never fired on this "
        f"dataset; retained as protection for future permissive QC",
        "",
        "Thresholds belong on the RAW matrix. Run "
        "dedup_threshold_sweep.py before setting max_snp_distance.",
    ]
    if drift:
        lines.insert(0, drift)
        print(drift, file=sys.stderr)
    msg = "\n".join(lines)
    print(msg, file=sys.stderr)
    Path(snakemake.log[0]).write_text(msg + "\n")


if __name__ == "__main__":
    main()
