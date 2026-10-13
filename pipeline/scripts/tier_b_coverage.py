#!/usr/bin/env python3
"""
Tier B coverage filter, sourced from snippy rather than dnadiff.

snippy's snps.aligned.fa is reference-length, so the fraction of union sites
that a sample genuinely calls is computed for free by genotype_matrix.py and
stored as `called_frac` in genotypes.npz. Running dnadiff on ~3,200 draft
assemblies to learn the same thing was redundant work in v1.

A fragmented draft that calls few sites is not evidence of a novel lineage, it
is a QC failure. Dropping those before dedup also removes exactly the samples
that would otherwise chain unrelated clusters together through single linkage.
"""

import csv
from pathlib import Path

import numpy as np

snakemake = globals().get("snakemake")  # noqa: F821


def main():
    min_aligned = float(snakemake.params.min_aligned)
    whitelist = {v["accession"] for v in snakemake.params.whitelist.values()}

    z = np.load(snakemake.input.npz, allow_pickle=True)
    samples = [str(s) for s in z["samples"]]
    if "called_frac" not in z:
        raise SystemExit(
            "FATAL: genotypes.npz has no called_frac array. Rebuild it with "
            "the v2 genotype_matrix.py, which reads snps.aligned.fa."
        )
    called = {s: float(c) for s, c in zip(samples, z["called_frac"])}

    rows = list(csv.DictReader(open(snakemake.input.cand), delimiter="\t"))

    kept, dropped, missing = [], [], []
    for r in rows:
        acc = r["accession"]
        # Whitelist is the OUTERMOST filter: every reduction point must exempt
        # anchors, or one of them silently deletes a phylogenetic reference.
        # (KIM10+ was lost this way at the GCA/GCF pairing step.)
        if acc in whitelist:
            r["called_frac"] = round(called.get(acc, float("nan")), 4) if acc in called else ""
            r["coverage_note"] = "anchor_bypass"
            kept.append(r)
            continue
        if acc not in called:
            r["called_frac"] = ""
            r["coverage_note"] = "not_in_genotype_matrix"
            missing.append(r)
            continue
        frac = called[acc]
        r["called_frac"] = round(frac, 4)
        if frac >= min_aligned:
            r["coverage_note"] = ""
            kept.append(r)
        else:
            r["coverage_note"] = f"called_frac={frac:.4f}<{min_aligned}"
            dropped.append(r)

    def dump(path, rs):
        if not rs:
            Path(path).write_text("accession\n")
            return
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rs[0]), delimiter="\t")
            w.writeheader()
            w.writerows(rs)

    dump(snakemake.output.keep, kept)

    from _invariants import assert_anchors_present
    assert_anchors_present([r["accession"] for r in kept],
                           snakemake.params.whitelist, "tier_b_coverage")

    msg = "\n".join([
        f"[tier_b_coverage] candidates={len(rows)} kept={len(kept)} "
        f"dropped={len(dropped)} not_in_matrix={len(missing)} "
        f"anchors_bypassed={sum(1 for r in kept if r.get('coverage_note') == 'anchor_bypass')} "
        f"(min called fraction = {min_aligned})",
    ])
    print(msg)
    Path(snakemake.log[0]).write_text(msg + "\n")


if __name__ == "__main__":
    main()
