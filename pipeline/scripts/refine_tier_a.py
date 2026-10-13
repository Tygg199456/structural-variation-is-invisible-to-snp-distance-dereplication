#!/usr/bin/env python3
"""
GATE 1b -- apply the marker-free cover gate to the Tier A candidate list.

Why this is a separate step from assign_tiers.py:
the cover gate needs a dnadiff report per genome, and a report needs the
genome downloaded. Keeping it out of the metadata stage is what turns Gate 1a
into a minutes-long, download-free decision point.

Completeness here is AlignedBases of the CO92 chromosome in a dnadiff
alignment -- not CheckM. CheckM's marker-based completeness is recorded as a
covariate elsewhere and is never a gate on reductive-evolution genomes.

Anchors (config whitelist) bypass the gate: they are the phylogenetic
reference points and must survive every filter.
"""

import csv
import re
import sys
from pathlib import Path

snakemake = globals().get("snakemake")  # noqa: F821


def parse_dnadiff_aligned(path: Path) -> float:
    """AlignedBases percentage for the reference (first) column."""
    if not path.exists():
        return float("nan")
    pat = re.compile(r"^AlignedBases\s+(\S+)\s+(\S+)")
    for line in path.read_text().splitlines():
        m = pat.match(line.strip())
        if m:
            field = m.group(1)
            try:
                return float(field[field.find("(") + 1:field.find("%")])
            except ValueError:
                return float("nan")
    return float("nan")


def main():
    cfg = snakemake.params.cfg
    qc_cfg = cfg.get("qc", {})
    covariate_mode = bool(qc_cfg.get("reference_based_metrics_are_covariates", False))
    threshold = cfg["tiers"]["A"]["qc"]["min_dnadiff_aligned_pct"]
    whitelist = {v["accession"] for v in cfg["whitelist"].values()}

    rows = list(csv.DictReader(open(snakemake.input.cand), delimiter="\t"))
    dnadiff_dir = Path(snakemake.input.dnadiff[0]).parent

    kept, dropped = [], []
    for r in rows:
        acc = r["accession"]
        aligned = parse_dnadiff_aligned(dnadiff_dir / f"{acc}.report")
        r["dnadiff_aligned_pct"] = "" if aligned != aligned else round(aligned, 2)
        below = (aligned == aligned) and aligned < threshold
        if covariate_mode:
            # Reference-based metrics are covariate, never filter: this study
            # measures divergence and decay directly, so any quality score
            # computed against a single reference is collinear with the
            # signal. Measured 2026-10-09: the 5 genomes this gate used to
            # drop are all Russian, and their low AlignedBases reflects real
            # large deletions (whole pMT1 loss) rather than assembly defect.
            r["refine_note"] = "covariate_below_threshold" if below else ""
            kept.append(r)
        elif acc in whitelist:
            r["refine_note"] = "anchor_bypass"
            kept.append(r)
        elif aligned != aligned:
            r["refine_note"] = "dnadiff_report_missing"
            dropped.append(r)
        elif below:
            r["refine_note"] = f"dnadiff_aligned={aligned:.2f}<{threshold}"
            dropped.append(r)
        else:
            r["refine_note"] = ""
            kept.append(r)

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
                           cfg["whitelist"], "refine_tier_a")

    n_below = sum(1 for r in rows
                  if (r.get("dnadiff_aligned_pct") or "") != ""
                  and float(r["dnadiff_aligned_pct"]) < threshold)
    msg = "\n".join([
        f"[refine_tier_a] candidates={len(rows)} kept={len(kept)} "
        f"dropped={len(dropped)} (threshold AlignedBases >= {threshold}%)",
        f"[refine_tier_a] mode={'COVARIATE (no filtering)' if covariate_mode else 'FILTER'}",
        f"[refine_tier_a] {n_below} genomes below threshold, recorded as "
        + ("covariate, retained" if covariate_mode else "excluded"),
        f"[refine_tier_a] anchors kept: "
        f"{sum(1 for r in kept if r.get('anchor') in ('True', True, 'true'))}",
    ])
    print(msg)
    Path(snakemake.log[0]).write_text(msg + "\n")


if __name__ == "__main__":
    main()
