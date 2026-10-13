#!/usr/bin/env python3
"""
Assign assemblies to Tier A / Tier B, or exclude with a recorded reason.

THE CENTRAL DESIGN DECISION, stated here because this is where it is enforced:

    CheckM completeness is NOT used as a filter.

CheckM scores a genome against Enterobacteriaceae single-copy marker genes.
Y. pestis has pseudogenised or deleted a substantial fraction of exactly those
markers -- that is the phenomenon this study measures. On the measured Tier A
candidate set, 110 of 257 genomes (42.8%) fall below 98% completeness, and the
strains that drop out are the deep-branching anchors the analysis depends on
(Angola 96.76, Antiqua 97.70, Pestoides F ~97-98). Filtering on completeness
would delete the signal and call it quality control.

Contamination IS retained: it detects mixed or chimeric assemblies, a real and
independent problem. Completeness is written out as a covariate and tested
against pseudogene count in Fig.S1.

Assembly completeness is instead measured marker-free, as the fraction of the
CO92 chromosome covered in a dnadiff alignment.
"""

import csv
import re
import sys
from pathlib import Path

snakemake = globals().get("snakemake")  # noqa: F821


def parse_dnadiff(path: Path) -> float:
    """AlignedBases percentage for the reference (first) column."""
    pat = re.compile(r"^AlignedBases\s+(\S+)\s+(\S+)")
    for line in path.read_text().splitlines():
        m = pat.match(line.strip())
        if m:
            field = m.group(1)                       # e.g. "4521000(97.15%)"
            inner = field[field.find("(") + 1:field.find("%")]
            try:
                return float(inner)
            except ValueError:
                return 0.0
    return 0.0


def norm_acc(acc: str) -> str:
    """GCF_000009065.1 -> GCF_000009065 (drop version)."""
    return acc.rsplit(".", 1)[0] if acc else ""


def resolve_gca_gcf_pairs(meta, whitelist=frozenset()):
    """Drop the GenBank half of every GenBank/RefSeq pair before download.

    paired_accession is a metadata field, so this costs nothing and saves
    ~1,132 downloads on the real Y. pestis taxon. A GCA whose partner GCF is
    present in the same pull is a mirror of an assembly we already have; the
    RefSeq copy is preferred because it carries NCBI annotation.

    Whitelisted anchors are exempt. KIM10+ is deposited as GCA_000006645.1,
    whose RefSeq mirror is present in the same taxon pull, so an unguarded
    pairing pass silently deletes a phylogenetic anchor -- observed as
    "anchors retained 7/8" on the first Gate 1a run.

    Returns (kept_records, n_dropped).
    """
    gcf_present = {norm_acc(r.get("accession", ""))
                   for r in meta
                   if str(r.get("accession", "")).startswith("GCF_")}
    kept, dropped = [], 0
    for r in meta:
        acc = r.get("accession", "")
        if str(acc).startswith("GCA_") and acc not in whitelist:
            partner = norm_acc(r.get("paired_accession", ""))
            if partner and partner in gcf_present:
                dropped += 1
                continue
        kept.append(r)
    return kept, dropped


def to_float(x, default=None):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def main():
    cfg = snakemake.params.cfg
    # v2: this script runs twice.
    #   stage="metadata"  -- Gate 1a, no downloads, no dnadiff input at all
    #   stage="full"      -- Gate 1b, dnadiff report directory supplied
    # Detected from the presence of the dnadiff input so the rule does not
    # need a second params entry.
    dnadiff_input = getattr(snakemake.input, "dnadiff", None)
    stage = "full" if dnadiff_input else "metadata"
    tA, tB = cfg["tiers"]["A"]["qc"], cfg["tiers"]["B"]["qc"]
    whitelist = {v["accession"] for v in cfg["whitelist"].values()}

    meta_all = list(csv.DictReader(open(snakemake.input.meta), delimiter="\t"))
    if stage == "metadata":
        meta, n_pair_dropped = resolve_gca_gcf_pairs(meta_all, whitelist)
    else:
        meta, n_pair_dropped = meta_all, 0

    dnadiff_dir = Path(dnadiff_input[0]).parent if dnadiff_input else None

    rows_A, rows_B, excluded, report = [], [], [], []

    for r in meta:
        acc = r["accession"]
        level = (r.get("assembly_level") or "").lower()
        total = to_float(r.get("total_sequence_length"), 0) or 0
        contam = to_float(r.get("checkm_contamination"))
        complete = to_float(r.get("checkm_completeness"))
        contigs = to_float(r.get("number_of_contigs"), 1) or 1
        n50 = to_float(r.get("contig_n50"), 0) or 0
        amb = to_float(r.get("ambiguous_fraction"), 0.0) or 0.0

        dd = 0.0
        if dnadiff_dir is not None:
            p = dnadiff_dir / f"{acc}.report"
            if p.exists():
                dd = parse_dnadiff(p)

        is_anchor = acc in whitelist
        reasons: list[str] = []

        # ---- shared gates (contamination + length window + marker-free cover)
        # In the metadata stage there are no dnadiff reports yet, so the
        # marker-free cover gate is silently skipped (dd stays 0.0) and is
        # applied later by refine_tier_a.py on the Tier A candidate list.
        def shared(q):
            bad = []
            if contam is not None and contam > q["max_contamination"]:
                bad.append(f"contamination={contam:.2f}>{q['max_contamination']}")
            if not (q["min_total_bp"] <= total <= q["max_total_bp"]):
                bad.append(f"total_bp={int(total)}_outside_window")
            if amb > q["max_ambiguous_frac"]:
                bad.append(f"ambiguous={amb:.4f}>{q['max_ambiguous_frac']}")
            if stage == "full" and dd and dd < q["min_dnadiff_aligned_pct"]:
                bad.append(f"dnadiff_aligned={dd:.2f}<{q['min_dnadiff_aligned_pct']}")
            return bad

        tier = None
        if level in ("complete genome", "complete", "chromosome"):
            bad = shared(tA)
            if not bad or is_anchor:
                tier = "A"
                reasons = ["anchor_bypass"] if (bad and is_anchor) else []
            else:
                reasons = bad

        if tier is None:
            bad = shared(tB)
            if contigs > tB["max_contigs"]:
                bad.append(f"contigs={int(contigs)}>{tB['max_contigs']}")
            if n50 < tB["min_contig_n50"]:
                bad.append(f"n50={int(n50)}<{tB['min_contig_n50']}")
            if not bad or is_anchor:
                tier = "B"
                reasons = ["anchor_bypass"] if (bad and is_anchor) else []
            else:
                reasons = bad

        rec = {
            "accession": acc,
            "tier": tier or "excluded",
            "assembly_level": level,
            "total_bp": int(total),
            "contigs": int(contigs),
            "contig_n50": int(n50),
            "checkm_contamination": contam,
            # recorded, never filtered on -- see module docstring
            "checkm_completeness_COVARIATE_ONLY": complete,
            "dnadiff_aligned_pct": round(dd, 2),
            "anchor": is_anchor,
            "collection_date": r.get("collection_date", ""),
            "geo_loc_name": r.get("geo_loc_name", ""),
            "host": r.get("host", ""),
            "strain": r.get("strain", ""),
            "notes": ";".join(reasons),
        }
        report.append(rec)

        if tier == "A":
            rows_A.append(rec)
            rows_B.append(rec)          # Tier A is a subset of Tier B
        elif tier == "B":
            rows_B.append(rec)
        else:
            excluded.append(rec)

    def dump(path, rows):
        if not rows:
            Path(path).write_text("accession\n")
            return
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(rows)

    dump(snakemake.output.tierA, rows_A)
    dump(snakemake.output.tierB, rows_B)
    dump(snakemake.output.report, report)
    dump(snakemake.output.excl, excluded)

    # ---- gate-1 console summary -------------------------------------------
    below98 = sum(1 for r in rows_A
                  if (r["checkm_completeness_COVARIATE_ONLY"] or 100) < 98)
    msg = [
        "",
        "=" * 72,
        f"GATE 1 -- tier assignment  (stage={stage})",
        "=" * 72,
        f"metadata records read             : {len(meta_all)}",
        f"GCA/GCF mirrors dropped           : {n_pair_dropped}",
        f"Tier A candidates (pre-SNP-dedup) : {len(rows_A)}",
        f"Tier B candidates                 : {len(rows_B)}",
        f"excluded                          : {len(excluded)}",
        f"anchors retained                  : {sum(1 for r in rows_A if r['anchor'])}/8",
        "",
        f"Tier A with CheckM completeness <98%: {below98}"
        f" ({100 * below98 / max(len(rows_A), 1):.1f}%)",
        "  -> these are RETAINED by design. Filtering them out would remove"
        "     the deep lineages and confound reductive evolution with QC.",
    ]
    if stage == "metadata":
        msg.append("")
        msg.append("  (metadata stage: dnadiff cover gate deferred to refine_tier_a)")
    msg.append("=" * 72)
    msg.append("")
    print("\n".join(msg), file=sys.stderr)
    Path(snakemake.log[0]).write_text("\n".join(msg))


if __name__ == "__main__":
    main()
