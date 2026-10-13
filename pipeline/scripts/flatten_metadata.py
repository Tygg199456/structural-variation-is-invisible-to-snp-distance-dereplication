#!/usr/bin/env python3
"""
NCBI Datasets JSONL -> flat TSV + field completeness report.

DECISION THAT MUST NOT BE LOST:
Normalise missing tokens from config (metadata.missing_tokens). Measured
targets to reproduce: Tier A collection_date 89.1%, geo_loc_name 94.2%,
host 60.1%, isolation_source 39.1%; all-set 94.7 / 94.0 / 29.4 / 66.6%.

Field layout is dictated by scripts/assign_tiers.py, which consumes
assembly_metadata.tsv. Any column rename here silently breaks tier assignment.

Inputs / outputs are supplied by Snakemake via the `snakemake` object.
"""

import csv
import json
from pathlib import Path

snakemake = globals().get("snakemake")  # noqa: F821

# Columns assign_tiers.py reads. Order is part of the contract.
COLUMNS = [
    "accession",
    "assembly_level",
    "assembly_status",
    "total_sequence_length",
    "number_of_contigs",
    "contig_n50",
    "gc_percent",
    "ambiguous_fraction",
    "checkm_completeness",
    "checkm_contamination",
    "checkm_version",
    "source_database",
    "paired_accession",
    "reference_genome",
    "strain",
    "collection_date",
    "geo_loc_name",
    "host",
    "isolation_source",
    "bioproject",
    "release_date",
]


def load_jsonl(path):
    """Read NCBI Datasets --as-json-lines, tolerating leading non-JSON lines
    (the CLI can emit progress text; we never want that parsed as a record)."""
    out = []
    p = Path(path)
    if not p.exists():
        return out
    with p.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line.startswith("{"):
                continue
            try:
                obj = json.loads(line)
            except json.JSONDecodeError:
                continue
            if isinstance(obj, dict) and "assembly_info" in obj:
                out.append(obj)
            elif isinstance(obj, dict) and "reports" in obj:
                out.extend(obj["reports"])
    return out


def is_missing(value, missing_tokens):
    if value is None:
        return True
    if isinstance(value, str):
        return value.strip().lower() in missing_tokens
    if isinstance(value, (list, dict)):
        return len(value) == 0
    return False


def biosample_field(rec, field, missing_tokens):
    """Host lives in biosample.attributes; the rest are top-level biosample keys."""
    bs = ((rec.get("assembly_info") or {}).get("biosample")) or {}
    top = {"collection_date", "geo_loc_name", "isolation_source", "strain", "host"}
    value = bs.get(field) if field in top else None
    if is_missing(value, missing_tokens):
        value = None
        for attr in (bs.get("attributes") or []):
            if not isinstance(attr, dict):
                continue
            if str(attr.get("name", "")).lower() == field:
                cand = attr.get("value")
                if not is_missing(cand, missing_tokens):
                    value = cand
                break
    if is_missing(value, missing_tokens):
        return ""
    return str(value).strip()


def ambiguous_fraction(rec):
    """N bases as a fraction of the reported total length.

    NCBI reports atgc_count (unambiguous) alongside total_sequence_length;
    the difference is the ambiguous residue count.
    """
    st = rec.get("assembly_stats") or {}
    try:
        total = int(st["total_sequence_length"])
        atgc = int(st["atgc_count"])
    except (KeyError, TypeError, ValueError):
        return ""
    if total <= 0:
        return ""
    frac = max(0.0, (total - atgc) / total)
    return f"{frac:.6f}"


def row_for(rec, missing_tokens):
    ai = rec.get("assembly_info") or {}
    st = rec.get("assembly_stats") or {}
    ci = rec.get("checkm_info") or {}
    org = rec.get("organism") or {}
    inf = org.get("infraspecific_names") or {}

    strain = inf.get("strain") or inf.get("isolate") or biosample_field(rec, "strain", missing_tokens)

    return {
        "accession": rec.get("accession", ""),
        "assembly_level": ai.get("assembly_level", ""),
        "assembly_status": ai.get("assembly_status", ""),
        "total_sequence_length": st.get("total_sequence_length", ""),
        "number_of_contigs": st.get("number_of_contigs", ""),
        "contig_n50": st.get("contig_n50", ""),
        "gc_percent": st.get("gc_percent", ""),
        "ambiguous_fraction": ambiguous_fraction(rec),
        "checkm_completeness": ci.get("completeness", ""),
        "checkm_contamination": ci.get("contamination", ""),
        "checkm_version": ci.get("checkm_version", ""),
        "source_database": (rec.get("source_database") or "").replace("SOURCE_DATABASE_", ""),
        "paired_accession": rec.get("paired_accession", ""),
        "reference_genome": "1" if ai.get("reference_genome") else "",
        "strain": strain,
        "collection_date": biosample_field(rec, "collection_date", missing_tokens),
        "geo_loc_name": biosample_field(rec, "geo_loc_name", missing_tokens),
        "host": biosample_field(rec, "host", missing_tokens),
        "isolation_source": biosample_field(rec, "isolation_source", missing_tokens),
        "bioproject": ai.get("bioproject_accession", ""),
        "release_date": ai.get("release_date", ""),
    }


FIELDS_FOR_COMPLETENESS = ["collection_date", "geo_loc_name", "host", "isolation_source"]


def completeness_rows(rows, label):
    out = []
    for field in FIELDS_FOR_COMPLETENESS:
        nonnull = sum(1 for r in rows if r[field])
        total = len(rows)
        out.append({
            "dataset": label,
            "field": field,
            "nonnull": nonnull,
            "total": total,
            "fraction": round(nonnull / total, 4) if total else 0.0,
        })
    return out


def main():
    missing_tokens = {t.strip().lower() for t in snakemake.params.missing}

    all_recs = load_jsonl(snakemake.input.allj)
    comp_recs = load_jsonl(snakemake.input.compj)

    rows_all = [row_for(r, missing_tokens) for r in all_recs]
    rows_comp = [row_for(r, missing_tokens) for r in comp_recs]

    # assembly_metadata.tsv feeds assign_tiers.py; it needs every candidate,
    # so it is built from the full taxon pull rather than the complete subset.
    tsv_path = Path(snakemake.output.tsv)
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    with tsv_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows_all)

    comp_rows = completeness_rows(rows_all, "all") + completeness_rows(rows_comp, "complete_chromosome")
    comp_path = Path(snakemake.output.comp)
    comp_path.parent.mkdir(parents=True, exist_ok=True)
    with comp_path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["dataset", "field", "nonnull", "total", "fraction"], delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(comp_rows)

    log_lines = [
        f"records: all={len(rows_all)} complete_chromosome={len(rows_comp)}",
        f"wrote {tsv_path} ({len(rows_all)} rows x {len(COLUMNS)} cols)",
        f"wrote {comp_path}",
        "field completeness:",
    ]
    for r in comp_rows:
        log_lines.append(
            f"  {r['dataset']:20s} {r['field']:18s} {r['nonnull']}/{r['total']} = {r['fraction']:.3f}"
        )
    text = "\n".join(log_lines)
    print(text)
    log = getattr(snakemake, "log", None)
    if log:
        Path(log[0]).write_text(text + "\n")


if __name__ == "__main__":
    main()
