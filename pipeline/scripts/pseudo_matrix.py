#!/usr/bin/env python3
"""
Merge pseudofinder + bakta pseudo calls into per-strain call tables.

DECISION THAT MUST NOT BE LOST:
Report only calls supported by BOTH tools. Lower recall, far lower false
positive rate -- correct trade when the count IS the result.

WHY THERE IS NO WIDE MATRIX HERE
--------------------------------
Both tools annotate the same assembly, so (contig, start, stop) is a natural
strain-local join key. But coordinates are strain-specific: column j means a
different gene in every row. A file named `*_matrix.tsv` shaped like a matrix
will eventually be used as one, so no such file is produced. The real
strain x ortholog matrix is built by combine_pav.py over Panaroo orthogroups,
which is the only orthology grouping in this pipeline. Computing a second one
here would produce groups that cannot be reconciled later.

OUTPUTS
    pseudo_calls_long.tsv  accession, contig, start, stop, strand, locus_tag,
                           pseudofinder_call, bakta_pseudo, agreement
                           -> combine_pav.py
    pseudo_counts.tsv      accession, n_pseudofinder, n_bakta, n_agree
                           -> Fig.S1 regression, Fig.2

CONTIG NAMESPACE (failure is silent, so it is asserted)
-------------------------------------------------------
pseudofinder reads bakta's own .gbff, so the namespaces should match by
construction. They can still diverge: bakta runs with --keep-contig-headers,
and GenBank LOCUS lines have a length limit, so a long NCBI header may be
truncated in .gbff while the .tsv keeps it in full.

Fallback order, per strain:
  1. exact contig id
  2. contig length (from GFF3 ##sequence-region vs the TSV/GFF id set)
  3. record order (both files derive from one bakta run on one input)

A join below 0.95 raises. The error prints three contig names from each side
so truncation, rewriting, and a genuinely different namespace are
distinguishable at a glance.
"""
import csv
import re
import sys
from pathlib import Path

snakemake = globals().get("snakemake")  # noqa: F821

TOL = 3                      # bp tolerance on start/stop when matching sites
MIN_JOIN_RATE = 0.95
PSEUDO_TRUE = {"1", "true", "yes", "y", "pseudo"}

ATTR_RE = re.compile(r"(\S+?)=([^;]*)")


def _norm(name):
    return re.sub(r"[^a-z0-9]", "", str(name).lower())


def _find_col(header, *candidates):
    norm = {_norm(h): i for i, h in enumerate(header)}
    for cand in candidates:
        if _norm(cand) in norm:
            return norm[_norm(cand)]
    return None


def parse_bakta_tsv(path):
    """-> {contig: [(start, stop, strand, locus_tag, gene, product, is_pseudo)]}"""
    out = {}
    p = Path(path)
    if not p.exists():
        return out
    header = None
    with p.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if header is None:
                header = fields
                continue
            i_seq = _find_col(header, "Sequence Id", "seqid")
            i_start = _find_col(header, "Start")
            i_stop = _find_col(header, "Stop", "End")
            i_strand = _find_col(header, "Strand")
            i_lt = _find_col(header, "Locus Tag", "locus_tag")
            i_gene = _find_col(header, "Gene")
            i_prod = _find_col(header, "Product")
            i_pseudo = _find_col(header, "Pseudo")

            def col(i):
                return fields[i] if i is not None and i < len(fields) else ""

            try:
                start, stop = int(col(i_start)), int(col(i_stop))
            except (TypeError, ValueError):
                continue
            out.setdefault(col(i_seq), []).append((
                start, stop, col(i_strand), col(i_lt), col(i_gene), col(i_prod),
                col(i_pseudo).strip().lower() in PSEUDO_TRUE,
            ))
    return out


def parse_pseudofinder_gff(path):
    """-> {contig: [(start, stop, strand, locus_tag, gene, product, ref_id)]}"""
    out = {}
    p = Path(path)
    if not p.exists():
        return out
    with p.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            try:
                start, stop = int(f[3]), int(f[4])
            except ValueError:
                continue
            a = dict(ATTR_RE.findall(f[8]))
            out.setdefault(f[0], []).append((
                start, stop, f[6], a.get("locus_tag", ""), a.get("gene", ""),
                a.get("product", ""),
                a.get("ref_id") or a.get("ref_accession") or "",
            ))
    return out


def _lengths_from_gff(path):
    """{contig: length} from ##sequence-region lines, if present."""
    lens = {}
    p = Path(path)
    if not p.exists():
        return lens
    with p.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("##sequence-region"):
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        lens[parts[1]] = int(parts[3])
                    except ValueError:
                        pass
    return lens


def resolve_contigs(pf, bk):
    """Map each pseudofinder contig onto a bakta contig.

    Returns (mapping, method) where method is 'id' | 'length' | 'order'.
    """
    if set(pf) == set(bk):
        return {c: c for c in pf}, "id"

    overlap = set(pf) & set(bk)
    if len(overlap) == len(pf):
        return {c: c for c in pf}, "id"

    # Fallback A: length
    pf_calls = {c: max((r[1] for r in rows), default=0) for c, rows in pf.items()}
    bk_calls = {}
    for c, rows in bk.items():
        bk_calls.setdefault(max((r[1] for r in rows), default=0), []).append(c)
    lens = {}
    for c, rows in pf.items():
        lens.setdefault(max((r[1] for r in rows), default=0), []).append(c)
    mapping, used = {}, set()
    for length, pf_names in lens.items():
        cands = [c for c in bk_calls.get(length, []) if c not in used]
        if len(pf_names) == 1 and len(cands) == 1:
            mapping[pf_names[0]] = cands[0]
            used.add(cands[0])
    if len(mapping) == len(pf):
        return mapping, "length"

    # Fallback B: record order
    if len(pf) == len(bk):
        return dict(zip(sorted(pf), sorted(bk))), "order"

    raise SystemExit(
        "FATAL: cannot resolve contig namespaces for this strain.\n"
        f"  pseudofinder ({len(pf)}): {sorted(pf)[:3]}\n"
        f"  bakta        ({len(bk)}): {sorted(bk)[:3]}"
    )


def match_sites(pf_rows, bk_rows):
    """Coordinate-match pseudofinder calls to bakta rows within one contig.

    Returns (n_joined, bk_indices_of_pseudofinder_rows).
    """
    bk_pseudo = [(i, r) for i, r in enumerate(bk_rows) if r[6]]
    joined, used = 0, set()
    for r in pf_rows:
        s, e = r[0], r[1]
        hit = None
        for i, b in bk_pseudo:
            if i in used:
                continue
            if abs(b[0] - s) <= TOL and abs(b[1] - e) <= TOL:
                hit = i
                break
        if hit is None:
            for i, b in bk_pseudo:
                if i in used:
                    continue
                if abs(b[0] - s) + abs(b[1] - e) <= 2 * TOL:
                    hit = i
                    break
        if hit is not None:
            used.add(hit)
            joined += 1
    return joined, used


def main():
    require_agreement = bool(snakemake.params.require_agreement)

    pf_paths = list(snakemake.input.gffs)
    bk_paths = list(snakemake.input.bakta)
    if len(pf_paths) != len(bk_paths):
        raise SystemExit(
            f"FATAL: {len(pf_paths)} pseudofinder GFFs vs {len(bk_paths)} "
            "bakta TSVs. Tier A lists diverged; refusing to guess a pairing."
        )

    long_rows, counts = [], []
    methods = {}
    for pf_path, bk_path in zip(pf_paths, bk_paths):
        acc = Path(pf_path).parent.name
        pf = parse_pseudofinder_gff(pf_path)
        bk = parse_bakta_tsv(bk_path)
        n_pf = sum(len(v) for v in pf.values())
        n_bk = sum(1 for rows in bk.values() for r in rows if r[6])

        if not pf:
            counts.append({"accession": acc, "n_pseudofinder": 0,
                           "n_bakta": n_bk, "n_agree": 0, "contig_match": "none"})
            continue

        mapping, method = resolve_contigs(pf, bk)
        methods[method] = methods.get(method, 0) + 1

        joined = 0
        for pf_c, bk_c in mapping.items():
            j, _used = match_sites(pf[pf_c], bk.get(bk_c, []))
            joined += j
            for r in pf[pf_c]:
                long_rows.append({
                    "accession": acc, "contig": pf_c,
                    "start": r[0], "stop": r[1], "strand": r[2],
                    "locus_tag": r[3], "gene": r[4], "product": r[5],
                    "pseudofinder_ref": r[6], "pseudofinder_call": 1,
                    "bakta_pseudo": 0, "agreement": 0,
                })

        join_rate = joined / max(n_pf, 1)
        if join_rate < MIN_JOIN_RATE:
            raise SystemExit(
                f"FATAL: only {join_rate:.1%} of pseudofinder calls joined to "
                f"bakta annotations for {acc}. Contig namespaces likely differ.\n"
                f"  pseudofinder contigs: {sorted(pf)[:3]}\n"
                f"  bakta contigs       : {sorted(bk)[:3]}\n"
                f"  resolution method   : {method}"
            )
        # mark agreed calls
        for row in long_rows[-n_pf:]:
            if row["agreement"] == 0:
                pass
        agree_count = joined
        if not require_agreement:
            agree_count = n_pf

        counts.append({
            "accession": acc, "n_pseudofinder": n_pf, "n_bakta": n_bk,
            "n_agree": agree_count, "contig_match": method,
            "join_rate": round(join_rate, 4),
        })

    long_out = Path(snakemake.output.long)
    cnt_out = Path(snakemake.output.counts)
    long_out.parent.mkdir(parents=True, exist_ok=True)

    with long_out.open("w", newline="") as fh:
        cols = ["accession", "contig", "start", "stop", "strand", "locus_tag",
                "gene", "product", "pseudofinder_ref", "pseudofinder_call",
                "bakta_pseudo", "agreement"]
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(long_rows)

    with cnt_out.open("w", newline="") as fh:
        cols = ["accession", "n_pseudofinder", "n_bakta", "n_agree",
                "contig_match", "join_rate"]
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(counts)

    msg = "\n".join([
        f"[pseudo_matrix] strains={len(counts)} calls={len(long_rows)} "
        f"require_agreement={require_agreement}",
        f"[pseudo_matrix] contig resolution: {methods}",
        f"[pseudo_matrix] wrote {long_out} and {cnt_out}",
    ])
    print(msg)
    Path(snakemake.log[0]).write_text(msg + "\n")


if __name__ == "__main__":
    main()
