#!/usr/bin/env python3
"""
How entangled are "gene absent" and "gene pseudogenised" in this pipeline?

THE FINDING THAT FORCES THIS AUDIT
Measured on the CO92 smoke test:

    bakta pseudo=True CDS                        5
    bakta .faa sequences with an internal stop   0
    NCBI RefSeq CDS carrying pseudo=           217  (of 4322)

The zero is structural, not incidental. bakta predicts CDSs with pyrodigal,
and an ORF is by definition start-codon-to-stop-codon with no in-frame stop
inside it. A gene carrying a premature stop therefore CANNOT appear as a
bakta CDS with an internal stop -- it is truncated, split, or not called at
all. So bakta's `pseudo` tag is not a pseudogene detector, and 5 is not a
low detection rate; it is a different quantity. NCBI PGAP runs a dedicated
pseudogene step comparing predicted proteins against references, which is
where 217 comes from. CO92 has closed replicons and no contig ends, so those
217 cannot be edge-truncation artefacts.

THE CONSEQUENCE NOBODY COSTED
If a pseudogenised gene's residual ORF is too short to be called, it is
ABSENT from the bakta GFF3 -- so Panaroo records it as gene LOSS, not as
"present but pseudogenised". Fig.2's absence panel and pseudogene panel are
therefore not independent observations of the same genomes; they are two
recording modes for one phenomenon. Chain 2004 summed both (149 pseudogenes
plus 317 absent genes, 13% inactivated), so combining them has precedent --
but the combination must be declared and the overlap quantified, not left
for the reader to assume two independent lines of evidence.

This script quantifies it, on CO92 only, from two annotations of the SAME
assembly:

    NCBI pseudo region with NO bakta CDS        -> fully entangled; Panaroo
                                                   will score it as absence
    NCBI pseudo region with a SHORTER bakta CDS -> partly entangled
    NCBI pseudo region with a FULL-LENGTH CDS   -> the annotators disagree

The three proportions belong in Methods and in the Fig.2 caption.

It also sets the expectation for pseudofinder: the fully-entangled fraction
is roughly what pseudofinder must recover from sequence comparison, because
those genes are invisible in the bakta CDS set.

SCOPE
PGAP is used here as a one-off BENCHMARK on a single genome, never as
cross-strain measurement data. Mixing PGAP vintages across strains remains
forbidden -- it is the largest source of spurious pseudogene differences.
A benchmark on one curated genome does not reintroduce that problem.

Usage:
    python3 scripts/pseudogene_annotation_audit.py \
        --ncbi-gff data/ref/gff/ncbi_dataset/data/GCF_000009065.1/genomic.gff \
        --bakta-gff /tmp/bakta_smoke/smoke.gff3 \
        --out results/qc/pseudogene_entanglement.tsv

Expected magnitude: ~217 NCBI pseudo regions examined; the three categories
should sum to that. A "full-length CDS" category near 100% would mean the
coordinates are not joining -- check the join-rate warning.
"""
import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

SHORTER_FRAC = 0.80          # bakta CDS below this fraction of the NCBI span
MIN_OVERLAP = 0.50           # of the NCBI span, to count as the same locus


def norm_rep(rep):
    """NCBI writes NC_003143.1; other tables write NC_003143."""
    if "." in rep and rep.rsplit(".", 1)[1].isdigit():
        return rep.rsplit(".", 1)[0]
    return rep


def parse_gff(path, want_types=("CDS",)):
    feats = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in want_types:
                continue
            attrs = {}
            for kv in f[8].split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    attrs[k] = v
            feats[norm_rep(f[0])].append({
                "start": int(f[3]), "end": int(f[4]), "strand": f[6],
                "attrs": attrs,
                "pseudo": ("pseudo" in attrs
                           or attrs.get("pseudogene") is not None),
                "name": (attrs.get("gene") or attrs.get("locus_tag")
                         or attrs.get("ID", "")),
                "product": attrs.get("product", ""),
            })
    for k in feats:
        feats[k].sort(key=lambda x: x["start"])
    return feats


def overlap(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0) + 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ncbi-gff", required=True)
    ap.add_argument("--bakta-gff", required=True)
    ap.add_argument("--out", default="results/qc/pseudogene_entanglement.tsv")
    a = ap.parse_args()

    # ONLY "pseudogene". NCBI emits BOTH a CDS line and a pseudogene line for
    # the same locus, each carrying `pseudo`, so taking both counts every
    # pseudogene twice. Measured 2026-10-11 on CO92: ("CDS","pseudogene") gave
    # 434 pseudo regions against this script's own stated expectation of 217,
    # and the per-category counts were exactly 2x (130/98/206 -> 65/49/103).
    # The PROPORTIONS were unaffected because the doubling was uniform, but the
    # absolute counts -- and the pseudofinder recall target derived from them --
    # were wrong. The sanity check that caught it is the expected-magnitude line
    # in this docstring.
    ncbi = parse_gff(a.ncbi_gff, ("pseudogene",))
    bakta = parse_gff(a.bakta_gff, ("CDS",))

    if not ncbi or not bakta:
        sys.exit("FATAL: one of the GFFs parsed to nothing")

    shared = set(ncbi) & set(bakta)
    if not shared:
        sys.exit(f"FATAL: no shared replicon names.\n"
                 f"  NCBI : {sorted(ncbi)[:4]}\n"
                 f"  bakta: {sorted(bakta)[:4]}\n"
                 f"  Normalisation failed -- the join would be empty, which "
                 f"is a legal result and therefore a silent failure.")

    ncbi_pseudo = [(rep, f) for rep in shared for f in ncbi[rep] if f["pseudo"]]
    bakta_pseudo = [(rep, f) for rep in shared for f in bakta[rep]
                    if f["pseudo"]]

    print("=" * 84)
    print("Pseudogene annotation audit -- CO92, two annotations of one assembly")
    print("=" * 84)
    print(f"shared replicons        : {sorted(shared)}")
    for rep in sorted(shared):
        print(f"  {rep:<16}NCBI CDS {len(ncbi[rep]):>5}   "
              f"bakta CDS {len(bakta[rep]):>5}")
    print(f"NCBI pseudo-flagged     : {len(ncbi_pseudo)}")
    print(f"bakta pseudo-flagged    : {len(bakta_pseudo)}")
    print()

    rows = []
    cat_counts = defaultdict(int)
    for rep, nf in ncbi_pseudo:
        span = nf["end"] - nf["start"] + 1
        best, best_ov = None, 0
        for bf in bakta[rep]:
            if bf["end"] < nf["start"]:
                continue
            if bf["start"] > nf["end"]:
                break
            ov = overlap(nf["start"], nf["end"], bf["start"], bf["end"])
            if ov > best_ov:
                best, best_ov = bf, ov

        if best is None or best_ov < MIN_OVERLAP * span:
            cat = "no_bakta_CDS"
            blen = 0
        else:
            blen = best["end"] - best["start"] + 1
            cat = ("shorter_bakta_CDS" if blen < SHORTER_FRAC * span
                   else "full_length_bakta_CDS")
        cat_counts[cat] += 1
        rows.append({
            "replicon": rep, "ncbi_start": nf["start"], "ncbi_end": nf["end"],
            "ncbi_span": span, "ncbi_gene": nf["name"],
            "ncbi_product": nf["product"][:120],
            "bakta_span": blen,
            "bakta_frac_of_ncbi": round(blen / span, 3) if span else 0,
            "bakta_gene": best["name"] if best else "",
            "category": cat,
        })

    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    if rows:
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(rows)

    total = max(len(rows), 1)
    print("ENTANGLEMENT OF 'ABSENT' AND 'PSEUDOGENISED'")
    print(f"  {'category':<26}{'n':>7}{'pct':>8}  meaning")
    for cat, meaning in (
        ("no_bakta_CDS",
         "fully entangled -- Panaroo will score this as ABSENCE"),
        ("shorter_bakta_CDS",
         "partly entangled -- present in Panaroo, sequence truncated"),
        ("full_length_bakta_CDS",
         "the two annotators disagree about this gene"),
    ):
        n = cat_counts[cat]
        print(f"  {cat:<26}{n:>7}{100 * n / total:>7.1f}%  {meaning}")

    if cat_counts["full_length_bakta_CDS"] > 0.9 * total:
        print()
        print("  WARNING: almost everything landed in 'full_length', which is")
        print("  what a failed coordinate join looks like. Check the replicon")
        print("  names and the overlap threshold before believing this.")

    ent = cat_counts["no_bakta_CDS"] + cat_counts["shorter_bakta_CDS"]
    print()
    print(f"combined entanglement: {ent}/{total} ({100 * ent / total:.1f}%)")
    print()
    print("FOR METHODS")
    print("  Report pseudogenisation and gene absence as ONE inactivation")
    print("  measure, with this overlap stated. Chain 2004 summed them (149")
    print("  pseudogenes plus 317 absent genes, 13% inactivated), so the")
    print("  combination has precedent -- but presenting the two panels as")
    print("  independent evidence would be wrong.")
    print()
    print("FOR PSEUDOFINDER")
    print(f"  {cat_counts['no_bakta_CDS']} genes are invisible in the bakta")
    print("  CDS set, so pseudofinder must recover them from sequence")
    print("  comparison alone. That count is the recall target for the")
    print("  CO92 benchmark.")
    print()
    print("SCOPE: PGAP is used here as a one-off benchmark on ONE genome.")
    print("Using PGAP pseudo calls as cross-strain data stays forbidden.")
    print(f"written: {a.out}")
    print("=" * 84)


if __name__ == "__main__":
    main()
