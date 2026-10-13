#!/usr/bin/env python3
"""
Why is this genome's called fraction low? Look at WHERE the uncalled
positions are, not at how many there are.

THE QUESTION THIS SETTLES
-------------------------
Three QC metrics flagged the same genomes, and the regression that was
supposed to separate "lineage effect" from "real defect" turned out to be
circular: distance to CO92 is counted only over callable sites, so a
divergent genome loses callability AND loses the SNPs that would have shown
its divergence. Both axes are suppressed by one mechanism.

No continuous divergence measure escapes this. Every between-genome distance
fails in exactly the regions that are most divergent.

So stop looking for a statistical outlier and ask a physical question
instead. Uncalled positions have two distinct causes with two distinct
spatial signatures:

    scattered small blocks (0.7-2.5 kb, dozens to hundreds)
        -> IS elements (IS100/IS285/IS1541/IS1661) defeating unique mapping
        -> TECHNICAL

    a few large contiguous blocks (>= 10 kb, notably one of ~102 kb)
        -> real deletion. The pgm locus (pigmentation segment plus the
           high-pathogenicity island) is ~102 kb and is absent from many
           attenuated and laboratory strains
        -> BIOLOGICAL

    an entire replicon missing
        -> plasmid loss
        -> BIOLOGICAL

That distinction matters more than any threshold: a genome flagged for "low
quality" that is actually carrying a 102 kb deletion is not a quality
problem, it is the exact phenomenon this study measures. Excluding it would
delete the signal.

HYPOTHESIS TO TEST, NOT A CONCLUSION
Several flagged genomes are attenuated or laboratory strains (EV NIIEG is a
Russian live-vaccine lineage; A1122 is a reference/diagnostic strain;
Nepal516 is a sequenced Antiqua representative). Whether each is genuinely
pgm-negative is exactly what this script decides from the data.

Usage:
    python3 scripts/uncalled_region_profile.py \
        results/snippy \
        --accessions GCF_000013805.1 GCF_000222975.1 GCF_000590535.2 \
        --out results/qc/uncalled_profile.tsv
    # or --from-tsv results/qc/reference_bias.tsv --lowest 20
"""
import argparse
import csv
import sys
from pathlib import Path

MISSING = set("Nn-.*")
IS_MIN, IS_MAX = 500, 3000          # IS element size window
LARGE_BLOCK = 10_000                # candidate real deletion
PGM_LIKE_MIN, PGM_LIKE_MAX = 80_000, 130_000   # pgm locus is ~102 kb


def read_fasta(path):
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def runs_of_missing(seq):
    """Yield (start_1based, end_1based, length) for each maximal missing run."""
    start = None
    for i, ch in enumerate(seq):
        if ch in MISSING:
            if start is None:
                start = i
        elif start is not None:
            yield (start + 1, i, i - start)
            start = None
    if start is not None:
        yield (start + 1, len(seq), len(seq) - start)


def classify(acc, aligned_fa):
    per_replicon, blocks = {}, []
    for rep, seq in read_fasta(aligned_fa):
        L = len(seq)
        miss = sum(1 for c in seq if c in MISSING)
        per_replicon[rep] = {"length": L, "missing": miss,
                             "frac_missing": miss / L if L else 1.0}
        for s, e, n in runs_of_missing(seq):
            if n >= IS_MIN:
                blocks.append({"replicon": rep, "start": s, "end": e,
                               "length": n})

    total_len = sum(v["length"] for v in per_replicon.values())
    total_miss = sum(v["missing"] for v in per_replicon.values())

    is_like = [b for b in blocks if IS_MIN <= b["length"] <= IS_MAX]
    large = [b for b in blocks if b["length"] >= LARGE_BLOCK]
    pgm_like = [b for b in large if PGM_LIKE_MIN <= b["length"] <= PGM_LIKE_MAX]
    lost_replicons = [r for r, v in per_replicon.items()
                      if v["frac_missing"] > 0.90]

    bp_large = sum(b["length"] for b in large)
    bp_is = sum(b["length"] for b in is_like)

    if lost_replicons or pgm_like:
        verdict = "BIOLOGICAL_deletion"
    elif bp_large > 0.5 * total_miss:
        verdict = "BIOLOGICAL_large_blocks"
    elif bp_is > 0.3 * total_miss and len(is_like) >= 20:
        verdict = "TECHNICAL_IS_mapping"
    else:
        verdict = "MIXED_inspect"

    return {
        "accession": acc,
        "total_bp": total_len,
        "uncalled_bp": total_miss,
        "called_frac": round(1 - total_miss / total_len, 4) if total_len else 0,
        "n_blocks_ge500bp": len(blocks),
        "n_IS_like_500_3000bp": len(is_like),
        "bp_in_IS_like": bp_is,
        "n_blocks_ge10kb": len(large),
        "bp_in_blocks_ge10kb": bp_large,
        "n_pgm_sized_80_130kb": len(pgm_like),
        "pgm_sized_coords": ";".join(
            f"{b['replicon']}:{b['start']}-{b['end']}" for b in pgm_like),
        "replicons_absent": ";".join(lost_replicons),
        "largest_block_bp": max((b["length"] for b in blocks), default=0),
        "verdict": verdict,
    }, blocks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("snippy_dir")
    ap.add_argument("--accessions", nargs="*", default=[])
    ap.add_argument("--from-tsv", help="reference_bias.tsv or tier_report.tsv")
    ap.add_argument("--lowest", type=int, default=20,
                    help="take the N lowest called_frac rows from --from-tsv")
    ap.add_argument("--out", default="results/qc/uncalled_profile.tsv")
    ap.add_argument("--blocks-out", default="results/qc/uncalled_blocks.tsv")
    a = ap.parse_args()

    accs = list(a.accessions)
    if a.from_tsv:
        rows = [r for r in csv.DictReader(open(a.from_tsv), delimiter="\t")
                if r.get("called_frac")]
        rows.sort(key=lambda r: float(r["called_frac"]))
        accs += [r["accession"] for r in rows[:a.lowest]]
    accs = list(dict.fromkeys(accs))
    if not accs:
        sys.exit("FATAL: no accessions given")

    results, all_blocks = [], []
    for acc in accs:
        fa = Path(a.snippy_dir) / acc / "snps.aligned.fa"
        if not fa.exists():
            print(f"  skip {acc}: {fa} missing", file=sys.stderr)
            continue
        res, blocks = classify(acc, fa)
        results.append(res)
        for b in blocks:
            b["accession"] = acc
            all_blocks.append(b)

    if not results:
        sys.exit("FATAL: no aligned.fa files read")

    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(results[0]), delimiter="\t")
        w.writeheader()
        w.writerows(results)
    with open(a.blocks_out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["accession", "replicon", "start",
                                           "end", "length"], delimiter="\t")
        w.writeheader()
        w.writerows(sorted(all_blocks, key=lambda b: -b["length"]))

    print("=" * 96)
    print("Uncalled-region profile -- spatial cause of low called_frac")
    print("=" * 96)
    print(f"{'accession':<20}{'cf':>7}{'uncalled':>10}{'IS-like':>9}"
          f"{'>=10kb':>8}{'bp>=10kb':>10}{'pgm?':>6}{'largest':>9}  verdict")
    for r in sorted(results, key=lambda r: r["called_frac"]):
        print(f"{r['accession']:<20}{r['called_frac']:>7.3f}"
              f"{r['uncalled_bp']:>10,}{r['n_IS_like_500_3000bp']:>9}"
              f"{r['n_blocks_ge10kb']:>8}{r['bp_in_blocks_ge10kb']:>10,}"
              f"{r['n_pgm_sized_80_130kb']:>6}{r['largest_block_bp']:>9,}"
              f"  {r['verdict']}")

    bio = [r for r in results if r["verdict"].startswith("BIOLOGICAL")]
    tech = [r for r in results if r["verdict"].startswith("TECHNICAL")]
    mixed = [r for r in results if r["verdict"].startswith("MIXED")]

    print()
    print(f"BIOLOGICAL (real deletion / plasmid loss) : {len(bio)}")
    print(f"TECHNICAL  (IS-element mapping failure)   : {len(tech)}")
    print(f"MIXED      (inspect individually)         : {len(mixed)}")
    print()
    if bio:
        print("The BIOLOGICAL group must NOT be excluded on quality grounds.")
        print("A large deletion is the phenomenon this study measures; "
              "removing it")
        print("deletes the signal. Record these as biological observations "
              "and check")
        print("the gene content of the deleted blocks "
              "(results/qc/uncalled_blocks.tsv).")
        for r in bio:
            print(f"    {r['accession']}  largest block "
                  f"{r['largest_block_bp']:,} bp"
                  + (f"  pgm-sized at {r['pgm_sized_coords']}"
                     if r['pgm_sized_coords'] else "")
                  + (f"  replicons absent: {r['replicons_absent']}"
                     if r['replicons_absent'] else ""))
    print()
    print(f"written: {a.out}")
    print(f"         {a.blocks_out}  (every block >=500 bp, for annotation)")
    print("Next: intersect the large blocks with the CO92 annotation to name")
    print("the deleted genes. A ~102 kb block is consistent with the pgm")
    print("locus, but confirm from gene content rather than size alone.")
    print("=" * 96)


if __name__ == "__main__":
    main()
