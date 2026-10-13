#!/usr/bin/env python3
"""
Expand deletion_profile_all.tsv into the per-block TSV that
recurrent_block_polarity.py expects.

WHY THIS EXISTS
---------------
results/qc/uncalled_blocks.tsv happens to carry the same column layout
(accession / replicon / start / end / length), so pointing the polarity
script at it looks right and runs without complaint. It is the wrong file:
it is produced by uncalled_region_profile.py and covers only the **50
lowest-callability genomes**. Those 50 are enriched for real deletions by
construction (the profile exists precisely to inspect them), so any
recurrence or polarity statistic computed on them is biased and does not
describe the dataset.

Measured, 2026-10-10, same binning as recurrent_block_polarity.py (BIN=5000,
blocks >= 10 kb):

    input                          genomes  distinct coords  recurring >=10x
    uncalled_blocks.tsv                 50               22                5
    deletion_profile_all.tsv          3441              161               18

The 3441 / 161 / 18 row is the one that matches the dataset-level block
inventory (top block NC_003143:795,000 +15 kb absent from 1,127 genomes =
32.8%).

INPUT   results/qc/deletion_profile_all.tsv
        columns: accession, replicons_absent, big_blocks_5kb, n_big, n_absent
        big_blocks_5kb holds "replicon:start:length" entries joined by ";"

OUTPUT  accession, replicon, start, end, length   (one row per block)

Usage:
    python3 scripts/blocks_from_deletion_profile.py \
        results/qc/deletion_profile_all.tsv \
        --out results/qc/uncalled_blocks_all.tsv \
        --min-length 10000
"""
import argparse
import csv
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("profile", help="results/qc/deletion_profile_all.tsv")
    ap.add_argument("--out", default="results/qc/uncalled_blocks_all.tsv")
    ap.add_argument("--min-length", type=int, default=10_000)
    a = ap.parse_args()

    n_rows = n_blocks = n_genomes = 0
    malformed = 0
    genomes = set()
    blocks = []

    with open(a.profile) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            n_rows += 1
            acc = r["accession"]
            genomes.add(acc)
            field = (r.get("big_blocks_5kb") or "").strip()
            if not field:
                continue
            for blk in field.split(";"):
                blk = blk.strip()
                if not blk:
                    continue
                parts = blk.split(":")
                if len(parts) != 3:
                    malformed += 1
                    continue
                rep, start, length = parts[0], int(parts[1]), int(parts[2])
                if length < a.min_length:
                    continue
                blocks.append((acc, rep, start, start + length, length))
                n_blocks += 1

    n_genomes = len(genomes)
    if not blocks:
        sys.exit("FATAL: no blocks at or above the length cutoff -- "
                 "check the column name (big_blocks_5kb)")
    if malformed:
        print(f"WARNING: {malformed} malformed block entries skipped",
              file=sys.stderr)

    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["accession", "replicon", "start", "end", "length"])
        w.writerows(blocks)

    print(f"profile rows      : {n_rows}")
    print(f"genomes seen      : {n_genomes}")
    print(f"blocks >= {a.min_length:,} bp : {n_blocks}")
    print(f"written           : {a.out}")


if __name__ == "__main__":
    main()
