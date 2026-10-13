#!/usr/bin/env python3
"""
Apply the combined mask BED to a full-length alignment.

THE DECISION THIS SCRIPT ENFORCES
---------------------------------
Masked columns are rewritten as N. They are NOT deleted.

Downstream, `snp-sites -C` counts the constant sites that feed both
`iqtree -fconst` and the BEAST2 `constantSiteWeights`. Deleting masked
columns would silently shrink the constant-site total, and the constant-site
total is exactly what stops BEAST inflating every branch length by ~186x.
A mask step that deletes columns therefore reintroduces the worst failure
mode in the whole pipeline, two rules downstream, with nothing to show it.

snp-sites treats N as missing, so masked positions contribute to neither the
constant nor the variable count. Writing N removes their influence while
keeping the coordinate system intact.

Reference coordinates are 1-based inclusive; BED is 0-based half-open.
"""

import sys
from collections import defaultdict
from pathlib import Path

snakemake = globals().get("snakemake")  # noqa: F821


def read_bed(path):
    regions = defaultdict(list)
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        f = line.split("\t")
        if len(f) < 3:
            continue
        regions[f[0]].append((int(f[1]), int(f[2])))   # 0-based half-open
    merged = {}
    for chrom, ivs in regions.items():
        ivs.sort()
        out = []
        for s, e in ivs:
            if out and s <= out[-1][1]:
                out[-1] = (out[-1][0], max(out[-1][1], e))
            else:
                out.append((s, e))
        merged[chrom] = out
    return merged


def iter_fasta(path):
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:].rstrip("\n"), []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def main():
    regions = read_bed(snakemake.input.bed)

    # snippy-core writes one record per taxon, each the full concatenated
    # reference length, so a single flat coordinate list applies to all of
    # them. If the alignment is per-replicon, the BED chrom names must match.
    flat = []
    for chrom, ivs in regions.items():
        flat.extend(ivs)
    flat.sort()

    total_masked = sum(e - s for s, e in flat)
    n_written = 0
    aln_len = None

    with open(snakemake.output.aln, "w") as out:
        for name, seq in iter_fasta(snakemake.input.aln):
            if aln_len is None:
                aln_len = len(seq)
            elif len(seq) != aln_len:
                raise SystemExit(
                    f"FATAL: {name} has length {len(seq)}, expected {aln_len}. "
                    f"Alignment is not rectangular; masking by coordinate is "
                    f"unsafe."
                )
            chars = list(seq)
            for s, e in flat:
                if s >= aln_len:
                    continue
                e = min(e, aln_len)
                for i in range(s, e):
                    chars[i] = "N"
            out.write(f">{name}\n")
            masked = "".join(chars)
            for i in range(0, len(masked), 60):
                out.write(masked[i:i + 60] + "\n")
            n_written += 1

    if aln_len and total_masked > 0.5 * aln_len:
        print(f"[mask_alignment] WARNING: masking {total_masked:,} of "
              f"{aln_len:,} columns ({100 * total_masked / aln_len:.1f}%). "
              f"Check the repeat/IS BED is not over-broad.", file=sys.stderr)

    msg = (f"[mask_alignment] taxa={n_written}  alignment_length={aln_len:,}  "
           f"masked_columns={total_masked:,} "
           f"({100 * total_masked / max(aln_len or 1, 1):.2f}%)  "
           f"columns written as N, none deleted")
    print(msg, file=sys.stderr)
    Path(snakemake.log[0]).write_text(msg + "\n")


if __name__ == "__main__":
    main()
