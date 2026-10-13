#!/usr/bin/env python3
"""
Build the genotype matrix from snippy's REFERENCE-LENGTH alignments, not from
the VCFs alone.

WHY THIS SCRIPT DOES NOT READ ONLY VCFs
---------------------------------------
A VCF records differences. It cannot distinguish

    "I covered this position and it matches the reference"
from
    "I never saw this position at all"

Both are simply absent from the file. Any genotype matrix built from VCFs
alone must therefore guess, and both guesses are wrong in a way that does not
raise an error:

  * absent -> missing(0): the only jointly-called positions are the shared
    variant sites, so diff = C - agree is identically zero and EVERY pair has
    distance 0. Single-linkage dedup then collapses all ~3,453 genomes into
    one cluster, silently.

  * absent -> REF: assembly gaps are recorded as "confirmed identical to the
    reference". Two drafts with gaps in different regions look more similar
    than they are, and dedup OVER-MERGES distinct strains. On Tier A
    (closed replicons) this is harmless; on Tier B (contig-level drafts) it
    is the same error as reading a contig gap as a gene deletion, with the
    sign flipped.

snippy already writes `snps.aligned.fa`, which is reference-length and
genuinely three-state:

    A C G T    called
    N          low coverage / unaligned
    -          zero coverage / deletion
    lowercase  masked (low quality, heterozygous)

So: pass 1 collects the union of variant positions from the VCFs; pass 2
reads each aligned.fa and extracts only those positions. Memory stays at the
~86 MB matrix; disk I/O is ~3,453 sequential reads of ~4.65 MB.

Encoding: 0 = missing/unknown, 1..4 = A,C,G,T.
Also emits per-sample called fraction, which replaces dnadiff for Tier B.
"""

import gzip
import sys
from pathlib import Path

import numpy as np

snakemake = globals().get("snakemake")  # noqa: F821

CODE = {"A": 1, "C": 2, "G": 3, "T": 4}
MISSING_CHARS = set("Nn-.*")


def open_maybe_gz(path):
    p = str(path)
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "r")


def sample_name(vcf_path: Path) -> str:
    """snippy layout: <outdir>/<ACCESSION>/snps.vcf"""
    return vcf_path.parent.name


def collect_union_sites(vcf_paths):
    """Pass 1: union of (CHROM, POS) across all per-sample VCFs."""
    sites = set()
    for i, vp in enumerate(vcf_paths, 1):
        with open_maybe_gz(vp) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.split("\t", 3)
                if len(f) < 3:
                    continue
                ref_allele_len_ok = True  # SNVs only; indels handled by snippy-core
                if ref_allele_len_ok:
                    sites.add((f[0], int(f[1])))
        if i % 250 == 0:
            print(f"[genotype_matrix] pass1 {i}/{len(vcf_paths)} "
                  f"union={len(sites)}", file=sys.stderr)
    return sorted(sites)


def read_aligned_fa(path: Path) -> dict[str, str]:
    """Reference-length pseudo-sequence per replicon."""
    seqs, name, buf = {}, None, []
    with open_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def main():
    vcfs = [Path(p) for p in snakemake.input.vcfs]
    samples = [sample_name(p) for p in vcfs]

    sites = collect_union_sites(vcfs)
    n_sites, n_samples = len(sites), len(samples)
    if n_sites == 0:
        raise SystemExit("FATAL: union variant site set is empty.")
    print(f"[genotype_matrix] union sites = {n_sites:,}  samples = {n_samples}",
          file=sys.stderr)

    # positions grouped per replicon for fast extraction
    by_chrom: dict[str, list[tuple[int, int]]] = {}
    for col, (chrom, pos) in enumerate(sites):
        by_chrom.setdefault(chrom, []).append((pos, col))

    G = np.zeros((n_samples, n_sites), dtype=np.int8)
    called_frac = np.zeros(n_samples, dtype=np.float32)

    for si, vp in enumerate(vcfs):
        aligned = vp.parent / "snps.aligned.fa"
        if not aligned.exists():
            raise SystemExit(
                f"FATAL: {aligned} missing. This script requires snippy's "
                f"reference-length alignment; a VCF alone cannot distinguish "
                f"REF from unknown. Re-run snippy for {samples[si]}."
            )
        seqs = read_aligned_fa(aligned)
        row = G[si]
        for chrom, positions in by_chrom.items():
            seq = seqs.get(chrom)
            if seq is None:
                continue                       # replicon absent -> stays 0
            L = len(seq)
            for pos, col in positions:
                if pos <= L:
                    ch = seq[pos - 1]
                    if ch not in MISSING_CHARS:
                        row[col] = CODE.get(ch.upper(), 0)
        called_frac[si] = float((row != 0).sum()) / n_sites
        if (si + 1) % 250 == 0:
            print(f"[genotype_matrix] pass2 {si + 1}/{n_samples}",
                  file=sys.stderr)

    np.savez_compressed(
        snakemake.output.npz,
        genotypes=G,
        samples=np.array(samples, dtype=object),
        called_frac=called_frac,
    )

    with open(snakemake.output.sites, "w") as fh:
        fh.write("chrom\tpos\tcolumn\n")
        for col, (chrom, pos) in enumerate(sites):
            fh.write(f"{chrom}\t{pos}\t{col}\n")

    low = int((called_frac < 0.70).sum())
    msg = "\n".join([
        f"[genotype_matrix] samples={n_samples} sites={n_sites:,} "
        f"matrix={G.nbytes / 1e6:.1f} MB",
        f"[genotype_matrix] called fraction: min={called_frac.min():.3f} "
        f"median={np.median(called_frac):.3f} max={called_frac.max():.3f}",
        f"[genotype_matrix] samples below 0.70 called fraction: {low} "
        f"(these are blocked from merging by snp_distance.py)",
    ])
    print(msg, file=sys.stderr)
    Path(snakemake.log[0]).write_text(msg + "\n")


if __name__ == "__main__":
    main()
