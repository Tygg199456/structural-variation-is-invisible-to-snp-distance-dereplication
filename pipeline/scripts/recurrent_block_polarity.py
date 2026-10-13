#!/usr/bin/env python3
"""
Annotate the recurrent absent-in-query blocks and determine which genome
changed.

v2 fixes two things in v1.

FIX 1 -- THE OUTGROUP IS A DIFFERENT SPECIES
v1 split on a single outgroup-presence threshold. But IP32953 is
Y. pseudotuberculosis, so its callable fraction against CO92 is far lower
than any Y. pestis genome's. Wherever species-level divergence makes CO92
unmappable, the outgroup reads as "absent" and v1 would call a
CO92-lineage gain that never happened.

v2 reports the outgroup's GLOBAL callable fraction as a baseline and then
uses three bands instead of two:

    outgroup frac >= 0.70   -> derived LOSS in the query
    outgroup frac <= 0.20   -> CO92-lineage GAIN (query is ancestral)
    in between              -> AMBIGUOUS, no call

AMBIGUOUS blocks enter neither the exclusion list nor the Q2 main analysis;
they are reported separately. Forcing a call on them is how a
species-divergence artefact becomes a published polarity claim.

FIX 2 -- USE NCBI'S CO92 ANNOTATION, NOT BAKTA'S
v1's example pointed --gff at the Bakta output directory, which made this
script look blocked on the Bakta install. It never was: --gff is optional,
and the right annotation for naming genes in CO92 coordinates is NCBI's
curated GFF3 for GCF_000009065.1. Bakta's job is uniform re-annotation of
all 206 Tier A genomes so that pseudogene and PAV calls are comparable
across strains -- a different task in a different coordinate space.

WHY POLARITY MATTERS AT ALL
CO92 is Orientalis, a derived lineage carrying prophages, IS insertions and
islands other lineages never had. A CO92-specific 15 kb prophage appears,
under mapping to CO92, as a "15 kb deletion" in 1,127 genomes. Scored as
loss, it becomes the most massively parallel loss in the dataset and the
first row of Table 2 at an overwhelming FDR -- entirely an artefact of
reference choice.

SCOPE
The reference-relative deletion profile is valid for DEREPLICATION whatever
the polarity: a presence/absence difference proves two assemblies are not
the same genome. It must NOT be reused as loss events for ASR or the
parallel-loss null.

Usage:
    python3 scripts/recurrent_block_polarity.py \
        results/qc/uncalled_blocks.tsv \
        --gff data/ref/gff/ncbi_dataset/data/GCF_000009065.1/genomic.gff \
        --outgroup-aligned results/snippy/GCA_000047365.1/snps.aligned.fa \
        --min-recurrence 10 \
        --out results/qc/recurrent_block_polarity.tsv
"""
import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

BIN = 5000
MISSING = set("Nn-.*")


def read_blocks(path, min_len):
    rows = []
    for r in csv.DictReader(open(path), delimiter="\t"):
        n = int(r["length"])
        if n >= min_len:
            rows.append({"accession": r["accession"],
                         "replicon": r["replicon"],
                         "start": int(r["start"]), "length": n})
    return rows


def norm_rep(rep):
    """NCBI GFF uses NC_003143.1; deletion profiles use NC_003143.
    Joining them unnormalised yields an EMPTY gene table, which is a legal
    result and therefore a silent failure."""
    return rep.rsplit(".", 1)[0] if "." in rep and rep.rsplit(".", 1)[1].isdigit() \
        else rep


def key_of(b):
    return (norm_rep(b["replicon"]), b["start"] // BIN * BIN,
            b["length"] // BIN * BIN)


def read_fasta(path):
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = norm_rep(line[1:].split()[0]), []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def outgroup_profile(aligned_fa, regions):
    """Per-region callable fraction, plus the genome-wide baseline."""
    seqs = {r: s for r, s in read_fasta(aligned_fa)}
    total = called_total = 0
    per_replicon = {}
    for rep, seq in seqs.items():
        c = sum(1 for ch in seq if ch not in MISSING)
        per_replicon[rep] = {"length": len(seq), "called": c,
                             "frac": c / len(seq) if seq else 0.0}
        total += len(seq)
        called_total += c
    global_frac = called_total / total if total else 0.0

    per_region = {}
    for (rep, start, length) in regions:
        seq = seqs.get(rep)
        if seq is None:
            per_region[(rep, start, length)] = None
            continue
        s = max(start - 1, 0)
        e = min(s + max(length, BIN), len(seq))
        win = seq[s:e]
        if not win:
            per_region[(rep, start, length)] = None
            continue
        per_region[(rep, start, length)] = \
            sum(1 for ch in win if ch not in MISSING) / len(win)
    return global_frac, per_replicon, per_region


def read_gff_genes(path):
    """NCBI or Bakta GFF3 -> replicon -> [(start, end, name, product)]."""
    genes = defaultdict(list)
    if not path or not Path(path).exists():
        return genes
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("gene", "CDS", "pseudogene"):
                continue
            attrs = {}
            for kv in f[8].split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    attrs[k] = v
            name = (attrs.get("gene") or attrs.get("Name")
                    or attrs.get("locus_tag") or attrs.get("ID", ""))
            genes[norm_rep(f[0])].append((int(f[3]), int(f[4]), name,
                                          attrs.get("product", "")))
    for k in genes:
        genes[k].sort()
    return genes


def genes_in(genes, rep, start, length):
    s, e = start, start + length
    return [(n, p) for (gs, ge, n, p) in genes.get(rep, [])
            if ge >= s and gs <= e]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("blocks")
    ap.add_argument("--gff", help="NCBI CO92 GFF3 (not Bakta's output)")
    ap.add_argument("--outgroup-aligned", nargs="+", default=[],
                    help="one or more outgroup snps.aligned.fa; the consensus "
                         "rule below is ASYMMETRIC")
    ap.add_argument("--min-length", type=int, default=10_000)
    ap.add_argument("--min-recurrence", type=int, default=10)
    ap.add_argument("--loss-frac", type=float, default=0.70,
                    help="outgroup callable fraction at or above this -> LOSS")
    ap.add_argument("--gain-frac", type=float, default=0.20,
                    help="outgroup callable fraction at or below this -> GAIN")
    ap.add_argument("--replicon-present-frac", type=float, default=0.10,
                    help="outgroup callable fraction of a WHOLE replicon below "
                         "this means the outgroup lacks that replicon")
    ap.add_argument("--out",
                    default="results/qc/recurrent_block_polarity.tsv")
    a = ap.parse_args()

    blocks = read_blocks(a.blocks, a.min_length)
    if not blocks:
        sys.exit("FATAL: no blocks at or above the length cutoff")

    counts = Counter(key_of(b) for b in blocks)
    n_genomes = len({b["accession"] for b in blocks})
    recurrent = sorted(((k, c) for k, c in counts.items()
                        if c >= a.min_recurrence), key=lambda kc: -kc[1])
    singletons = sum(1 for c in counts.values() if c == 1)

    # ---- outgroups, one profile each -------------------------------------
    # ASYMMETRIC CONSENSUS. Ancestral PRESENCE needs only one witness; a
    # claim of ancestral ABSENCE needs unanimity, because independent loss is
    # far more common than independent gain. So:
    #     any outgroup retains it   -> the ancestor had it -> derived LOSS
    #     all outgroups lack it     -> CO92-lineage GAIN (weaker evidence)
    #     outgroups disagree        -> AMBIGUOUS, and say so: "two outgroups
    #                                  disagree" is a more informative
    #                                  limitation than "one outgroup was
    #                                  half-callable"
    og_profiles = []
    for path in a.outgroup_aligned:
        if not Path(path).exists():
            print(f"WARNING: {path} missing, skipped", file=sys.stderr)
            continue
        gf, per_rep, per_reg = outgroup_profile(
            path, [k for k, _ in recurrent])
        og_profiles.append({"path": path, "global": gf,
                            "per_replicon": per_rep, "per_region": per_reg})

    if not og_profiles:
        print("WARNING: no outgroup alignment -- every row will read "
              "POLARITY_UNKNOWN.\n"
              "         Run snippy on GCA_000047365.1 "
              "(Y. pseudotuberculosis IP32953) against CO92 first.",
              file=sys.stderr)
        global_frac, per_replicon, og = 0.0, {}, {}
    else:
        global_frac = og_profiles[0]["global"]
        per_replicon = og_profiles[0]["per_replicon"]
        og = og_profiles[0]["per_region"]

    genes = read_gff_genes(a.gff)
    if not genes:
        print("NOTE: no GFF supplied or none parsed -- polarity still works, "
              "gene names will be empty. Use NCBI's GFF3 for "
              "GCF_000009065.1.", file=sys.stderr)

    rows = []
    for (rep, start, length), count in recurrent:
        frac = og.get((rep, start, length))
        fracs = [p["per_region"].get((rep, start, length))
                 for p in og_profiles]
        fracs = [f for f in fracs if f is not None]
        rep_fracs = [(p["per_replicon"].get(rep) or {}).get("frac")
                     for p in og_profiles]
        rep_fracs = [f for f in rep_fracs if f is not None]

        frac = max(fracs) if fracs else None          # the best witness
        rep_frac = max(rep_fracs) if rep_fracs else None
        n_loss_votes = sum(1 for f in fracs if f >= a.loss_frac)
        n_gain_votes = sum(1 for f in fracs if f <= a.gain_frac)
        n_informative = len(fracs)

        if rep_frac is not None and rep_frac < a.replicon_present_frac:
            # The outgroup does not carry this replicon at all. pMT1 and pPCP1
            # are Y. pestis-specific, acquired at speciation, so the outgroup
            # is STRUCTURALLY uninformative here -- calling absence a
            # "CO92-lineage gain" would be a category error. Polarity for
            # these must come from within-species ASR using basal lineages
            # (Pestoides / 0.PE) as the ancestral proxy.
            polarity = "OUTGROUP_UNINFORMATIVE"
            note = (f"outgroup carries only {rep_frac:.1%} of this replicon; "
                    f"it is absent from the outgroup, so no polarity can be "
                    f"read from it")
        elif frac is None:
            polarity = "POLARITY_UNKNOWN"
            note = "no outgroup data for this replicon"
        elif n_loss_votes >= 1:
            polarity = "derived_LOSS"
            note = (f"{n_loss_votes}/{n_informative} outgroup(s) retain the "
                    f"region (best {frac:.2f}); one witness to ancestral "
                    f"presence is enough, so absence in the query is a loss")
        elif n_gain_votes == n_informative and n_informative >= 1:
            polarity = "CO92_lineage_GAIN"
            note = (f"all {n_informative} outgroup(s) lack it -> CO92 "
                    f"acquired it and the query is ancestral, NOT a loss"
                    + ("" if n_informative > 1 else
                       "; single outgroup, so this is the weaker direction"))
        else:
            polarity = "AMBIGUOUS"
            if n_informative > 1:
                note = (f"outgroups disagree "
                        f"({', '.join(f'{f:.2f}' for f in fracs)}); a real "
                        f"interspecies structural difference, not a single "
                        f"half-callable outgroup")
            else:
                note = (f"outgroup partially callable ({frac:.2f}); "
                        f"species-level divergence cannot be separated from "
                        f"true absence")

        g = genes_in(genes, rep, start, length)
        rows.append({
            "replicon": rep, "start_bin": start, "length_bin": length,
            "n_genomes_absent": count,
            "pct_genomes_absent": round(100 * count / n_genomes, 1),
            "outgroup_called_frac": "" if frac is None else round(frac, 3),
            "outgroup_fracs_all": ";".join(f"{f:.3f}" for f in fracs),
            "n_outgroups_informative": n_informative,
            "n_loss_votes": n_loss_votes,
            "n_gain_votes": n_gain_votes,
            "polarity": polarity,
            "n_genes": len(g),
            "gene_names": ",".join(n for n, _ in g if n)[:400],
            "products": " | ".join(sorted({p for _, p in g if p}))[:400],
            "note": note,
        })

    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    if rows:
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(rows)

    print("=" * 88)
    print("Recurrent absent-in-query blocks, with event polarity")
    print("=" * 88)
    print(f"genomes with any block >= {a.min_length:,} bp : {n_genomes}")
    print(f"distinct binned coordinates                : {len(counts)}")
    print(f"  recurring >= {a.min_recurrence}x                          : "
          f"{len(recurrent)}")
    print(f"  singletons (lineage/strain specific)     : {singletons}")
    print()
    print("OUTGROUP BASELINE -- read this before trusting any polarity call")
    for p in og_profiles:
        print(f"  {Path(p['path']).parent.name}  global {p['global']:.3f}")
        for rep, v in sorted(p["per_replicon"].items()):
            print(f"    {rep:<16}{v['frac']:.3f}  "
                  f"({v['called']:,}/{v['length']:,} bp)")
    if len(og_profiles) == 1:
        print("  NOTE: one outgroup only. GAIN calls rest on a single "
              "witness; add a second (IP 31758 = GCF_000016945.1) to "
              "strengthen them and to resolve AMBIGUOUS blocks.")
    if global_frac and global_frac < 0.60:
        print("  WARNING: the outgroup calls under 60% of CO92. At this level")
        print("  'absent in outgroup' is dominated by species divergence, and")
        print("  GAIN calls are unreliable. Consider a closer outgroup (a")
        print("  second Y. pseudotuberculosis) or report polarity only for")
        print("  blocks where the outgroup is confidently callable.")
    print()

    print(f"{'replicon':<16}{'start':>11}{'len':>9}{'n':>7}{'pct':>7}"
          f"{'og':>7}  polarity")
    for r in rows:
        print(f"{r['replicon']:<16}{r['start_bin']:>11,}{r['length_bin']:>9,}"
              f"{r['n_genomes_absent']:>7}{r['pct_genomes_absent']:>6.1f}%"
              f"{str(r['outgroup_called_frac']):>7}  {r['polarity']}")
        if r["gene_names"]:
            print(f"{'':<16}genes: {r['gene_names'][:150]}")

    # Join-rate assertion. An empty gene table is a legal result, so without
    # this the replicon-naming mismatch that produced n_genes=0 everywhere
    # fails silently.
    if genes:
        annotated = sum(1 for r in rows if r["n_genes"] > 0)
        rate = annotated / max(len(rows), 1)
        print()
        print(f"gene-annotation join rate: {annotated}/{len(rows)} "
              f"({rate:.0%})")
        if rate < 0.80:
            print("  WARNING: most blocks got no genes. Check that the GFF "
                  "seqids match the block replicon names -- NCBI writes "
                  "NC_003143.1 while the deletion profiles use NC_003143.")
            print(f"  GFF replicons : {sorted(genes)[:4]}")
            print(f"  block replicons: "
                  f"{sorted({r['replicon'] for r in rows})[:4]}")

    gains = [r for r in rows if r["polarity"] == "CO92_lineage_GAIN"]
    losses = [r for r in rows if r["polarity"] == "derived_LOSS"]
    amb = [r for r in rows if r["polarity"] == "AMBIGUOUS"]
    unk = [r for r in rows if r["polarity"] == "POLARITY_UNKNOWN"]
    uninf = [r for r in rows if r["polarity"] == "OUTGROUP_UNINFORMATIVE"]

    print()
    print(f"derived losses     : {len(losses)}")
    print(f"CO92-lineage gains : {len(gains)}")
    print(f"ambiguous          : {len(amb)}")
    print(f"polarity unknown   : {len(unk)}")
    print(f"outgroup uninformative (replicon absent from outgroup): "
          f"{len(uninf)}")
    if uninf:
        print("  These sit on Y. pestis-specific replicons (pMT1, pPCP1).")
        print("  Polarise them by within-species ASR with basal lineages as")
        print("  the ancestral proxy, not by the outgroup.")

    if gains or amb:
        print()
        print("Exclude BOTH the gains and the ambiguous blocks from the")
        print("parallel-loss test. Gains are not losses; ambiguous blocks")
        print("cannot be shown to be either, and a forced call here becomes a")
        print("published polarity claim resting on species divergence.")
        print()
        print("config.yaml:")
        print("  parallel_loss:")
        if gains:
            print("    exclude_regions:          # determined: CO92-lineage "
                  "gains, not losses")
            for r in gains:
                print(f"      - {{replicon: {r['replicon']}, "
                      f"start: {r['start_bin']}, length: {r['length_bin']}, "
                      f"reason: CO92_lineage_gain}}")
        if amb:
            print("    exclude_regions_ambiguous:  # undetermined: not "
                  "polarisable with this outgroup")
            for r in amb:
                print(f"      - {{replicon: {r['replicon']}, "
                      f"start: {r['start_bin']}, length: {r['length_bin']}, "
                      f"reason: ambiguous_polarity}}")
    if uninf:
        # NOT excluded. Plasmid loss is real biological signal for Q2 --
        # pMT1 carries ymt and caf1, so losing it means losing flea
        # transmissibility. These blocks need polarity from within-species
        # ASR, not removal from the analysis.
        print("    needs_internal_polarity:   # plasmid blocks; polarise by "
              "within-species ASR, do NOT exclude")
        for r in uninf:
            print(f"      - {{replicon: {r['replicon']}, "
                  f"start: {r['start_bin']}, length: {r['length_bin']}, "
                  f"reason: outgroup_lacks_replicon}}")

    print()
    print("SCOPE: the reference-relative deletion profile stays valid for")
    print("DEREPLICATION regardless of polarity. It must not be reused as")
    print("loss events for ASR or the parallel-loss null.")
    print(f"written: {a.out}")
    print("=" * 88)


if __name__ == "__main__":
    main()
