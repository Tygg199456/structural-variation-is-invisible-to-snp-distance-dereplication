#!/usr/bin/env python3
"""Close-out checks requested by the external review (2026-10-11, 2nd round).

1. pgm-only S3 with the UNPRUNED arm2 (all other-block carriers, 154) -- no
   pruning anywhere, no plasmid criterion, no ancestral-state dispute.
2. The 2x2 that closes SS4.2: within A256, chromosome-only (no plasmid sequence
   submitted) vs plasmid-carrying, crossed with dnadiff REF <95% vs >=95%.
3. Identify the two OLD-absent -> NEW-CARRIED genomes (old criterion's false
   positives).
"""
import csv, re, collections
from math import comb


def fisher_two_sided(table):
    (a, b), (c, d) = table
    n = a + b + c + d
    r1, r2, c1, c2 = a + b, c + d, a + c, b + d
    lo, hi = max(0, c1 - r2), min(r1, c1)

    def hyp(x):
        if x < lo or x > hi:
            return 0.0
        return comb(r1, x) * comb(r2, c1 - x) / comb(n, c1)
    obs = hyp(a)
    p = sum(hyp(x) for x in range(lo, hi + 1) if hyp(x) <= obs * (1 + 1e-9))
    aa, bb, cc, dd = a, b, c, d
    if min(aa, bb, cc, dd) == 0:
        aa, bb, cc, dd = aa + .5, bb + .5, cc + .5, dd + .5
    return (aa * dd) / (bb * cc), p


def wilson(k, n, z=1.959964):
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * ((p * (1 - p) / n + z * z / (4 * n * n)) ** .5) / d
    return max(0, c - h), min(1, c + h)


norm = lambda a: re.sub(r"\.\d+$", "", (a or "").strip())
LAB = re.compile(r"\b(EV\s*76|EV\s*NIIEG|EV76|KIM|KIM10|A1122|CO92|Pasteurella|"
                 r"Pestoides|Antiqua|Nepal516|Angola|C-?792|CMCC|SCPM|Tjiwidej|"
                 r"Java9|Harbin35|ATCC|NCTC|Biotype|vaccine|attenuat)\b", re.I)
TAG2P = {"NC_003131": "pCD1", "NC_003132": "pPCP1", "NC_003134": "pMT1"}
PGM = ("NC_003143", 2135000, 105000)
REF_TOTAL = 4829855

prof = {norm(r["accession"]): r for r in
        csv.DictReader(open("results/qc/deletion_profile_all.tsv"), delimiter="\t")}
inv = {norm(r["accession"]): r for r in
       csv.DictReader(open("results/qc/replicon_three_state.tsv"), delimiter="\t")}
meta = {}
for r in csv.DictReader(open("results/qc/assembly_metadata.tsv"), delimiter="\t"):
    meta[norm(r["accession"])] = r
    p = norm(r.get("paired_accession") or "")
    if p:
        meta.setdefault(p, r)
c4 = {norm(r["accession"]): r for r in
      csv.DictReader(open("results/qc/c4_decomposition.tsv"), delimiter="\t")}
pol = {r["polarity"] for r in csv.DictReader(open("results/qc/recurrent_block_polarity.tsv"),
                                             delimiter="\t")}
loss_keys = {(r["replicon"], int(r["start_bin"]), int(r["length_bin"]))
             for r in csv.DictReader(open("results/qc/recurrent_block_polarity.tsv"),
                                     delimiter="\t") if r["polarity"] == "derived_LOSS"}


def blks(a):
    return [b for b in (prof[a].get("big_blocks_5kb") or "").split(";")
            if len(b.split(":")) == 3]


def bk(b):
    p = b.split(":")
    return (p[0], int(p[1]), int(p[2]))


def lab(a):
    return bool(LAB.search(meta.get(a, {}).get("strain") or ""))


def is_named(a):
    ab = [x.strip() for x in (prof[a].get("replicons_absent") or "").split(";") if x.strip()]
    return any(bk(b) == PGM for b in blks(a)) or any(t in ab for t in TAG2P)


P258 = [a for a in prof
        if (meta.get(a, {}).get("assembly_level") or "").lower()
        .startswith(("complete", "chromosome"))]
pgm = [a for a in P258 if any(bk(b) == PGM for b in blks(a))]
other = [a for a in P258 if blks(a) and not is_named(a)]

print("=" * 90)
print("1. pgm-only S3, arm2 UNPRUNED (all other-block carriers)")
print("=" * 90)
k1, n1 = sum(lab(a) for a in pgm), len(pgm)
k2, n2 = sum(lab(a) for a in other), len(other)
orv, p = fisher_two_sided([[k1, n1 - k1], [k2, n2 - k2]])
print("  arm1 pgm        : %d/%d = %.1f%%  Wilson %.1f-%.1f" % (k1, n1, 100 * k1 / n1, *[100 * x for x in wilson(k1, n1)]))
print("  arm2 other(upr) : %d/%d = %.1f%%  Wilson %.1f-%.1f" % (k2, n2, 100 * k2 / n2, *[100 * x for x in wilson(k2, n2)]))
print("  Fisher p = %.4f   OR = %.2f" % (p, orv))

print()
print("=" * 90)
print("2. A256 2x2: chromosome-only vs plasmid-carrying  X  dnadiff REF <95% vs >=95%")
print("=" * 90)
tab = collections.Counter()
low_accs, chrom_accs = [], []
for a in c4:
    un = int(c4[a]["unaligned_report"])
    low = un > 0.05 * REF_TOTAL
    st = inv.get(a)
    chrom = st is not None and st["pCD1_state"] == "UNRESOLVED" \
        and st["pPCP1_state"] == "UNRESOLVED" and st["pMT1_state"] == "UNRESOLVED"
    tab[("chromosome-only" if chrom else "plasmid seqs", "<95" if low else ">=95")] += 1
    if low:
        low_accs.append(a)
    if chrom:
        chrom_accs.append(a)
print("  %-18s %6s %6s" % ("", "<95", ">=95"))
for row in ("chromosome-only", "plasmid seqs"):
    print("  %-18s %6d %6d" % (row, tab[(row, "<95")], tab[(row, ">=95")]))
print("  A256 n=%d  chromosome-only=%d  REF<95 n=%d" % (len(c4), len(chrom_accs), len(low_accs)))
fib = [[tab[("chromosome-only", "<95")], tab[("chromosome-only", ">=95")]],
       [tab[("plasmid seqs", "<95")], tab[("plasmid seqs", ">=95")]]]
orv, p = fisher_two_sided(fib)
print("  Fisher p = %.4f  OR = %.2f" % (p, orv))
print("  REF<95 genomes: %s" % ", ".join("%s(%d)" % (a, int(c4[a]["unaligned_report"])) for a in sorted(low_accs)))

print()
print("=" * 90)
print("3. OLD-absent -> NEW-CARRIED false positives")
print("=" * 90)
for a in sorted(inv):
    if a not in prof:
        continue
    oa = set((prof[a].get("replicons_absent") or "").split(";"))
    for tag, p_ in TAG2P.items():
        if tag in oa and inv[a][p_ + "_state"] == "CARRIED":
            print("  %-18s %s  old=absent  new=CARRIED  n_seq=%s n_pl=%s  ev=%s"
                  % (a, p_, inv[a]["n_seq"], inv[a]["n_plasmid_seq"], inv[a]["evidence"]))
