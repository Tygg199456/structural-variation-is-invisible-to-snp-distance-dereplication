#!/usr/bin/env python3
"""Final external-review checks (2026-10-11, third round)."""
import csv, re, collections
from math import comb
import numpy as np

RAW = "results/dedup/snp_distance_raw_tierA.npy"
ORD = "results/dedup/snp_distance_raw_tierA.order.txt"
TS = "results/qc/replicon_three_state_seqreport.tsv"
POLAR = "results/qc/recurrent_block_polarity.tsv"
PGM = ("NC_003143", 2135000, 105000)
TAG2P = {"NC_003131": "pCD1", "NC_003132": "pPCP1", "NC_003134": "pMT1"}
LAB = re.compile(r"\b(EV\s*76|EV\s*NIIEG|EV76|KIM|KIM10|A1122|CO92|Pasteurella|"
                 r"Pestoides|Antiqua|Nepal516|Angola|C-?792|CMCC|SCPM|Tjiwidej|"
                 r"Java9|Harbin35|ATCC|NCTC|Biotype|vaccine|attenuat)\b", re.I)
norm = lambda a: re.sub(r"\.\d+$", "", str(a).strip())


def fisher(t):
    (a, b), (c, d) = t
    n = a + b + c + d
    r1, r2, c1, c2 = a + b, c + d, a + c, b + d
    lo, hi = max(0, c1 - r2), min(r1, c1)

    def h(x):
        return comb(r1, x) * comb(r2, c1 - x) / comb(n, c1) if lo <= x <= hi else 0.0
    obs = h(a)
    p = sum(h(x) for x in range(lo, hi + 1) if h(x) <= obs * (1 + 1e-9))
    aa, bb, cc, dd = a, b, c, d
    if min(a, b, c, d) == 0:
        aa, bb, cc, dd = a + .5, b + .5, c + .5, d + .5
    return (aa * dd) / (bb * cc), p


ts = {norm(r["accession"]): (r["pCD1_state"], r["pPCP1_state"], r["pMT1_state"])
      for r in csv.DictReader(open(TS), delimiter="\t")}
prof = {r["accession"]: ((r.get("replicons_absent") or ""), (r.get("big_blocks_5kb") or ""))
        for r in csv.DictReader(open("results/qc/deletion_profile_all.tsv"), delimiter="\t")}
old = {norm(r["accession"]): (r.get("replicons_absent") or "").strip()
       for r in csv.DictReader(open("results/qc/deletion_profile_all.tsv"), delimiter="\t")}
meta = {}
for r in csv.DictReader(open("results/qc/assembly_metadata.tsv"), delimiter="\t"):
    meta[norm(r["accession"])] = r
    p = norm(r.get("paired_accession") or "")
    if p:
        meta.setdefault(p, r)
loss_keys = {(r["replicon"], int(r["start_bin"]), int(r["length_bin"]))
             for r in csv.DictReader(open(POLAR), delimiter="\t")
             if r["polarity"] == "derived_LOSS"}


def blks(a):
    return [b for b in prof[a][1].split(";") if len(b.split(":")) == 3]


def bk(b):
    p = b.split(":")
    return (p[0], int(p[1]), int(p[2]))


def chrom(a):
    return frozenset(b for b in prof[a][1].split(";") if b.startswith("NC_003143:"))


D = np.load(RAW).astype(float)
order = [l.strip() for l in open(ORD) if l.strip()]
idx = {a: i for i, a in enumerate(order)}


def cons_block(ax, ay):
    if chrom(ax) != chrom(ay):
        return True
    ta, tb = ts.get(norm(ax)), ts.get(norm(ay))
    return ta is not None and tb is not None and ta != tb


def classfile(f):
    c = collections.Counter()
    for r in csv.DictReader(open(f), delimiter="\t"):
        if int(r["size"]) <= 1:
            continue
        mem = [m for m in r["members"].split(",") if m]
        mixed = len({prof.get(m) for m in mem}) > 1
        if not mixed:
            c["not_mixed"] += 1
            continue
        if len({chrom(m) for m in mem}) > 1:
            c["a_chromosome"] += 1
            continue
        unres = res = False
        for k in range(3):
            st = [ts.get(norm(m), ("?",) * 3)[k] for m in mem]
            if len(set(st)) > 1:
                if any(s == "UNRESOLVED" for s in st):
                    unres = True
                else:
                    res = True
        c["b_plasmid_unresolved" if unres else "c_plasmid_resolved" if res else "other"] += 1
    return c


for f in ("threestate", "conservative"):
    c = classfile("results/dedup/tierA_clusters.%s.tsv" % f)
    nkeep = sum(1 for _ in open("results/dedup/tierA_keep.%s.tsv" % f)) - 1
    multi = sum(c.values())
    print("(A) %-12s kept=%-4d multi=%-3d mixed=%-3d a=%d b=%d c=%d not_mixed=%d true(a+c)=%d"
          % (f, nkeep, multi, multi - c["not_mixed"], c["a_chromosome"],
             c["b_plasmid_unresolved"], c["c_plasmid_resolved"], c["not_mixed"],
             c["a_chromosome"] + c["c_plasmid_resolved"]))

# conservative d<=20 and CO92
accs = [a for a in order if a in prof]
ii = [idx[a] for a in accs]
Ds = D[np.ix_(ii, ii)]
iu = np.triu_indices(len(accs), 1)
allnear = int((Ds[iu] <= 20).sum())
bl = [(x, y) for x in range(len(accs)) for y in range(x + 1, len(accs))
      if cons_block(accs[x], accs[y])]
nbnear = sum(1 for (x, y) in bl if Ds[x, y] <= 20)
print()
print("(B) conservative: blocked=%d/%d ; blocked with raw d<=20 = %d / %d (%.1f%%)"
      % (len(bl), len(accs) * (len(accs) - 1) // 2, nbnear, allnear, 100 * nbnear / allnear))
print("    (CO92, GCF_001293415.1) blocked (conservative) = %s"
      % cons_block("GCF_000009065.1", "GCF_001293415.1"))

# three-state partition, pgm-only unpruned
P258 = [a for a in prof if norm(a) in ts
        and (meta.get(norm(a), {}).get("assembly_level") or "").lower()
        .startswith(("complete", "chromosome"))]
pgm = [a for a in P258 if any(bk(b) == PGM for b in blks(a))]
lost_any = [a for a in P258 if any(ts[norm(a)][k] == "LOST" for k in range(3))]
named = set(pgm) | set(lost_any)
other = [a for a in P258 if blks(a) and a not in named]


def lab(a):
    return bool(LAB.search(meta.get(norm(a), {}).get("strain") or ""))


k1 = sum(lab(a) for a in pgm)
k2 = sum(lab(a) for a in other)
orv, pv = fisher([[k1, len(pgm) - k1], [k2, len(other) - k2]])
print()
print("(C) three-state pgm-only UNPRUNED: arm1 %d/%d=%.1f%% vs arm2 %d/%d=%.1f%% Fisher p=%.5f OR=%.2f"
      % (k1, len(pgm), 100 * k1 / len(pgm), k2, len(other), 100 * k2 / len(other), pv, orv))

# 5 same-strain pairs
print()
print("(D) 5 same-strain pairs under seq-report three-state:")
for a, b in (("GCF_000018805.1", "GCF_000834845.1"), ("GCF_000016445.1", "GCF_000834315.1"),
             ("GCF_000013825.1", "GCF_000834825.1"), ("GCF_000009065.1", "GCF_001293415.1"),
             ("GCF_000222975.1", "GCF_000834755.1")):
    ta, tb = ts.get(norm(a)), ts.get(norm(b))
    same = ta == tb
    print("    %-18s %s | %-18s %s  -> %s" % (a, ta, b, tb, "CONCORDANT" if same else "DISCORDANT(blocked)"))

# corrected 19-list counts
print()
cnt = collections.Counter()
accs19 = set()
for a in sorted(P258):
    oa = set(old.get(norm(a), "").split(";")) if old.get(norm(a), "") else set()
    if norm(a) not in old:
        continue
    for tag, pp in TAG2P.items():
        if tag not in oa and ts[norm(a)][{"pCD1": 0, "pPCP1": 1, "pMT1": 2}[pp]] == "LOST":
            cnt[pp] += 1
            accs19.add(a)
print("(E) old-present -> new-LOST: rows=%d  pMT1=%d pPCP1=%d pCD1=%d  distinct accessions=%d"
      % (sum(cnt.values()), cnt["pMT1"], cnt["pPCP1"], cnt["pCD1"], len(accs19)))
print("    pMT1 total LOST in P258 = %d ; of which missed by old criterion = %d/13"
      % (sum(1 for a in P258 if ts[norm(a)][2] == "LOST"), cnt["pMT1"]))
