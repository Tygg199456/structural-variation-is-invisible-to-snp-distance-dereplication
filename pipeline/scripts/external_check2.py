#!/usr/bin/env python3
"""Second round of external-review checks (2026-10-11): five requested numbers.

(1) classify the threestate run's multi-member clusters a/b/c
(2) confirm CO92 / GCF_001293415.1 remain blocked under threestate
(3) rebuild the P258 partition under the three-state and rerun pgm-only Fisher
(4) list the 19 old-present -> new-LOST (accession x replicon)
(5) blocked pairs among raw d <= 20, for full and threestate
"""
import csv, re, os, collections
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


def fisher(table):
    (a, b), (c, d) = table
    n = a + b + c + d
    r1, r2, c1, c2 = a + b, c + d, a + c, b + d
    lo, hi = max(0, c1 - r2), min(r1, c1)

    def h(x):
        return comb(r1, x) * comb(r2, c1 - x) / comb(n, c1) if lo <= x <= hi else 0.0
    obs = h(a)
    p = sum(h(x) for x in range(lo, hi + 1) if h(x) <= obs * (1 + 1e-9))
    if min(table[0][0], table[0][1], table[1][0], table[1][1]) == 0:
        a, b, c, d = a + .5, b + .5, c + .5, d + .5
    return (a * d) / (b * c), p


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
    return [b for b in (prof[a][1]).split(";") if len(b.split(":")) == 3]


def bk(b):
    p = b.split(":")
    return (p[0], int(p[1]), int(p[2]))


# ---- (1) classify a run's multi-member clusters ---------------------------
def classify(cluster_file):
    rows = [r for r in csv.DictReader(open(cluster_file), delimiter="\t")
            if int(r["size"]) > 1]
    cnt = collections.Counter()
    details = []
    for r in rows:
        mem = [m for m in r["members"].split(",") if m]
        mixed = len({prof.get(m) for m in mem}) > 1
        if not mixed:
            cat = "not_mixed"
        else:
            chrom = [frozenset(b for b in prof[m][1].split(";") if b.startswith("NC_003143:"))
                     for m in mem]
            if len(set(chrom)) > 1:
                cat = "a_chromosome"
            else:
                unres = res = False
                for k in range(3):
                    st = [ts.get(norm(m), ("?",) * 3)[k] for m in mem]
                    if len(set(st)) > 1:
                        (unres := True) if any(s == "UNRESOLVED" for s in st) else (res := True)
                cat = ("b_plasmid_unresolved" if unres else
                       "c_plasmid_resolved" if res else "other")
        cnt[cat] += 1
        details.append((r["cluster_id"], r["size"], cat))
    return cnt, details


for f in ("off", "full", "threestate"):
    c, d = classify("results/dedup/tierA_clusters.%s.tsv" % f)
    nm = sum(v for k, v in c.items() if k != "not_mixed")
    print("(1) %-11s multi=%d mixed=%d  a=%d b=%d c=%d not_mixed=%d  true(a+c)=%d"
          % (f, sum(c.values()), nm, c["a_chromosome"], c["b_plasmid_unresolved"],
             c["c_plasmid_resolved"], c["not_mixed"],
             c["a_chromosome"] + c["c_plasmid_resolved"]))

# ---- (2) CO92 pair still blocked under threestate --------------------------
print()
print("(2) threestate blocked pairs involving CO92 (recomputed directly):")


D = np.load(RAW).astype(float)
order = [l.strip() for l in open(ORD) if l.strip()]
idx = {a: i for i, a in enumerate(order)}


def ts_block(ax, ay):
    cx = frozenset(b for b in prof[ax][1].split(";") if b.startswith("NC_003143:"))
    cy = frozenset(b for b in prof[ay][1].split(";") if b.startswith("NC_003143:"))
    if cx != cy:
        return True
    tx, ty = ts.get(norm(ax)), ts.get(norm(ay))
    if tx and ty:
        return any(tx[k] != "UNRESOLVED" and ty[k] != "UNRESOLVED" and tx[k] != ty[k]
                   for k in range(3))
    return False


co92 = "GCF_000009065.1"
cand = [a for a in order if a in prof]
n_co = 0
pair = None
for a in cand:
    if a == co92:
        continue
    if co92 in (a,) or True:
        if ts_block(co92, a):
            n_co += 1
        if a == "GCF_001293415.1":
            pair = ts_block(co92, a)
print("   CO92 blocked pairs (threestate) = %d" % n_co)
print("   (CO92, GCF_001293415.1) blocked under threestate = %s" % pair)

# ---- (3) P258 partition under three-state + pgm-only Fisher ----------------
print()
P258 = [a for a in prof if norm(a) in ts
        and (meta.get(norm(a), {}).get("assembly_level") or "").lower()
        .startswith(("complete", "chromosome"))]
pgm = [a for a in P258 if any(bk(b) == PGM for b in blks(a))]
lost_any = [a for a in P258 if any(ts[norm(a)][k] == "LOST" for k in range(3))]
named = set(pgm) | set(lost_any)
other = [a for a in P258 if blks(a) and a not in named]
novar = [a for a in P258 if not blks(a) and a not in named]
arm2r = [a for a in other if any(bk(b) in loss_keys for b in blks(a))]


def lab(a):
    return bool(LAB.search(meta.get(norm(a), {}).get("strain") or ""))


k1 = sum(lab(a) for a in pgm)
k2 = sum(lab(a) for a in arm2r)
orv, pv = fisher([[k1, len(pgm) - k1], [k2, len(arm2r) - k2]])
print("(3) P258 three-state partition: pgm=%d lost_any=%d named=%d other=%d novar=%d (sum %d)"
      % (len(pgm), len(lost_any), len(named), len(other), len(novar),
         len(named) + len(other) + len(novar)))
print("    arm2_refined(other, >=1 derived_LOSS block) = %d" % len(arm2r))
print("    pgm-only: arm1 %d/%d=%.1f%%  arm2r %d/%d=%.1f%%  Fisher p=%.5f OR=%.2f"
      % (k1, len(pgm), 100 * k1 / len(pgm), k2, len(arm2r), 100 * k2 / len(arm2r), pv, orv))

# ---- (4) 19 old-present -> new-LOST ---------------------------------------
print()
print("(4) OLD_present -> NEW_LOST (P258):")
n = 0
for a in sorted(P258):
    oa = set(old.get(norm(a), "").split(";")) if old.get(norm(a), "") else set()
    if norm(a) not in old:
        continue
    for tag, pp in TAG2P.items():
        if tag not in oa and ts[norm(a)][{"pCD1": 0, "pPCP1": 1, "pMT1": 2}[pp]] == "LOST":
            n += 1
            print("   %-18s %s" % (a, pp))
print("   total = %d" % n)

# ---- (5) blocked pairs among raw d <= 20 ----------------------------------
print()


def blocked_sets():
    accs = [a for a in order if a in prof]
    out = {}
    fb = set()
    for x in range(len(accs)):
        for y in range(x + 1, len(accs)):
            if prof[accs[x]] != prof[accs[y]]:
                fb.add((x, y))
    out["full"] = (accs, fb)
    tb = set()
    for x in range(len(accs)):
        for y in range(x + 1, len(accs)):
            if ts_block(accs[x], accs[y]):
                tb.add((x, y))
    out["threestate"] = (accs, tb)
    return out


for mode, (accs, bs) in blocked_sets().items():
    ii = [idx[a] for a in accs]
    Ds = D[np.ix_(ii, ii)]
    iu = np.triu_indices(len(accs), 1)
    dvals = Ds[iu]
    allnear = int((dvals <= 20).sum())
    nbnear = sum(1 for (x, y) in bs if Ds[x, y] <= 20)
    print("(5) %-11s blocked=%d / %d pairs ; blocked with raw d<=20 = %d ; all pairs d<=20 = %d"
          % (mode, len(bs), len(accs) * (len(accs) - 1) // 2, nbnear, allnear))
