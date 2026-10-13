#!/usr/bin/env python3
"""Next numbers: (1) three-state stratified table (Table 5), (2) conservative
max cluster size, (3) merged CO92_lineage_GAIN region carrier count."""
import csv, re, collections

norm = lambda a: re.sub(r"\.\d+$", "", str(a).strip())
TS = "results/qc/replicon_three_state_seqreport.tsv"
ts = {norm(r["accession"]): (r["pCD1_state"], r["pPCP1_state"], r["pMT1_state"])
      for r in csv.DictReader(open(TS), delimiter="\t")}
prof = {}
for r in csv.DictReader(open("results/qc/deletion_profile_all.tsv"), delimiter="\t"):
    prof[r["accession"]] = ((r.get("replicons_absent") or ""), (r.get("big_blocks_5kb") or ""))
meta = {}
for r in csv.DictReader(open("results/qc/assembly_metadata.tsv"), delimiter="\t"):
    meta[norm(r["accession"])] = r
PGM = ("NC_003143", 2135000, 105000)
GAIN = ("NC_003143", 2365000, 40000)          # merged 2,365,000-2,405,000


def blks(a):
    return [b for b in prof[a][1].split(";") if len(b.split(":")) == 3]


def bk(b):
    p = b.split(":")
    return (p[0], int(p[1]), int(p[2]))


def level(a):
    return (meta.get(norm(a), {}).get("assembly_level") or "").lower()


P258 = [a for a in prof if norm(a) in ts and level(a).startswith(("complete", "chromosome"))]

print("=" * 84)
print("(1) three-state stratified table, P258 (Table 5)")
print("=" * 84)
print("  %-12s %5s %7s %7s %7s %7s" % ("level", "n", "pgm-", "pCD1-", "pPCP1-", "pMT1-"))
for lv in ("complete", "chromosome"):
    S = [a for a in P258 if level(a).startswith(lv)]
    pg = sum(1 for a in S if any(bk(b) == PGM for b in blks(a)))
    c1 = sum(1 for a in S if ts[norm(a)][0] == "LOST")
    c2 = sum(1 for a in S if ts[norm(a)][1] == "LOST")
    c3 = sum(1 for a in S if ts[norm(a)][2] == "LOST")
    print("  %-12s %5d %7d %7d %7d %7d" % (lv, len(S), pg, c1, c2, c3))
S = P258
pg = sum(1 for a in S if any(bk(b) == PGM for b in blks(a)))
print("  %-12s %5d %7d %7d %7d %7d" % ("TOTAL", len(S), pg,
      sum(1 for a in S if ts[norm(a)][0] == "LOST"),
      sum(1 for a in S if ts[norm(a)][1] == "LOST"),
      sum(1 for a in S if ts[norm(a)][2] == "LOST")))
print("  pct of P258: pgm- %.1f%%  pCD1- %.1f%%  pPCP1- %.1f%%  pMT1- %.1f%%"
      % (100 * pg / len(S), 100 * sum(1 for a in S if ts[norm(a)][0] == "LOST") / len(S),
         100 * sum(1 for a in S if ts[norm(a)][1] == "LOST") / len(S),
         100 * sum(1 for a in S if ts[norm(a)][2] == "LOST") / len(S)))

print()
print("=" * 84)
print("(2) conservative max cluster size")
print("=" * 84)
rows = list(csv.DictReader(open("results/dedup/tierA_clusters.conservative.tsv"), delimiter="\t"))
multi = [r for r in rows if int(r["size"]) > 1]
mx = max(rows, key=lambda r: int(r["size"]))
print("  clusters=%d multi=%d  max size=%s (cluster %s, representative %s)"
      % (len(rows), len(multi), mx["size"], mx["cluster_id"], mx["representative"]))
for r in sorted(rows, key=lambda r: -int(r["size"]))[:3]:
    print("    size %-3s cluster %-4s members: %s" % (r["size"], r["cluster_id"], r["members"][:120]))

print()
print("=" * 84)
print("(3) merged CO92_lineage_GAIN region NC_003143:2,365,000-2,405,000")
print("=" * 84)


def overlaps_region(a, reg):
    for b in blks(a):
        k = bk(b)
        if k[0] != reg[0]:
            continue
        if max(k[1], reg[1]) < min(k[1] + k[2], reg[1] + reg[2]):
            return True
    return False


for name, S in (("P3441", list(prof)),
                ("P258", P258)):
    n = sum(1 for a in S if overlaps_region(a, GAIN))
    print("  %-6s n=%d  carrier of merged region = %d (%.1f%%)" % (name, len(S), n, 100 * n / len(S)))
    # also the two raw bins separately
    b1 = ("NC_003143", 2365000, 25000)   # 2,365,000-2,390,000
    b2 = ("NC_003143", 2370000, 35000)   # 2,370,000-2,405,000
    n1 = sum(1 for a in S if overlaps_region(a, b1))
    n2 = sum(1 for a in S if overlaps_region(a, b2))
    print("         raw bins: 2,365,000+25k = %d ; 2,370,000+35k = %d" % (n1, n2))
