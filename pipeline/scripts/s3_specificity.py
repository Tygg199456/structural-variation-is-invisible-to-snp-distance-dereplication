#!/usr/bin/env python3
"""S3 specificity test, pgm-only arm (external review 2026-10-11).

The published S3 claim: laboratory naming is associated with the four
virulence-related loss classes (pgm, pCD1, pPCP1, pMT1), not with large
structural variants in general.  Because the three plasmid classes are exactly
the ones contaminated by genome-submission practice, this reruns the test with
arm1 restricted to the **pgm** class only -- a chromosome deletion that any
complete/chromosome assembly can resolve and that no plasmid-submission
convention can fake.  If the association survives on pgm alone, the headline
does not depend on any plasmid criterion.

Partition (P258 = complete/chromosome in the deletion profile):
   named-class carriers             arm1 (all four classes)
   other >=10kb block carriers      arm2
   no variant                       (excluded)
arm2_refined = carriers whose every block is a polarity=LOSS block (drops
AMBIGUOUS/GAIN/OUTGROUP_UNINFORMATIVE blocks).

Run with the cobra python (scipy):
  python scripts/s3_specificity.py
"""
import csv, re, sys, collections, os
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _numbers import record

from math import comb


def fisher_exact(table):
    """Two-sided Fisher exact p and sample odds ratio (no scipy needed)."""
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
    if min(aa, bb, cc, dd) == 0:          # Haldane-Anscombe correction
        aa, bb, cc, dd = aa + 0.5, bb + 0.5, cc + 0.5, dd + 0.5
    return (aa * dd) / (bb * cc), p

PROFILE = "results/qc/deletion_profile_all.tsv"
META = "results/qc/assembly_metadata.tsv"
POLAR = "results/qc/recurrent_block_polarity.tsv"
OUT = "results/qc/s3_pgm_only.tsv"

PGM = ("NC_003143", 2135000, 105000)          # 105 kb block at 2,135,000
TAG2P = {"NC_003131": "pCD1", "NC_003132": "pPCP1", "NC_003134": "pMT1"}
LAB_PAT = re.compile(
    r"\b(EV\s*76|EV\s*NIIEG|EV76|KIM|KIM10|A1122|CO92|Pasteurella|Pestoides|"
    r"Antiqua|Nepal516|Angola|C-?792|CMCC|SCPM|Tjiwidej|Java9|Harbin35|"
    r"ATCC|NCTC|Biotype|vaccine|attenuat)\b", re.I)
norm = lambda a: re.sub(r"\.\d+$", "", (a or "").strip())


def wilson(k, n, z=1.959964):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * ((p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5) / d
    return (max(0.0, c - h), min(1.0, c + h))


# ---- polarity -------------------------------------------------------------
loss_keys = set()
polvals = collections.Counter()
for r in csv.DictReader(open(POLAR), delimiter="\t"):
    polvals[r["polarity"]] += 1
    if r["polarity"] == "derived_LOSS":
        loss_keys.add((r["replicon"], int(r["start_bin"]), int(r["length_bin"])))
print("polarity value counts:", dict(polvals))
print("derived_LOSS blocks    : %d" % len(loss_keys))

# ---- profile + metadata ---------------------------------------------------
prof = {}
for r in csv.DictReader(open(PROFILE), delimiter="\t"):
    prof[norm(r["accession"])] = r
meta = {}
for r in csv.DictReader(open(META), delimiter="\t"):
    meta[norm(r["accession"])] = r
    p = norm(r.get("paired_accession") or "")
    if p:
        meta.setdefault(p, r)


def blocks(r):
    return [b for b in (r.get("big_blocks_5kb") or "").split(";") if len(b.split(":")) == 3]


def bkey(b):
    p = b.split(":")
    return (p[0], int(p[1]), int(p[2]))


def absent(r):
    return [x.strip() for x in (r.get("replicons_absent") or "").split(";") if x.strip()]


P258 = [a for a in prof
        if (meta.get(a, {}).get("assembly_level") or "").lower()
        .startswith(("complete", "chromosome"))]
named, other, novar = [], [], []
arm2_refined, pgm_arm = [], []
for a in P258:
    r = prof[a]
    bl = blocks(r)
    ab = absent(r)
    is_pgm = any(bkey(b) == PGM for b in bl)
    is_named = is_pgm or any(t in ab for t in TAG2P)
    if is_named:
        named.append(a)
    elif bl:
        other.append(a)
    else:
        novar.append(a)
    if is_pgm:
        pgm_arm.append(a)
    if bl and not is_named and any(bkey(b) in loss_keys for b in bl):
        arm2_refined.append(a)

print()
print("P258 partition: named=%d other_block=%d no_variant=%d  (sum=%d)"
      % (len(named), len(other), len(novar), len(named) + len(other) + len(novar)))


def labnamed(a):
    st = (meta.get(a, {}).get("strain") or "")
    return bool(LAB_PAT.search(st))


def test(name, arm1):
    k1 = sum(labnamed(a) for a in arm1)
    k2 = sum(labnamed(a) for a in arm2_refined)
    n1, n2 = len(arm1), len(arm2_refined)
    table = [[k1, n1 - k1], [k2, n2 - k2]]
    orv, p = fisher_exact(table)
    lo1, hi1 = wilson(k1, n1)
    lo2, hi2 = wilson(k2, n2)
    print()
    print("%s" % name)
    print("  arm1 : %3d/%3d = %5.1f%%  Wilson %.1f-%.1f"
          % (k1, n1, 100 * k1 / n1 if n1 else 0, 100 * lo1, 100 * hi1))
    print("  arm2r: %3d/%3d = %5.1f%%  Wilson %.1f-%.1f"
          % (k2, n2, 100 * k2 / n2 if n2 else 0, 100 * lo2, 100 * hi2))
    print("  Fisher p = %.5f   OR = %.2f" % (p, orv))
    return (name, k1, n1, k2, n2, p, orv)


r_all = test("S3 REFINED  (arm1 = all four named classes)", named)
r_pgm = test("S3 PGM-ONLY (arm1 = pgm class only)", pgm_arm)

with open(OUT, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["test", "arm1_labnamed", "arm1_n", "arm2_labnamed", "arm2_n",
                "fisher_p", "odds_ratio"])
    for r in (r_all, r_pgm):
        w.writerow([r[0], r[1], r[2], r[3], r[4], "%.6g" % r[5], "%.3f" % r[6]])
    record("s3_specificity.py", {
        "p258_named": len(named), "p258_other_block": len(other),
        "p258_no_variant": len(novar),
        "refined_arm1": [r_all[1], r_all[2]], "refined_arm2": [r_all[3], r_all[4]],
        "refined_p": r_all[5], "refined_or": round(r_all[6], 3),
        "pgm_arm1": [r_pgm[1], r_pgm[2]], "pgm_arm2": [r_pgm[3], r_pgm[4]],
        "pgm_p": r_pgm[5], "pgm_or": round(r_pgm[6], 3)})
print()
print("WROTE " + OUT)
