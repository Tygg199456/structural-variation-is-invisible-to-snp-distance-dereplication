#!/usr/bin/env python3
"""
Is the archive's structural-variant frequency a natural frequency?

WHY
---
The deletion profile reports that 476/3441 genomes (13.8%) carry the ~105 kb pgm
deletion, 129 lack pPCP1, 66 lack pCD1, 28 lack pMT1. Those numbers cannot be
read as natural frequencies: pgm and plasmid loss happens spontaneously during
laboratory passage, and public archives are full of laboratory, vaccine and
reference strains (EV NIIEG is a vaccine strain, A1122 a CDC diagnostic
reference, KIM10+ a pCD1-deleted derivative). Reading "13.8% of genomes" as
"13.8% of natural isolates" is a serious over-claim and a plague reviewer would
strike it immediately.

So test it instead of hedging. If the structural variants are concentrated in
laboratory-provenance records, the archive frequency is an upper bound and the
bias becomes a quantitative result rather than a qualitative warning.

METHOD
------
Four ORTHOGONAL laboratory-provenance indicators, deliberately not a single
hand-curated strain list (a curated list that has been tuned until the answer
appears is not evidence):

  A  assembly_level == Complete Genome   closed genomes are usually references
  B  strain name matches a curated lab/vaccine/reference pattern
  C  strain name is a bare number        field isolates are numbered, lab strains named
  D  isolation_source empty              field isolates carry a host/source

Each is reported separately, and their enrichment in each structural-variant
group is tested against the no-variant group with a two-sided Fisher exact test.
Agreement across the four is the point; no single indicator is trusted.

Usage:  python3 scripts/archive_bias_check.py
"""
import csv
import re
import sys
from collections import Counter

PROFILE = "results/qc/deletion_profile_all.tsv"
META = "results/qc/assembly_metadata.tsv"
OUT = "results/qc/archive_bias.tsv"

# The pgm locus block, from recurrent_block_polarity.tsv (102.5 kb at 2,135,000).
PGM = ("NC_003143", 2135000, 105000)

# Curated laboratory / vaccine / reference strain designations. Deliberately
# short and mechanical; the point of indicators C and D is to check it.
LAB_PAT = re.compile(
    r"\b(EV\s*76|EV\s*NIIEG|EV76|KIM|KIM10|A1122|CO92|Pasteurella|Pestoides|"
    r"Antiqua|Nepal516|Angola|C-?792|CMCC|SCPM|Tjiwidej|Java9|Harbin35|"
    r"ATCC|NCTC|Biotype|vaccine|attenuat)\b", re.I)


def norm(a):
    return a.strip()


def main():
    # ---- structural variant groups -------------------------------------
    pgm_neg, pPCP1_neg, pCD1_neg, pMT1_neg, clean = set(), set(), set(), set(), set()
    any_block = set()
    for r in csv.DictReader(open(PROFILE), delimiter="\t"):
        acc = norm(r["accession"])
        absent = [x.strip() for x in (r.get("replicons_absent") or "").split(";") if x.strip()]
        blocks = [b for b in (r.get("big_blocks_5kb") or "").split(";") if b.strip()]
        has_pgm = any(
            b.split(":")[0] == PGM[0]
            and abs(int(b.split(":")[1]) - PGM[1]) <= 5000
            and int(b.split(":")[2]) >= 100000
            for b in blocks if len(b.split(":")) == 3
        )
        if has_pgm:
            pgm_neg.add(acc)
        if "NC_003132" in absent:
            pPCP1_neg.add(acc)
        if "NC_003131" in absent:
            pCD1_neg.add(acc)
        if "NC_003134" in absent:
            pMT1_neg.add(acc)
        if not has_pgm and not absent and not blocks:
            clean.add(acc)
        if blocks or absent:
            any_block.add(acc)

    prof_accs = pgm_neg | pPCP1_neg | pCD1_neg | pMT1_neg | clean | any_block
    print("profile accessions touched : %d" % len(prof_accs))
    print("  pgm- (105 kb block)      : %d" % len(pgm_neg))
    print("  pPCP1-                   : %d" % len(pPCP1_neg))
    print("  pCD1-                    : %d" % len(pCD1_neg))
    print("  pMT1-                    : %d" % len(pMT1_neg))
    print("  no block and no replicon loss : %d" % len(clean))

    # ---- metadata ------------------------------------------------------
    meta = {}
    for r in csv.DictReader(open(META), delimiter="\t"):
        meta[norm(r["accession"])] = r
    # mirror accessions map to the same record
    for acc, r in list(meta.items()):
        p = (r.get("paired_accession") or "").strip()
        if p and p not in meta:
            meta[p] = r

    def flags(acc):
        r = meta.get(acc)
        if not r:
            return None
        strain = (r.get("strain") or "").strip()
        return {
            "A_complete": (r.get("assembly_level") or "").strip().lower().startswith("complete"),
            "B_labname": bool(LAB_PAT.search(strain)),
            "C_barenumber": bool(re.fullmatch(r"[0-9]{1,6}[A-Za-z]?", strain)),
            "D_nosource": not (r.get("isolation_source") or "").strip(),
            "strain": strain,
        }

    groups = [
        ("pgm-105kb", pgm_neg),
        ("pPCP1-", pPCP1_neg),
        ("pCD1-", pCD1_neg),
        ("pMT1-", pMT1_neg),
        ("no_variant", clean),
    ]
    base = {a: flags(a) for a in groups[4][1]}
    base = {k: v for k, v in base.items() if v}

    def fisher(a, b, c, d):
        """two-sided Fisher exact, 2x2; exact via hypergeometric sum."""
        from math import comb

        n = a + b + c + d
        if n == 0:
            return 1.0
        r1, c1 = a + b, a + c
        lo, hi = max(0, r1 + c1 - n), min(r1, c1)
        p0 = comb(c1, a) * comb(n - c1, r1 - a) / comb(n, r1)
        p = 0.0
        for x in range(lo, hi + 1):
            px = comb(c1, x) * comb(n - c1, r1 - x) / comb(n, r1)
            if px <= p0 * (1 + 1e-9):
                p += px
        return min(1.0, p)

    print()
    print("%-12s %6s %8s %8s %8s %8s %8s" %
          ("group", "n", "A_comp", "B_lab", "C_num", "D_nosrc", "p(A)"))
    rows = []
    indicators = ["A_complete", "B_labname", "C_barenumber", "D_nosource"]
    base_pos = {k: sum(1 for v in base.values() if v[k]) for k in indicators}
    nb = len(base)
    for name, accs in groups:
        fl = [flags(a) for a in accs]
        fl = [f for f in fl if f]
        if not fl:
            continue
        n = len(fl)
        pos = {k: sum(1 for f in fl if f[k]) for k in indicators}
        pa = fisher(pos["A_complete"], n - pos["A_complete"],
                    base_pos["A_complete"], nb - base_pos["A_complete"])
        print("%-12s %6d %7d%% %7d%% %7d%% %7d%% %8.2e" %
              (name, n, 100 * pos["A_complete"] / n, 100 * pos["B_labname"] / n,
               100 * pos["C_barenumber"] / n, 100 * pos["D_nosource"] / n, pa))
        for k in indicators:
            rows.append({
                "group": name, "n": n, "indicator": k,
                "n_positive": pos[k], "pct_positive": round(100 * pos[k] / n, 2),
                "baseline_n": nb, "baseline_n_positive": base_pos[k],
                "baseline_pct": round(100 * base_pos[k] / nb, 2),
                "fisher_p": "%.3g" % fisher(pos[k], n - pos[k],
                                            base_pos[k], nb - base_pos[k]),
            })

    # how many lab-looking strains in each group, named
    for name, accs in groups[:4]:
        hits = sorted({flags(a)["strain"] for a in accs
                       if flags(a) and flags(a)["B_labname"]})
        print("\n%-12s lab-named strains (%d): %s" %
              (name, len(hits), ", ".join(hits[:25])))

    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print("\nwritten: %s" % OUT)
    print("baseline = no_variant group (no >=10 kb block and no replicon loss), n=%d" % nb)


if __name__ == "__main__":
    sys.exit(main())
