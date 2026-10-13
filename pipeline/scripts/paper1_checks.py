#!/usr/bin/env python3
"""
The four analyses the external review put at P0/P1, in one pass.

P0-1  What are the genomes that fall OUTSIDE both arms of the pooled test?
      The pooled test compared 72 variant carriers against 32 no-variant
      genomes, but the complete/chromosome set is much larger. If the missing
      genomes carry >=10 kb blocks outside the four named classes, then the
      test is mislabelled: it ran "four named classes vs none", not
      "any variant vs none".

P0-2  Cluster 349's within-cluster distances and its temporal signal.
      Comparing the cluster diameter (bounded at t by construction, since
      complete linkage guarantees it) against a 5-pair noise floor was
      circular. Median/p95 within-cluster distance are NOT bounded by t, and a
      Mantel-type permutation on year differences answers whether the cluster
      carries time signal at all.

P1-3  The five dnadiff-excluded strains: within-group vs nearest-outside
      distance, to decide whether "a single clonal group" may be written.

P1-4  Structural-variant frequencies after excluding lab-named genomes, i.e.
      how much of the archive frequency is laboratory.

Run: python3 scripts/paper1_checks.py
"""
import csv
import itertools
import os
import re
from collections import Counter, defaultdict

import numpy as np

LAB = re.compile(r"\b(EV\s*76|EV\s*NIIEG|EV76|KIM|KIM10|A1122|CO92|Pasteurella|"
                 r"Pestoides|Antiqua|Nepal516|Angola|C-?792|CMCC|SCPM|Tjiwidej|"
                 r"Java9|Harbin35|ATCC|NCTC|vaccine|attenuat)\b", re.I)

REPL = {"NC_003143": "chromosome", "NC_003134": "pMT1",
        "NC_003131": "pCD1", "NC_003132": "pPCP1"}
NAMED_CLASS = {"pgm": ("NC_003143", 2135000),
               "pPCP1": ("NC_003132", None),
               "pCD1": ("NC_003131", None),
               "pMT1": ("NC_003134", None)}
SENT = 10 ** 6


def load():
    prof = {}
    for r in csv.DictReader(open("results/qc/deletion_profile_all.tsv"), delimiter="\t"):
        prof[r["accession"].strip()] = r
    meta = {}
    for r in csv.DictReader(open("results/qc/assembly_metadata.tsv"), delimiter="\t"):
        meta[r["accession"].strip()] = r
        p = (r.get("paired_accession") or "").strip()
        if p:
            meta.setdefault(p, r)
    return prof, meta


def blocks_of(r):
    return [b for b in (r.get("big_blocks_5kb") or "").split(";")
            if len(b.split(":")) == 3]


def is_pgm(b):
    p = b.split(":")
    return (p[0] == "NC_003143" and abs(int(p[1]) - 2135000) <= 5000
            and int(p[2]) >= 100000)


def complete_level(prof, meta):
    return [a for a in prof
            if (meta.get(a, {}).get("assembly_level") or "").lower()
            .startswith(("complete", "chromosome"))]


def main():
    prof, meta = load()
    comp = complete_level(prof, meta)
    print("=" * 96)
    print("P0-1  WHY THE TWO ARMS DO NOT SUM TO THE DENOMINATOR")
    print("=" * 96)
    print("  complete/chromosome genomes in the profile : %d" % len(comp))

    arms = defaultdict(set)
    for a in comp:
        r = prof[a]
        absent = [x.strip() for x in (r.get("replicons_absent") or "").split(";") if x.strip()]
        bl = blocks_of(r)
        pgm = any(is_pgm(b) for b in bl)
        if pgm:
            arms["pgm"].add(a)
        if "NC_003132" in absent:
            arms["pPCP1"].add(a)
        if "NC_003131" in absent:
            arms["pCD1"].add(a)
        if "NC_003134" in absent:
            arms["pMT1"].add(a)
    named = set().union(*arms.values())
    noloss = {a for a in comp
              if not blocks_of(prof[a])
              and not (prof[a].get("replicons_absent") or "").strip()}
    outside = [a for a in comp if a not in named and a not in noloss]

    print("  in a named class                : %d" % len(named))
    print("  no block and no replicon loss   : %d" % len(noloss))
    print("  OUTSIDE BOTH ARMS               : %d   <-- the review's 102-152" % len(outside))
    print()
    cat = Counter()
    for a in outside:
        bl = blocks_of(prof[a])
        absent = [x.strip() for x in (prof[a].get("replicons_absent") or "").split(";") if x.strip()]
        if bl and not absent:
            cat["has >=10kb block(s), no replicon loss"] += 1
        elif absent and not bl:
            cat["replicon loss, no >=10kb block"] += 1
        elif bl and absent:
            cat["both"] += 1
        else:
            cat["neither (should not happen)"] += 1
    for k, v in cat.most_common():
        print("    %-42s %4d" % (k, v))
    coords = Counter()
    for a in outside:
        for b in blocks_of(prof[a]):
            p = b.split(":")
            coords[(p[0], int(p[1]) // 5000 * 5000, int(p[2]) // 5000 * 5000)] += 1
    print()
    print("  blocks in the OUTSIDE group: %d distinct coordinates, top 5:" % len(coords))
    for k, v in coords.most_common(5):
        print("    %-12s %9d +%-8d %4d genomes" % (k[0], k[1], k[2], v))
    print("  -> if this is non-empty, the pooled test is 'four named classes vs none',")
    print("     not 'any variant vs none'. The p value stands; the LABEL does not.")

    # ------------------------------------------------------------------ P0-2
    print()
    print("=" * 96)
    print("P0-2  CLUSTER 349: within-cluster distances and temporal signal")
    print("=" * 96)
    D = np.load("results/dedup/snp_distance_raw_tierB.npy").astype(float)
    order = [l.strip() for l in
             open("results/dedup/snp_distance_raw_tierB.order.txt") if l.strip()]
    idx = {a: i for i, a in enumerate(order)}
    target = None
    for r in csv.DictReader(open("results/dedup/tierB_clusters.tsv"), delimiter="\t"):
        if int(r.get("size", 0)) >= 400:
            target = r
            break
    if target is None:
        print("  no large cluster found")
        return
    members = [m.strip() for m in target["members"].split(",") if m.strip()]
    print("  cluster_id=%s  size=%s  diameter=%s"
          % (target["cluster_id"], target["size"], target.get("diameter")))
    m = [a for a in members if a in idx]
    print("  members present in the matrix: %d" % len(m))
    sub = D[np.ix_([idx[a] for a in m], [idx[a] for a in m])]
    iu = np.triu_indices(len(m), 1)
    d = sub[iu]
    d = d[d < SENT]
    print("  within-cluster pairwise distance (n=%d pairs):" % d.size)
    for p in (0, 5, 25, 50, 75, 95, 99, 100):
        print("    p%-3d %6.1f" % (p, np.percentile(d, p)))
    print("    mean %.2f  sd %.2f" % (d.mean(), d.std(ddof=1)))
    print("  NOTE: median/p95 are NOT bounded by the threshold t; only the")
    print("        maximum is (complete linkage guarantees it <= t by construction).")

    # Mantel-type permutation on year differences
    years = {}
    for a in m:
        y = (meta.get(a, {}).get("collection_date") or "")[:4]
        if y.isdigit():
            years[a] = int(y)
    dated = [a for a in m if a in years]
    print("  dated members: %d  (years %d-%d)"
          % (len(dated), min(years[a] for a in dated), max(years[a] for a in dated)))
    if len(dated) > 10:
        di = [idx[a] for a in dated]
        subd = D[np.ix_(di, di)][np.triu_indices(len(di), 1)]
        yd = np.array([years[a] for a in dated])
        ydiff = np.abs(yd[:, None] - yd[None, :])[np.triu_indices(len(di), 1)]
        keep = subd < SENT
        subd, ydiff = subd[keep], ydiff[keep]

        def spearman(x, y):
            rx = np.argsort(np.argsort(x)).astype(float)
            ry = np.argsort(np.argsort(y)).astype(float)
            rx -= rx.mean(); ry -= ry.mean()
            den = np.sqrt((rx ** 2).sum() * (ry ** 2).sum())
            return float((rx * ry).sum() / den) if den else 0.0

        rho = spearman(subd, ydiff)
        rng = np.random.default_rng(20261011)
        null = np.empty(9999)
        for i in range(9999):
            perm = rng.permutation(yd)
            pd_ = np.abs(perm[:, None] - perm[None, :])[np.triu_indices(len(di), 1)]
            pd_ = pd_[keep]
            null[i] = spearman(subd, pd_)
        lo, hi = np.percentile(null, [2.5, 97.5])
        p_two = 2 * min((null >= rho).mean(), (null <= rho).mean())
        print("  Mantel-type permutation (9,999 permutations of sampling year):")
        print("    observed Spearman rho(within-cluster SNP distance, |year diff|) = %.4f" % rho)
        print("    null 95%% interval = [%.4f, %.4f]" % (lo, hi))
        print("    two-sided p = %.4f" % p_two)
        print("    verdict: %s" % ("NO detectable temporal signal" if p_two > 0.05
                                   else "detectable temporal signal"))
    with open("results/qc/cluster349_internal.tsv", "w") as fh:
        fh.write("statistic\tvalue\n")
        fh.write("cluster_id\t%s\n" % target["cluster_id"])
        fh.write("size\t%d\n" % len(m))
        fh.write("pairs\t%d\n" % d.size)
        for p in (0, 5, 25, 50, 75, 95, 99, 100):
            fh.write("p%d\t%.2f\n" % (p, np.percentile(d, p)))
        fh.write("mean\t%.2f\n" % d.mean())
        if len(dated) > 10:
            fh.write("mantel_rho\t%.4f\n" % rho)
            fh.write("mantel_ci_lo\t%.4f\n" % lo)
            fh.write("mantel_ci_hi\t%.4f\n" % hi)
            fh.write("mantel_p\t%.4f\n" % p_two)
    print("  wrote results/qc/cluster349_internal.tsv")

    # ------------------------------------------------------------------ P1-3
    print()
    print("=" * 96)
    print("P1-3  THE FIVE dnadiff-EXCLUDED STRAINS: one clonal group or five?")
    print("=" * 96)
    five = ["GCF_015337085.2", "GCF_015337445.1", "GCF_015338205.1",
            "GCF_015337285.1", "GCF_000590535.2"]
    DA = np.load("results/dedup/snp_distance_raw_tierA.npy").astype(float)
    oa = [l.strip() for l in
          open("results/dedup/snp_distance_tierA.order.txt") if l.strip()]
    ia = {a: i for i, a in enumerate(oa)}
    present = [a for a in five if a in ia]
    if len(present) >= 2:
        ii = [ia[a] for a in present]
        within = DA[np.ix_(ii, ii)][np.triu_indices(len(ii), 1)]
        out = [i for i in range(len(oa)) if i not in ii]
        between = DA[np.ix_(ii, out)]
        print("  within-group distances  : %s" % sorted(int(x) for x in within))
        print("    max within            : %d" % within.max())
        print("  distance to any other   : min %d   median %d"
              % (between.min(), np.median(between)))
        print("  ratio min(between)/max(within) = %.1f"
              % (between.min() / within.max() if within.max() else float("inf")))
        print("  -> %s" % ("a single clonal group (within << between)"
                           if between.min() > within.max()
                           else "NOT separable as one group; report source only"))
    else:
        print("  fewer than 2 of the five are in the Tier A matrix: %s" % present)

    # ------------------------------------------------------------------ P1-4
    print()
    print("=" * 96)
    print("P1-4  STRUCTURAL-VARIANT FREQUENCIES, ALL vs LAB-NAMED REMOVED")
    print("=" * 96)
    rows = []
    for scale, accs in (("all_profile", list(prof)),
                        ("complete_level", comp),
                        ("complete_minus_lab_named",
                         [a for a in comp
                          if not LAB.search(meta.get(a, {}).get("strain") or "")])):
        n = len(accs)
        cnt = Counter()
        for a in accs:
            r = prof[a]
            absent = [x.strip() for x in (r.get("replicons_absent") or "").split(";") if x.strip()]
            if any(is_pgm(b) for b in blocks_of(r)):
                cnt["pgm"] += 1
            for k, (rep, _) in NAMED_CLASS.items():
                if k != "pgm" and rep in absent:
                    cnt[k] += 1
        rows.append((scale, n, cnt))
    print("  %-26s %6s %8s %8s %8s %8s" % ("set", "n", "pgm-", "pPCP1-", "pCD1-", "pMT1-"))
    for scale, n, c in rows:
        print("  %-26s %6d %7.1f%% %7.1f%% %7.1f%% %7.1f%%"
              % (scale, n, 100 * c["pgm"] / n, 100 * c["pPCP1"] / n,
                 100 * c["pCD1"] / n, 100 * c["pMT1"] / n))
    with open("results/qc/variant_freq_strata.tsv", "w") as fh:
        fh.write("set\tn\tpgm_pct\tpPCP1_pct\tpCD1_pct\tpMT1_pct\n")
        for scale, n, c in rows:
            fh.write("%s\t%d\t%.2f\t%.2f\t%.2f\t%.2f\n"
                     % (scale, n, 100 * c["pgm"] / n, 100 * c["pPCP1"] / n,
                        100 * c["pCD1"] / n, 100 * c["pMT1"] / n))
    print("  wrote results/qc/variant_freq_strata.tsv")


if __name__ == "__main__":
    main()
