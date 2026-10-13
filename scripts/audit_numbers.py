#!/usr/bin/env python3
"""audit_numbers.py — producer side of results/numbers.json.

Re-derives, from the artefacts already on disk, every headline number the
manuscript reports in Results 3.1-3.8 and the figure-prep tables, and appends
one JSON record per section to results/numbers.json.

Why this exists (external review 2026-10-13): a manuscript number that no file
and no line can produce is a silent failure of exactly the kind this paper is
about.  The audit that found the 9/161 -> 9/154 drift was done by hand three
times; with this producer the fourth audit is a diff.

This script is READ-ONLY with respect to the analysis artefacts.  It raises on
a missing file rather than skipping: a partial audit is worse than none.

Run:  python3 scripts/audit_numbers.py
"""
import csv
import os
import re
import sys
from collections import Counter, OrderedDict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _numbers import record                                    # noqa: E402

QC = "results/qc"
DEDUP = "results/dedup"
FIG = "results/figures"
NUM = "results/numbers.json"

norm = lambda a: re.sub(r"\.\d+$", "", (a or "").strip())       # noqa: E731
PGM = ("NC_003143", 2135000, 105000)
TAG2P = {"NC_003131": "pCD1", "NC_003132": "pPCP1", "NC_003134": "pMT1"}
LAB = re.compile(r"\b(EV\s*76|EV\s*NIIEG|EV76|KIM|KIM10|A1122|CO92|Pasteurella|"
                 r"Pestoides|Antiqua|Nepal516|Angola|C-?792|CMCC|SCPM|Tjiwidej|"
                 r"Java9|Harbin35|ATCC|NCTC|Biotype|vaccine|attenuat)\b", re.I)


def rd(path, delim="\t"):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter=delim))


def fisher(t):
    from math import comb
    (a, b), (c, d) = t
    n, r1, r2, c1 = a + b + c + d, a + b, c + d, a + c
    lo, hi = max(0, c1 - r2), min(r1, c1)
    hyp = lambda x: (comb(r1, x) * comb(r2, c1 - x) / comb(n, c1)   # noqa: E731
                     if lo <= x <= hi else 0.0)
    obs = hyp(a)
    p = sum(hyp(x) for x in range(lo, hi + 1) if hyp(x) <= obs * (1 + 1e-9))
    aa, bb, cc, dd = a, b, c, d
    if min(aa, bb, cc, dd) == 0:
        aa, bb, cc, dd = aa + .5, bb + .5, cc + .5, dd + .5
    return p, (aa * dd) / (bb * cc)


def wilson(k, n, z=1.959964):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * ((p * (1 - p) / n + z * z / (4 * n * n)) ** .5) / d
    return (max(0.0, c - h), min(1.0, c + h))


def main():
    out = {}

    # ---------------------------------------------------------------- sets
    prof = {norm(r["accession"]): r for r in rd(f"{QC}/deletion_profile_all.tsv")}
    meta = {}
    for r in rd(f"{QC}/assembly_metadata.tsv"):
        meta[norm(r["accession"])] = r
        p = norm(r.get("paired_accession") or "")
        if p:
            meta.setdefault(p, r)
    tierA = [l.split("\t")[0] for l in open(f"{QC}/tierA_refined.tsv")][1:]
    tierB = [l.split("\t")[0] for l in open(f"{QC}/tierB_candidates.tsv")][1:]
    sitesA = sum(1 for _ in open(f"{DEDUP}/sites_tierA.tsv")) - 1
    sitesB = sum(1 for _ in open(f"{DEDUP}/sites_tierB.tsv")) - 1
    out["sets"] = {"p3441": len(prof), "tierA": len(tierA), "tierB": len(tierB),
                   "sites_tierA": sitesA, "sites_tierB": sitesB}

    # ------------------------------------------------- 3.1 noise floor
    _sb = rd(f"{FIG}/fig1c_self_blocks.csv", ",")
    selflen = [int(float(r["length_bp"])) for r in _sb]
    chrlen = [int(float(r["length_bp"])) for r in _sb if r["replicon"] == "NC_003143"]
    ss = rd(f"{FIG}/fig1b_same_strain.csv", ",")
    scal = {r["statistic"]: float(r["value"]) for r in rd(f"{FIG}/fig1b_scalars.csv", ",")}
    d = np.load(f"{DEDUP}/snp_distance_raw_tierB.npy")
    v = d[np.triu_indices_from(d, k=1)].astype(np.int64)
    out["noise_floor"] = {
        "self_runs_total_all_replicons": len(selflen),
        "self_runs_ge500bp_all_replicons": int(sum(1 for x in selflen if x >= 500)),
        "self_runs_total_chromosome": len(chrlen),
        "self_runs_ge500bp_chromosome": int(sum(1 for x in chrlen if x >= 500)),
        "self_runs_shortest_chromosome_bp": min(chrlen),
        "self_max_run_bp": max(selflen),
        "self_runs_ge10kb": int(sum(1 for x in selflen if x >= 10000)),
        "same_strain_tierA": {r["strain"]: int(r["raw_distance_tierA"]) for r in ss},
        "same_strain_tierB": {r["strain"]: int(r["raw_distance_tierB"]) for r in ss},
        "tierB_pairs": int(scal["n_pairs"]),
        "tierB_zero": int(scal["n_zero"]),
        "tierB_le25": int(scal["n_le_25"]),
        "tierB_le19": int((v <= 19).sum()),
        "tierB_pct_le19": round(100.0 * (v <= 19).sum() / v.size, 2),
        "tierB_percentiles": {("p%d" % q): float(np.percentile(v, q, method="linear"))
                              for q in (1, 5, 25, 50, 75, 95, 99)},
    }

    # ------------------------------------------------- 3.2 QC metrics
    flags = rd(f"{FIG}/fig2a_metric_flags.csv", ",")
    byflag = Counter(r["metric"] for r in flags if r["flagged"] == "1")
    checkm = rd(f"{FIG}/fig2d_checkm.csv", ",")
    out["quality_metrics"] = {
        "n_rows": len(flags), "n_genomes": len({r["accession"] for r in flags}),
        "flagged_checkm": byflag["checkm"], "flagged_called_frac": byflag["called_frac"],
        "flagged_dnadiff_ref": byflag["dnadiff_ref"],
        "checkm_min_anchor": min(float(r["completeness"]) for r in checkm),
    }
    strat = {r["set"]: {k: (float(r[k]) if k != "n" else int(r["n"]))
                        for k in ("n", "pgm_pct", "pPCP1_pct", "pCD1_pct", "pMT1_pct")}
             for r in rd(f"{QC}/variant_freq_strata.tsv")}
    sub = sum(1 for r in rd(f"{QC}/replicon_three_state_seqreport.tsv")
              if norm(r["accession"]) in set(map(norm, prof)) and int(r["n_plasmid_seq"]) > 0)
    out["archive_composition"] = {"strata_pct": strat,
                                  "p258_plasmid_submitted": sub, "p258_n": len(prof)}
    t2 = [r for r in rd(f"{FIG}/fig2c_dnadiff.csv", ",")]
    out["quality_metrics"]["dnadiff_rows"] = len(t2)

    # ------------------------------------------------- 3.3 spatial
    rec = rd(f"{FIG}/fig3b_recurrence.csv", ",")
    blocks = rd(f"{QC}/uncalled_blocks_all.tsv")
    pgm_locus = sum(1 for r in blocks
                    if r["replicon"] == PGM[0] and int(r["start"]) == PGM[1]
                    and int(r["length"]) >= 100000)
    pgm_bin = [r for r in rec if r["replicon"] == PGM[0]
               and int(r["start"]) == PGM[1] and int(r["length_bp"]) == PGM[2]]
    out["spatial"] = {
        "distinct_coordinates": len(rec),
        "seen_once": sum(1 for r in rec if int(r["n_genomes"]) == 1),
        "recur_ge10": sum(1 for r in rec if int(r["n_genomes"]) >= 10),
        "pgm_locus_n": pgm_locus,
        "pgm_bin_n": int(pgm_bin[0]["n_genomes"]) if pgm_bin else None,
    }

    # ------------------------------------------------- 3.5 / 3.6 polarity
    pol = rd(f"{FIG}/fig5a_region_polarity.csv", ",")
    bg = [float(r["delta"]) for r in rd(f"{FIG}/fig5b_background.csv", ",")]
    out["polarity"] = {
        "merged_regions": len(pol),
        "polarity_counts": dict(Counter(r["polarity"] for r in pol)),
        "raw_blocks": len(rd(f"{FIG}/fig5b_block_delta.csv", ",")),
        "background_n": len(bg),
        "background_median": round(float(np.median(bg)), 4),
        "background_p75": round(float(np.percentile(bg, 75)), 4),
        "background_p90": round(float(np.percentile(bg, 90)), 4),
        "background_p95": round(float(np.percentile(bg, 95)), 4),
        "background_p99": round(float(np.percentile(bg, 99)), 4),
    }

    # ------------------------------------------------- 3.4 dereplication
    D = np.load(f"{DEDUP}/snp_distance_raw_tierA.npy").astype(np.float64)
    np.fill_diagonal(D, 0)
    iu = np.triu_indices(D.shape[0], 1)
    eligible = int((D[iu] <= 20).sum())
    # blocking recomputed with the pipeline's own rules
    ts = {norm(r["accession"]): (r["pCD1_state"], r["pPCP1_state"], r["pMT1_state"])
          for r in rd(f"{QC}/replicon_three_state_seqreport.tsv")}
    A = [l.strip() for l in open(f"{DEDUP}/snp_distance_raw_tierA.order.txt")]
    P = [prof.get(norm(a)) for a in A]
    sig = [None if r is None else (r.get("replicons_absent") or "",
                                   r.get("big_blocks_5kb") or "") for r in P]
    chrom = [None if r is None else frozenset(
                 b for b in (r.get("big_blocks_5kb") or "").split(";")
                 if b.startswith("NC_003143:")) for r in P]
    blocked = {}
    for mode in ("full", "threestate", "conservative"):
        tot = el = 0
        for i in range(len(A)):
            if sig[i] is None:
                continue
            for j in range(i + 1, len(A)):
                if sig[j] is None:
                    continue
                if mode == "full":
                    blk = sig[i] != sig[j]
                else:
                    ta, tb = ts.get(norm(A[i])), ts.get(norm(A[j]))
                    blk = chrom[i] != chrom[j]
                    if mode == "threestate" and not blk and ta and tb:
                        blk = any(ta[k] != "UNRESOLVED" and tb[k] != "UNRESOLVED"
                                  and ta[k] != tb[k] for k in range(3))
                    elif mode == "conservative" and ta is not None and tb is not None:
                        blk = blk or (ta != tb)
                if blk:
                    tot += 1
                    if D[i, j] <= 20:
                        el += 1
        blocked[mode] = {"total": tot, "eligible": el}
    mm = {}
    for line in open(f"{QC}/dedup_modes_summary.txt"):
        p = line.split()
        if len(p) == 4 and p[0] in ("off", "full", "threestate", "conservative"):
            mm[p[0]] = {"kept": int(p[1]), "clusters": int(p[2]), "multi": int(p[3])}
    cls = rd(f"{QC}/dedup_mixed_classification.tsv")
    out["dereplication"] = {
        "constraints": mm, "blocked_pairs": blocked, "pairs_total": int(iu[0].size),
        "eligible_pairs": eligible,
        "mixed_clusters": dict(Counter(r["category"] for r in cls)),
    }

    # ------------------------------------------------- 3.7 the 419 set
    yrs = {r["statistic"]: r["value"] for r in rd(f"{FIG}/fig1d_cluster_years.csv", ",")}
    ci = rd(f"{FIG}/fig1d_cluster_internal.csv", ",")
    out["cluster_419"] = {
        "n_members": int(yrs["n_members"]), "year_min": int(yrs["year_min"]),
        "year_max": int(yrs["year_max"]), "n_dated": int(yrs["n_dated"]),
        "p50": int(yrs["p50"]), "p95": int(yrs["p95"]),
        "pairs": int(sum(int(r["n_pairs"]) for r in ci)),
    }

    # ------------------------------------------------- 3.8 directionality
    P258 = [a for a in prof
            if (meta.get(a, {}).get("assembly_level") or "").lower()
            .startswith(("complete", "chromosome"))]
    def blocks_of(a):
        return [b for b in (prof[a].get("big_blocks_5kb") or "").split(";") if b.count(":") == 2]
    def bkey(b):
        p = b.split(":")
        return (p[0], int(p[1]), int(p[2]))
    pgm_arm = [a for a in P258 if any(bkey(b) == PGM for b in blocks_of(a))]
    named = [a for a in P258
             if any(bkey(b) == PGM for b in blocks_of(a))
             or any(t in (prof[a].get("replicons_absent") or "") for t in TAG2P)]
    other = [a for a in P258 if a not in named and blocks_of(a)]
    novar = [a for a in P258 if a not in named and not blocks_of(a)]
    lab = lambda a: bool(LAB.search(meta.get(a, {}).get("strain") or ""))   # noqa: E731
    k1, n1 = sum(lab(a) for a in pgm_arm), len(pgm_arm)
    k2, n2 = sum(lab(a) for a in other), len(other)
    P, OR = fisher([[k1, n1 - k1], [k2, n2 - k2]])
    out["directionality"] = {
        "p258": len(P258), "named_class": len(named), "other_block": len(other),
        "no_variant": len(novar), "pgm_arm_n": n1, "pgm_arm_labnamed": k1,
        "arm2_n": n2, "arm2_labnamed": k2,
        "fisher_p": float("%.6g" % P), "odds_ratio": round(OR, 3),
        "wilson_arm1": [round(100 * x, 1) for x in wilson(k1, n1)],
        "wilson_arm2": [round(100 * x, 1) for x in wilson(k2, n2)],
        "labnamed_p258": sum(lab(a) for a in P258),
        "labnamed_with_named_class": sum(lab(a) for a in named),
    }

    # ------------------------------------------------- pseudogene entanglement
    pe = rd(f"{QC}/pseudogene_entanglement.tsv")
    cat = Counter(r["category"] for r in pe)
    out["pseudogene_entanglement"] = {
        "total": len(pe), "none": cat["no_bakta_CDS"],
        "shorter": cat["shorter_bakta_CDS"], "full": cat["full_length_bakta_CDS"],
        "pct_none": round(100 * cat["no_bakta_CDS"] / len(pe), 1),
        "pct_none_or_shorter": round(100 * (cat["no_bakta_CDS"] + cat["shorter_bakta_CDS"]) / len(pe), 1),
    }

    record("audit_numbers.py", out, path=NUM)
    import json
    print(json.dumps(out, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
