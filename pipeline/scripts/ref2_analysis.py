#!/usr/bin/env python3
"""ref2_analysis.py — second-reference sensitivity analysis (SY1011 / Paper 1).

Replicates the pipeline's own chain (union variant sites -> reference-length
genotype matrix -> raw distance -> complete-linkage dedup at t = 20 under the
four deletion-profile constraints) for the SAME 256 Tier A genomes mapped
against a second reference (GCF_000013825.1, Y. pestis Antiqua), and reports the
two endpoints pre-registered in Methods 2.9:

  E1  do the two records labelled CO92 merge at the rule-selected threshold?
  E2  is the proportion of multi-member clusters with discordant resolved
      deletion content above one half?

The threshold is re-selected in the new framework by the same rule (smallest t
with all five same-strain pairs merged, diameter <= 2t, largest cluster <= 25)
rather than carried over.

Reads:  results_ref2/snippy/<acc>/snps.{vcf,aligned.fa}   (written by sy_ref2snippy)
Writes: results_ref2/numbers.json  + stdout report.

Run:  cd pipeline/plague-decay && python3 scripts/ref2_analysis.py
"""
import csv
import gzip
import json
import os
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _numbers import record                                    # noqa: E402

REF = "GCF_000013825.1"
SNP = "results_ref2/snippy"
QC = "results/qc"
OUT = "results_ref2/numbers.json"
CODE = {"A": 1, "C": 2, "G": 3, "T": 4}
MISSING = set("Nn-.*")
SENTINEL = np.int32(10 ** 6)
T_GRID = [0, 3, 5, 10, 15, 20, 25, 30, 50, 75, 100]   # the pipeline's own sweep grid
norm = lambda a: re.sub(r"\.\d+$", "", (a or "").strip())       # noqa: E731
TAG2P = {"pCD1", "pPCP1", "pMT1"}


def open_maybe_gz(p):
    p = str(p)
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def read_aligned_fa(path):
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


def uncalled_blocks(seqs, min_bp=10000, bin_bp=5000):
    """Maximal [Nn-] runs >= min_bp, binned to 5 kb -> 'replicon:bin_start:bin_len'."""
    out = []
    for chrom, seq in seqs.items():
        n = len(seq)
        i = 0
        while i < n:
            if seq[i] in MISSING:
                j = i
                while j < n and seq[j] in MISSING:
                    j += 1
                if j - i >= min_bp:
                    lo = (i // bin_bp) * bin_bp
                    hi = ((j + bin_bp - 1) // bin_bp) * bin_bp
                    out.append("%s:%d:%d" % (chrom, lo, hi - lo))
                i = j
            else:
                i += 1
    return out


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


def main():
    accs = sorted(d.name for d in Path(SNP).iterdir() if (d / "snps.vcf").exists())
    print("[ref2] samples with snps.vcf: %d" % len(accs))
    missing_aln = [a for a in accs if not (Path(SNP) / a / "snps.aligned.fa").exists()]
    if missing_aln:
        raise SystemExit("FATAL: %d samples missing snps.aligned.fa (e.g. %s)"
                         % (len(missing_aln), missing_aln[:3]))

    # ---- pass 1: union of variant sites
    sites = set()
    for a in accs:
        with open_maybe_gz(Path(SNP) / a / "snps.vcf") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.split("\t", 3)
                if len(f) >= 3:
                    sites.add((f[0], int(f[1])))
    sites = sorted(sites)
    by_chrom = {}
    for col, (c, p) in enumerate(sites):
        by_chrom.setdefault(c, []).append((p, col))
    print("[ref2] union variant sites = %d" % len(sites))

    # ---- pass 2: genotype matrix + deletion profiles
    G = np.zeros((len(accs), len(sites)), dtype=np.int8)
    called = np.zeros(len(accs), dtype=np.float32)
    blocks = []
    for si, a in enumerate(accs):
        seqs = read_aligned_fa(Path(SNP) / a / "snps.aligned.fa")
        row = G[si]
        for chrom, poss in by_chrom.items():
            seq = seqs.get(chrom)
            if seq is None:
                continue
            L = len(seq)
            for pos, col in poss:
                if pos <= L:
                    ch = seq[pos - 1]
                    if ch not in MISSING:
                        row[col] = CODE.get(ch.upper(), 0)
        called[si] = float((row != 0).sum()) / len(sites)
        blocks.append(uncalled_blocks(seqs))
        if (si + 1) % 64 == 0:
            print("[ref2] pass2 %d/%d" % (si + 1, len(accs)))

    # ---- raw distance, exactly as snp_distance.py
    called_mask = (G != 0).astype(np.float32)
    C = called_mask @ called_mask.T
    onehot = np.zeros((G.shape[0], G.shape[1], 4), dtype=np.float32)
    for k in range(4):
        onehot[:, :, k] = (G == k + 1)
    agree = np.zeros_like(C)
    for k in range(4):
        X = onehot[:, :, k]
        agree += X @ X.T
    diff = C - agree
    D_raw = np.rint(diff).astype(np.int32)
    joint_frac = C / float(len(sites))
    blocked_lowcov = joint_frac < 0.70
    np.fill_diagonal(D_raw, 0)
    D = D_raw.copy()
    D[blocked_lowcov] = SENTINEL

    # ---- noise floor in the NEW framework (5 same-strain pairs)
    ss = list(csv.DictReader(open("results/figures/fig1b_same_strain.csv")))
    idx = {a: i for i, a in enumerate(accs)}
    floor = {}
    for r in ss:
        a, b = r["acc_a"], r["acc_b"]
        floor[r["strain"]] = int(D_raw[idx[a], idx[b]]) if a in idx and b in idx else None
    print("[ref2] same-strain distances in the new framework: %s" % floor)

    # ---- threshold selection by the same rule
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    Sm = squareform(np.where(D >= SENTINEL, 1e6, D).astype(np.float64), checks=False)
    Z = linkage(Sm, method="complete")
    selected, rule_tbl = None, []
    for t in T_GRID:
        lab = fcluster(Z, t=t, criterion="distance")
        sizes = np.bincount(lab)[1:]
        diam = 0
        for c in np.unique(lab):
            m = np.where(lab == c)[0]
            if m.size > 1:
                sub = D_raw[np.ix_(m, m)]
                diam = max(diam, int(sub.max()))
        merged = 0
        for r in ss:
            if r["acc_a"] in idx and r["acc_b"] in idx and lab[idx[r["acc_a"]]] == lab[idx[r["acc_b"]]]:
                merged += 1
        ok = (merged == 5) and (diam <= 2 * t) and (sizes.max() <= 25)
        rule_tbl.append({"t": t, "merged_pairs": merged, "max_diameter": diam,
                         "max_size": int(sizes.max()), "pass": bool(ok)})
        if ok and selected is None:
            selected = t
    print("[ref2] rule-selected threshold = %s" % selected)
    if selected is None:
        # no threshold satisfies the rule: report the smallest t that merges all 5
        for row in rule_tbl:
            if row["merged_pairs"] == 5:
                selected = row["t"]
                print("[ref2] NO threshold satisfies the rule; smallest merging t = %d "
                      "(max_size=%d)" % (selected, row["max_size"]))
                break
    lab = fcluster(Z, t=selected, criterion="distance")

    # ---- deletion-profile constraints (pipeline rules)
    prof = {}
    for r in csv.DictReader(open(f"{QC}/deletion_profile_all.tsv"), delimiter="\t"):
        prof[norm(r["accession"])] = (r.get("replicons_absent") or "",
                                      r.get("big_blocks_5kb") or "")
    ts = {norm(r["accession"]): (r["pCD1_state"], r["pPCP1_state"], r["pMT1_state"])
          for r in csv.DictReader(open(f"{QC}/replicon_three_state_seqreport.tsv"),
                                  delimiter="\t")}
    tstate = {a: ts.get(norm(a)) for a in accs}

    # block signature in the NEW reference: use the ref2 blocks, plus plasmid state
    sig2 = {}
    for i, a in enumerate(accs):
        pg = prof.get(norm(a))
        sig2[i] = (frozenset(blocks[i]),            # new-reference chromosomal blocks
                   (pg[0] if pg else ""))           # plasmid inventory (reference-independent)

    def blocked(i, j, mode):
        ta, tb = tstate[accs[i]], tstate[accs[j]]
        if mode == "full":
            return sig2[i] != sig2[j] or (ta != tb)
        chrom_diff = sig2[i][0] != sig2[j][0]
        if mode == "threestate":
            if chrom_diff:
                return True
            if ta and tb:
                return any(ta[k] != "UNRESOLVED" and tb[k] != "UNRESOLVED"
                           and ta[k] != tb[k] for k in range(3))
            return False
        # conservative
        if chrom_diff:
            return True
        return ta is not None and tb is not None and ta != tb

    MULTI = {}
    for mode in ("off", "full", "threestate", "conservative"):
        Dm = D.copy()
        if mode != "off":
            for i in range(len(accs)):
                for j in range(i + 1, len(accs)):
                    if Dm[i, j] >= SENTINEL:
                        continue
                    if blocked(i, j, mode):
                        Dm[i, j] = Dm[j, i] = SENTINEL
        Zm = linkage(squareform(np.where(Dm >= SENTINEL, 1e6, Dm).astype(np.float64),
                                checks=False), method="complete")
        lm = fcluster(Zm, t=selected, criterion="distance")
        sizes = np.bincount(lm)[1:]
        multi = [int(np.where(lm == c)[0].size) for c in np.unique(lm) if (lm == c).sum() > 1]
        # discordant clusters under the conservative rule
        disc = 0
        for c in np.unique(lm):
            m = np.where(lm == c)[0]
            if m.size < 2:
                continue
            bad = any(blocked(int(m[x]), int(m[y]), "conservative")
                      for x in range(m.size) for y in range(x + 1, m.size))
            if bad:
                disc += 1
        MULTI[mode] = {"clusters": int(len(sizes)), "multi_member": len(multi),
                       "discordant_conservative": disc,
                       "max_size": int(sizes.max())}
    # one representative is retained per cluster, so kept == clusters
    for m in MULTI:
        MULTI[m]["kept"] = MULTI[m]["clusters"]

    # ---- endpoints
    a1, a2 = "GCF_000009065.1", "GCF_001293415.1"
    e1_dist = int(D_raw[idx[a1], idx[a2]]) if a1 in idx and a2 in idx else None
    e1_merge_off = bool(e1_dist is not None and e1_dist <= selected)
    e1_merge_cons = bool(e1_merge_off and
                         not blocked(idx[a1], idx[a2], "conservative"))
    off = MULTI["off"]
    e2_frac = (off["discordant_conservative"] / off["multi_member"]
               if off["multi_member"] else None)
    rep = {
        "reference": REF,
        "n_samples": len(accs),
        "union_variant_sites": len(sites),
        "called_frac": {"min": round(float(called.min()), 3),
                        "median": round(float(np.median(called)), 3),
                        "max": round(float(called.max()), 3)},
        "same_strain_distances_new_framework": floor,
        "threshold_rule_table": rule_tbl,
        "selected_t": selected,
        "constraints": MULTI,
        "endpoint_E1_co92_pair": {"raw_distance": e1_dist,
                                  "merges_at_t_off": e1_merge_off,
                                  "merges_at_t_conservative": e1_merge_cons},
        "endpoint_E2_discordant_fraction_off": e2_frac,
        "endpoint_E2_discordant_of_multi_off": [off["discordant_conservative"],
                                                off["multi_member"]],
    }
    Path("results_ref2").mkdir(parents=True, exist_ok=True)
    record("ref2_analysis.py", rep, path=OUT)
    print(json.dumps(rep, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
