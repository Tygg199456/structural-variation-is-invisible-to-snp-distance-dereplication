#!/usr/bin/env python3
"""Run Tier A dedup under three profile-constraint modes and classify the
mixed clusters (external review 2026-10-11, execution steps 2-3-4).

Modes (env SY_DEDUP_PROFILE_MODE, implemented in snp_dedup.py):
  off         no profile constraint          -> baseline 175 / 36 / 25
  full        current signature constraint    -> reproduces the shipped 206
  threestate  three-state-aware constraint    -> new 175 -> ?

Then, from the OFF-run clusters, every multi-member cluster is classified:
  a  chromosome-block signatures differ            (real)
  b  only plasmid differences, >=1 side UNRESOLVED (submission artefact)
  c  only plasmid differences, both sides resolved (real plasmid loss)
True mixed = a + c.

Outputs
  results/dedup/tierA_keep.{off,full,threestate}.tsv
  results/dedup/tierA_clusters.{off,full,threestate}.tsv
  results/qc/dedup_mixed_classification.tsv
  results/qc/dedup_modes_summary.txt
"""
import csv, os, re, sys, importlib.util, collections
import os

ROOT = os.environ.get("PLAGUE_DECAY_ROOT",
                       os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.chdir(ROOT)
sys.path.insert(0, "scripts")

WHITELIST = {
    "CO92": {"accession": "GCF_000009065.1"},
    "KIM": {"accession": "GCA_000006645.1"},
    "Antiqua": {"accession": "GCF_000013825.1"},
    "Nepal516": {"accession": "GCF_000013805.1"},
    "str91001": {"accession": "GCF_000007885.1"},
    "PestoidesF": {"accession": "GCF_000016445.1"},
    "A1122": {"accession": "GCF_000222975.1"},
    "Angola": {"accession": "GCF_000018805.1"},
}
PRIORITY = ["refseq_gcf", "has_collection_date", "has_geo_loc_name",
            "closest_to_median_length"]
MAXD = 20
norm = lambda a: re.sub(r"\.\d+$", "", str(a).strip())

TS_PATH = "results/qc/replicon_three_state_seqreport.tsv"
tstate = {norm(r["accession"]): (r["pCD1_state"], r["pPCP1_state"], r["pMT1_state"])
          for r in csv.DictReader(open(TS_PATH), delimiter="\t")}
prof = {}
for r in csv.DictReader(open("results/qc/deletion_profile_all.tsv"), delimiter="\t"):
    prof[r["accession"]] = ((r.get("replicons_absent") or ""), (r.get("big_blocks_5kb") or ""))


class NS:
    pass


def ns(**kw):
    o = NS()
    for k, v in kw.items():
        setattr(o, k, v)
    return o


def run(mode):
    os.environ["SY_DEDUP_PROFILE_MODE"] = mode
    os.environ["SY_DEDUP_THREESTATE"] = TS_PATH
    spec = importlib.util.spec_from_file_location("sd_" + mode, "scripts/snp_dedup.py")
    sd = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(sd)
    sn = ns(input=ns(dist="results/dedup/snp_distance_raw_tierA.npy",
                     cand="results/qc/tierA_refined.tsv",
                     meta="results/qc/assembly_metadata.tsv",
                     profile="results/qc/deletion_profile_all.tsv"),
            output=ns(keep=f"results/dedup/tierA_keep.{mode}.tsv",
                      clusters=f"results/dedup/tierA_clusters.{mode}.tsv"),
            params=ns(max_dist=MAXD, priority=PRIORITY, whitelist=WHITELIST),
            log=[f"logs/02_tierA_dedup.{mode}.log"])
    sd.snakemake = sn
    sd.main()
    n = sum(1 for _ in open(sn.output.keep)) - 1
    rows = list(csv.DictReader(open(sn.output.clusters), delimiter="\t"))
    multi = [r for r in rows if int(r["size"]) > 1]
    return n, len(rows), len(multi), rows


out = []
for mode in ("off", "full", "threestate", "conservative"):
    n, ncl, nmulti, rows = run(mode)
    out.append((mode, n, ncl, nmulti, rows))
    print("MODE %-11s kept=%-4d clusters=%-4d multi=%d" % (mode, n, ncl, nmulti))

# ---- classify the OFF-run multi-member clusters -----------------------------
off_rows = out[0][4]


def chrom(a):
    p = prof.get(a)
    return frozenset(b for b in p[1].split(";") if b.startswith("NC_003143:")) if p else None


def full_sig(a):
    return prof.get(a)


cls_rows = []
cnt = collections.Counter()
for r in off_rows:
    if int(r["size"]) <= 1:
        continue
    mem = [m for m in r["members"].split(",") if m]
    mixed = len({full_sig(m) for m in mem}) > 1
    if not mixed:
        cat = "not_mixed"
    else:
        chrom_diff = len({chrom(m) for m in mem}) > 1
        if chrom_diff:
            cat = "a_chromosome"
        else:
            # only plasmid-level differences
            plasmid_diff_unres = False
            plasmid_diff_res = False
            for k in range(3):
                st = [tstate.get(norm(m), ("?", "?", "?"))[k] for m in mem]
                if len(set(st)) > 1:
                    if any(s == "UNRESOLVED" for s in st):
                        plasmid_diff_unres = True
                    else:
                        plasmid_diff_res = True
            cat = "b_plasmid_unresolved" if plasmid_diff_unres else (
                "c_plasmid_resolved" if plasmid_diff_res else "other")
    cnt[cat] += 1
    cls_rows.append({"cluster_id": r["cluster_id"], "size": r["size"],
                     "representative": r["representative"], "category": cat,
                     "members": r["members"]})

with open("results/qc/dedup_mixed_classification.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["cluster_id", "size", "representative",
                                       "category", "members"], delimiter="\t")
    w.writeheader()
    w.writerows(cls_rows)

print()
print("=== OFF-run multi-member cluster classification ===")
nmixed = sum(v for k, v in cnt.items() if k != "not_mixed")
for k in ("a_chromosome", "b_plasmid_unresolved", "c_plasmid_resolved",
          "not_mixed", "other"):
    if cnt[k]:
        print("  %-24s %3d" % (k, cnt[k]))
print("  mixed (old profile signature) = %d ; true mixed (a+c) = %d"
      % (nmixed, cnt["a_chromosome"] + cnt["c_plasmid_resolved"]))

with open("results/qc/dedup_modes_summary.txt", "w") as fh:
    fh.write("mode\tkept\tclusters\tmulti_member\n")
    for mode, n, ncl, nmulti, _ in out:
        fh.write("%s\t%d\t%d\t%d\n" % (mode, n, ncl, nmulti))
    fh.write("\ncategory\tcount\n")
    for k, v in sorted(cnt.items()):
        fh.write("%s\t%d\n" % (k, v))
print()
print("WROTE results/qc/dedup_mixed_classification.tsv and dedup_modes_summary.txt")
