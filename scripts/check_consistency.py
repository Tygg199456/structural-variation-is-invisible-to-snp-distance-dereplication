#!/usr/bin/env python3
"""Consistency assertions for the archive.  Non-zero exit on any failure.

Checks, in order:
  1  every results/figures/*.csv has a .prov provenance line, and vice versa
  2  figures/SHA256SUMS.txt hashes match the files on disk
  3  MANIFEST.sha256 hashes match the files on disk
  4  results/numbers.json parses and still holds the values the manuscript reports
     (FROZEN below — this is the assertion side of the number registry)

Run from the repository root:  python3 scripts/check_consistency.py
"""
import csv, glob, hashlib, json, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
fails = []


def sha256(path, buf=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(buf)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def check(cond, msg):
    print(("  ok   " if cond else "  FAIL ") + msg)
    if not cond:
        fails.append(msg)


# 1 ------------------------------------------------------------------ provenance
print("[1] results/figures CSV <-> .prov")
csvs = sorted(glob.glob(os.path.join(ROOT, "results/figures/*.csv")))
for c in csvs:
    check(os.path.exists(c + ".prov"), "prov exists for %s" % os.path.basename(c))
for p in sorted(glob.glob(os.path.join(ROOT, "results/figures/*.prov"))):
    check(os.path.exists(p[:-5]), "csv exists for %s" % os.path.basename(p))

# 2 ------------------------------------------------------------ figure manifest
print("[2] figures/SHA256SUMS.txt")
for line in open(os.path.join(ROOT, "figures/SHA256SUMS.txt")):
    line = line.strip()
    if not line or line.startswith("#"):
        continue
    h, rel = line.split("  ", 1)
    f = os.path.join(ROOT, rel)
    check(os.path.exists(f), "present: %s" % rel)
    if os.path.exists(f):
        check(sha256(f) == h, "hash matches: %s" % rel)

# 3 ------------------------------------------------------------------- manifest
print("[3] MANIFEST.sha256")
miss = bad = 0
for line in open(os.path.join(ROOT, "MANIFEST.sha256")):
    line = line.rstrip("\n")
    if not line:
        continue
    h, rel = line.split("  ", 1)
    f = os.path.join(ROOT, rel)
    if not os.path.exists(f):
        miss += 1
        continue
    if sha256(f) != h:
        bad += 1
check(miss == 0, "no missing files listed in MANIFEST (%d missing)" % miss)
check(bad == 0, "all MANIFEST hashes match (%d mismatched)" % bad)

# 4 --------------------------------------------------------- frozen manuscript values
print("[4] results/numbers.json frozen values")
latest = {}
for line in open(os.path.join(ROOT, "results/numbers.json")):
    line = line.strip()
    if line:
        rec = json.loads(line)
        for k, v in rec["numbers"].items():
            latest[k] = v

FROZEN = {
    ("sets", "tierA"): 256,
    ("sets", "tierB"): 3441,
    ("sets", "sites_tierA"): 12926,
    ("sets", "sites_tierB"): 67711,
    ("noise_floor", "self_runs_total_chromosome"): 201,
    ("noise_floor", "self_runs_ge500bp_chromosome"): 157,
    ("noise_floor", "self_runs_total_all_replicons"): 216,
    ("noise_floor", "self_runs_ge500bp_all_replicons"): 165,
    ("noise_floor", "self_runs_ge10kb"): 0,
    ("noise_floor", "self_max_run_bp"): 5553,
    ("noise_floor", "tierB_zero"): 2791,
    ("noise_floor", "tierB_le25"): 481230,
    ("noise_floor", "tierB_le19"): 301581,
    ("quality_metrics", "flagged_checkm"): 109,
    ("quality_metrics", "flagged_dnadiff_ref"): 5,
    ("quality_metrics", "flagged_called_frac"): 12,
    ("spatial", "distinct_coordinates"): 161,
    ("spatial", "pgm_locus_n"): 477,
    ("spatial", "pgm_bin_n"): 476,
    ("polarity", "merged_regions"): 13,
    ("polarity", "background_n"): 1000,
    ("polarity", "background_p95"): 0.0902,
    ("cluster_419", "n_members"): 419,
    ("cluster_419", "pairs"): 87571,
    ("directionality", "pgm_arm_n"): 41,
    ("directionality", "pgm_arm_labnamed"): 9,
    ("directionality", "arm2_n"): 154,
    ("directionality", "arm2_labnamed"): 9,
    ("directionality", "labnamed_p258"): 26,
    ("pseudogene_entanglement", "pct_none"): 30.0,
    ("pseudogene_entanglement", "pct_none_or_shorter"): 52.5,
}
for (sec, key), want in FROZEN.items():
    got = latest.get(sec, {}).get(key)
    check(got == want, "%s.%s == %r" % (sec, key, want) + ("" if got == want else "  (got %r)" % got))

# a few non-frozen but range-checked quantities
d = latest.get("dereplication", {})
check(d.get("constraints", {}).get("off", {}).get("kept") == 175, "dedup off kept == 175")
check(d.get("constraints", {}).get("conservative", {}).get("kept") == 209, "dedup conservative kept == 209")
check(d.get("eligible_pairs") == 289, "eligible pairs == 289")
check(d.get("blocked_pairs", {}).get("conservative", {}).get("eligible") == 163,
      "conservative blocked among eligible == 163")
fd = latest.get("directionality", {})
check(abs(fd.get("fisher_p", 0) - 0.00394942) < 1e-9, "fisher_p reproduces")
check(abs(fd.get("odds_ratio", 0) - 4.531) < 1e-6, "odds_ratio reproduces")

print()
if fails:
    print("%d CHECK(S) FAILED" % len(fails))
    sys.exit(1)
print("all consistency checks passed")
