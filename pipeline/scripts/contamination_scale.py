#!/usr/bin/env python3
"""P0 (SY1011): quantify the plasmid-loss / no-plasmid-sequence conflation.

§4.10 counts plasmid "losses" from the `replicons_absent` field.  But a genome
whose assembly contains ONLY a chromosome (n_seq <= 1) never had a plasmid
sequence to lose, so a "pCD1 absent" call on it is an assembly artefact, not a
loss.  This script measures how much of each plasmid class is such genomes.

Outputs
  results/qc/contamination_scale.txt          full report
  results/qc/plasmid_claim_vs_sequence.tsv    per-strain table
  results/qc/nseq_per_strain.tsv              cached n_seq per accession

Discipline (handoff §5): every column is asserted to exist and be non-empty
before use; no `or <plausible default>` coercion is allowed anywhere.
"""
import csv, glob, os, re, collections, subprocess, sys
import os

ROOT = os.environ.get("PLAGUE_DECAY_ROOT",
                       os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.chdir(ROOT)
norm = lambda a: re.sub(r"\.\d+$", "", (a or "").strip())


def read(path):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return rows


# ---- load deletion profile, asserting the columns exist -------------------
prows = read("results/qc/deletion_profile_all.tsv")
assert prows, "deletion_profile_all.tsv is EMPTY"
missing = [c for c in ("accession", "replicons_absent") if c not in prows[0]]
assert not missing, "deletion_profile_all.tsv missing columns: %s" % missing
prof = {}
for r in prows:
    prof[norm(r["accession"])] = r
assert all(norm(k) == k and k for k in prof), "empty/normalisation clash in accession"

# ---- load metadata --------------------------------------------------------
mrows = read("results/qc/assembly_metadata.tsv")
assert mrows, "assembly_metadata.tsv is EMPTY"
assert "assembly_level" in mrows[0], "assembly_metadata.tsv missing assembly_level"
meta = {}
for r in mrows:
    meta.setdefault(norm(r["accession"]), r)

# ---- n_seq with cache + grep (C speed) ------------------------------------
CACHE = "results/qc/nseq_per_strain.tsv"
cache = {}
if os.path.exists(CACHE):
    for r in csv.DictReader(open(CACHE), delimiter="\t"):
        try:
            cache[norm(r["accession"])] = int(r["n_seq"])
        except (KeyError, ValueError):
            pass


def nseq(k):
    if k in cache:
        return cache[k]
    g = sorted(glob.glob("data/genomes/%s*/*.fna" % k))
    if not g:
        cache[k] = -1
        return -1
    n = 0
    for f in g:
        out = subprocess.run(["grep", "-c", "^>", f], capture_output=True, text=True)
        if out.returncode not in (0, 1):
            raise SystemExit("grep failed rc=%d on %s" % (out.returncode, f))
        n += int(out.stdout.strip() or 0)
    cache[k] = n
    return n


def lv(k):
    return (meta.get(k, {}).get("assembly_level") or "")


RES = []


def run(name, S):
    c = collections.Counter()
    miss = 0
    for k in S:
        n = nseq(k)
        if n < 0:
            miss += 1
            RES.append((k, lv(k), n, prof[k].get("replicons_absent", ""),
                        "MISSING_FASTA"))
            continue
        ab = bool((prof[k].get("replicons_absent") or "").strip())
        if ab and n <= 1:
            cat = "claims plasmid loss, but n_seq<=1"
        elif ab and n > 1:
            cat = "claims plasmid loss, n_seq>1"
        elif not ab and n <= 1:
            cat = "no plasmid loss claimed, n_seq<=1"
        else:
            cat = "no plasmid loss, n_seq>1"
        c[cat] += 1
        RES.append((k, lv(k), n, prof[k].get("replicons_absent", ""), cat))
    denom = len(S) - miss
    print("  %-8s n=%-5d no_fasta=%-4d" % (name, len(S), miss))
    for k, v in c.most_common():
        print("       %-40s %5d  (%.1f%%)" % (k, v, 100.0 * v / denom if denom else 0))


print("=" * 96)
print("CONFLATION TEST: does replicons_absent track plasmid SEQUENCES or plasmid CLAIMS?")
print("=" * 96)
run("P3441", sorted(prof))
P258 = sorted(k for k in prof if lv(k).lower().startswith(("complete", "chromosome")))
assert P258, "P258 is empty -- assembly_level values do not start with complete/chromosome"
run("P258", P258)

print()
REPS = (("NC_003134", "pMT1"), ("NC_003131", "pCD1"), ("NC_003132", "pPCP1"))
for rep, nm in REPS:
    cls = [k for k in P258 if rep in (prof[k].get("replicons_absent") or "")]
    n1 = sum(1 for k in cls if (nseq(k) >= 0 and nseq(k) <= 1))
    print("  P258 %-6s class n=%-3d   of which n_seq<=1: %d (%.0f%%)"
          % (nm, len(cls), n1, (100.0 * n1 / len(cls)) if cls else 0))

print()
no = [k for k in P258 if not (prof[k].get("replicons_absent") or "").strip()]
n1 = sum(1 for k in no if (nseq(k) >= 0 and nseq(k) <= 1))
print("  reverse: genomes WITHOUT replicons_absent but n_seq<=1 : %d / %d (%.0f%%)"
      % (n1, len(no), (100.0 * n1 / len(no)) if no else 0))

# by assembly level, all three plasmids claimed absent
print()
print("  P258 genomes claiming ALL THREE plasmids absent, by level and n_seq:")
allc = [k for k in P258 if all(rep in (prof[k].get("replicons_absent") or "")
                              for rep, _ in REPS)]
by = collections.Counter((lv(k), "n_seq<=1" if (nseq(k) >= 0 and nseq(k) <= 1) else "n_seq>1")
                         for k in allc)
for k, v in sorted(by.items()):
    print("     %-40s %5d" % (str(k), v))
print("     total all-three-absent in P258: %d" % len(allc))

# ---- write tables ---------------------------------------------------------
with open(CACHE, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["accession", "n_seq"])
    for k in sorted(cache):
        w.writerow([k, cache[k]])

with open("results/qc/plasmid_claim_vs_sequence.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["accession", "assembly_level", "n_seq", "replicons_absent", "category"])
    for row in sorted(RES):
        w.writerow(row)

print()
print("WROTE results/qc/plasmid_claim_vs_sequence.tsv  (%d rows)" % len(RES))
print("WROTE results/qc/nseq_per_strain.tsv           (%d rows)" % len(cache))
