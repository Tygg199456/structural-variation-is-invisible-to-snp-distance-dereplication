#!/usr/bin/env python3
"""Three-state replicon inventory from NCBI's structured sequence_report.jsonl.

Replaces FASTA-header text parsing (external review 2026-10-11): the identity
of a plasmid is taken from the submission's own machine-readable fields --
assignedMoleculeLocationType == "Plasmid", chrName/sequenceName, and length --
never from a header string.  Length is used only as a fallback for an unnamed
plasmid and only within 25% of a reference plasmid's length.

States
  CARRIED    plasmid sequence present and matched to this replicon
  LOST       plasmid sequence(s) submitted, but not this replicon
  UNRESOLVED no plasmid sequence at all (loss indistinguishable from
             non-submission)

Input : out/ncbi_dataset/data/<acc>/sequence_report.jsonl (datasets --include seq-report)
Output: results/qc/replicon_three_state_seqreport.tsv  and .txt summary
"""
import json, os, re, csv, glob, collections
import os

ROOT = os.environ.get("PLAGUE_DECAY_ROOT",
                       os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.chdir(ROOT)
SR = "/tmp/p258_sr/out/ncbi_dataset/data"
REFP = {"pCD1": 70305, "pPCP1": 9612, "pMT1": 96210}
NAME = {"pCD1": re.compile(r"pCD1", re.I),
        "pPCP1": re.compile(r"pPCP1", re.I),
        "pMT1": re.compile(r"pMT[-_ ]?1", re.I)}
TAG2P = {"NC_003131": "pCD1", "NC_003132": "pPCP1", "NC_003134": "pMT1"}
TOL = 0.25
norm = lambda a: re.sub(r"\.\d+$", "", (a or "").strip())

old = {norm(r["accession"]): (r.get("replicons_absent") or "").strip()
       for r in csv.DictReader(open("results/qc/deletion_profile_all.tsv"), delimiter="\t")}
meta = {}
for r in csv.DictReader(open("results/qc/assembly_metadata.tsv"), delimiter="\t"):
    meta[norm(r["accession"])] = r
    p = norm(r.get("paired_accession") or "")
    if p:
        meta.setdefault(p, r)

rows = []
for f in sorted(glob.glob(SR + "/*/sequence_report.jsonl")):
    recs = [json.loads(l) for l in open(f) if l.strip()]
    if not recs:
        continue
    acc = recs[0]["assemblyAccession"]
    plasmid = [r for r in recs
               if (r.get("assignedMoleculeLocationType") or "").lower() == "plasmid"]
    state, ev = {}, {}
    if not plasmid:
        for p in REFP:
            state[p] = "UNRESOLVED"
            ev[p] = "no plasmid in submitted sequence report"
    else:
        assigned = collections.defaultdict(list)
        for r in plasmid:
            nm = (r.get("chrName") or r.get("sequenceName") or "")
            L = int(r["length"])
            hit = None
            for p, pat in NAME.items():
                if pat.search(nm):
                    hit = p
                    break
            if hit is None:
                cands = [(abs(L - REFP[p]), p) for p in REFP if abs(L - REFP[p]) <= TOL * REFP[p]]
                if cands:
                    hit = min(cands)[1]
            if hit:
                assigned[hit].append((nm, L))
        for p in REFP:
            if assigned.get(p):
                state[p] = "CARRIED"
                ev[p] = ";".join("%s(%d)" % (h, l) for h, l in assigned[p])
            else:
                state[p] = "LOST"
                ev[p] = "plasmids submitted, not " + p
    rows.append((acc, len(recs), len(plasmid), sum(int(r["length"]) for r in recs),
                 state, ev, plasmid))

with open("results/qc/replicon_three_state_seqreport.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["accession", "n_seq", "n_plasmid_seq", "total_bp",
                "pCD1_state", "pPCP1_state", "pMT1_state", "old_replicons_absent", "evidence"])
    for acc, n_seq, npl, tot, state, ev, plasmid in rows:
        allev = "; ".join("%s(%s)" % (r.get("chrName") or r.get("sequenceName"),
                                      r.get("assignedMoleculeLocationType")) for r in plasmid) or "(none)"
        w.writerow([acc, n_seq, npl, tot, state["pCD1"], state["pPCP1"], state["pMT1"],
                    old.get(norm(acc), ""), allev])

c = collections.Counter()
for acc, n_seq, npl, tot, state, ev, plasmid in rows:
    for p in REFP:
        c[(p, state[p])] += 1
print("genomes (P258):", len(rows))
for p in REFP:
    print("  %-6s CARRIED=%-4d LOST=%-4d UNRESOLVED=%-4d"
          % (p, c[(p, "CARRIED")], c[(p, "LOST")], c[(p, "UNRESOLVED")]))

print()
print("=== old `replicons_absent` -> new three-state (P258) ===")
mig = collections.Counter()
for acc, n_seq, npl, tot, state, ev, plasmid in rows:
    oa = set(old.get(norm(acc), "").split(";")) if old.get(norm(acc), "") else set()
    if norm(acc) not in old:
        continue
    for tag, p in TAG2P.items():
        mig[("absent" if tag in oa else "present", state[p])] += 1
for k in sorted(mig):
    print("  OLD_%-8s -> %-11s %5d" % (k[0], k[1], mig[k]))

print()
unres = [a for a, n, npl, t, st, e, pl in rows
         if st["pCD1"] == "UNRESOLVED" and st["pPCP1"] == "UNRESOLVED" and st["pMT1"] == "UNRESOLVED"]
print("=== no plasmid submitted (all three UNRESOLVED): %d / %d = %.1f%%"
      % (len(unres), len(rows), 100.0 * len(unres) / len(rows)))
print("  " + ", ".join(sorted(unres)))

print()
print("=== genuine LOST (plasmids submitted, this one absent) ===")
for tag, p in TAG2P.items():
    lost = [a for a, n, npl, t, st, e, pl in rows if st[p] == "LOST"]
    print("  %-6s n=%-3d %s" % (p, len(lost), ", ".join(sorted(lost))))

print()
print("=== strains of interest ===")
for d in ("GCF_015337085.2", "GCF_015338205.1", "GCF_015337285.1", "GCF_000590535.2",
          "GCF_015337645.1", "GCF_000009065.1", "GCA_000006645.1", "GCF_000006645.1"):
    for row in rows:
        if row[0] == d:
            print("  %-18s n_seq=%d n_pl=%d pCD1=%-10s pPCP1=%-10s pMT1=%-10s" %
                  (d, row[1], row[2], row[4]["pCD1"], row[4]["pPCP1"], row[4]["pMT1"]))

with open("results/qc/replicon_three_state_seqreport.txt", "w") as fh:
    pass
print()
print("WROTE results/qc/replicon_three_state_seqreport.tsv (%d rows)" % len(rows))
