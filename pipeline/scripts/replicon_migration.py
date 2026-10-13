#!/usr/bin/env python3
"""Three-state replicon summary + old->new migration for the P258 analysis set.

Reads the inventory already written by replicon_three_state.py and joins it to
the deletion profile (not to the directory names, which carry a version
suffix) and to assembly metadata to restrict to P258 (complete/chromosome).

Outputs results/qc/replicon_migration_P258.tsv and prints the summary.
"""
import csv, re, collections

norm = lambda a: re.sub(r"\.\d+$", "", (a or "").strip())
TAG2P = {"NC_003131": "pCD1", "NC_003132": "pPCP1", "NC_003134": "pMT1"}

inv = {norm(r["accession"]): r for r in
       csv.DictReader(open("results/qc/replicon_three_state.tsv"), delimiter="\t")}
prof = {norm(r["accession"]): (r.get("replicons_absent") or "").strip() for r in
        csv.DictReader(open("results/qc/deletion_profile_all.tsv"), delimiter="\t")}
meta = {}
for r in csv.DictReader(open("results/qc/assembly_metadata.tsv"), delimiter="\t"):
    meta[norm(r["accession"])] = r

P258 = [a for a in inv
        if a in prof
        and (meta.get(a, {}).get("assembly_level") or "").lower()
        .startswith(("complete", "chromosome"))]
print("P258 (in inventory and profile, complete/chromosome): %d" % len(P258))

print()
print("=== three-state counts in P258 ===")
c = collections.Counter()
for a in P258:
    for p in TAG2P.values():
        c[(p, inv[a][p + "_state"])] += 1
for p in TAG2P.values():
    tot = sum(c[(p, s)] for s in ("CARRIED", "LOST", "UNRESOLVED"))
    print("  %-6s CARRIED=%-4d LOST=%-4d UNRESOLVED=%-4d  (sum %d)"
          % (p, c[(p, "CARRIED")], c[(p, "LOST")], c[(p, "UNRESOLVED")], tot))

print()
print("=== old `replicons_absent` -> new three-state (P258) ===")
mig = collections.Counter()
for a in P258:
    oa = set(prof[a].split(";")) if prof[a] else set()
    for tag, p in TAG2P.items():
        mig[("absent" if tag in oa else "present", inv[a][p + "_state"])] += 1
for k in sorted(mig):
    print("  OLD_%-8s -> %-11s %5d" % (k[0], k[1], mig[k]))

print()
print("=== genomes with NO plasmid sequence submitted (all three UNRESOLVED) ===")
unres = [a for a in P258 if inv[a]["pCD1_state"] == "UNRESOLVED"
         and inv[a]["pPCP1_state"] == "UNRESOLVED" and inv[a]["pMT1_state"] == "UNRESOLVED"]
print("  count: %d / %d" % (len(unres), len(P258)))
print("  of these, old profile claimed a plasmid loss: %d"
      % sum(1 for a in unres if prof[a]))
print("  e.g.: " + ", ".join(sorted(unres)[:8]))

print()
print("=== genuine LOST in P258 (submitter packaged plasmids, this one absent) ===")
for tag, p in TAG2P.items():
    lost = [a for a in P258 if inv[a][p + "_state"] == "LOST"]
    print("  %-6s LOST n=%-3d  %s" % (p, len(lost), ", ".join(sorted(lost)[:12])))

with open("results/qc/replicon_migration_P258.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["accession", "n_seq", "n_plasmid_seq", "old_replicons_absent",
                "pCD1_state", "pPCP1_state", "pMT1_state"])
    for a in sorted(P258):
        r = inv[a]
        w.writerow([a, r["n_seq"], r["n_plasmid_seq"], prof[a],
                    r["pCD1_state"], r["pPCP1_state"], r["pMT1_state"]])
print()
print("WROTE results/qc/replicon_migration_P258.tsv (%d rows)" % len(P258))
