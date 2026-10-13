#!/usr/bin/env python3
"""Generate tables/S1..S8.tsv from the archived derived tables.

Each table is written as tab-separated text with a leading comment line recording
its source, so a reader can trace any cell back to a file in results/.
"""
import csv, json, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "tables")
os.makedirs(OUT, exist_ok=True)


def rd(rel, delim="\t"):
    with open(os.path.join(ROOT, rel), newline="") as fh:
        return list(csv.DictReader(fh, delimiter=delim))


def wr(name, header, rows, source):
    p = os.path.join(OUT, name)
    with open(p, "w", newline="") as fh:
        fh.write("# source: %s\n" % source)
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)
    print("  wrote tables/%-12s %5d rows" % (name, len(rows)))


print("generating supplementary tables")
# S1 software manifest
rows = [["snippy", "4.6.0", "--ctgs --mincov 10 --minfrac 0.9"],
        ["bakta", "1.9.4", "DB 5.1 full, MD5 f8823533b789dd315025fdcc46f1a8c1"],
        ["pyrodigal", "3.7.1", "-"], ["pyhmmer", "0.10.15", "pinned"],
        ["AMRFinderPlus", "4.2.7", "DB 2026-08-07.1"],
        ["NCBI Datasets CLI", "18.38.0", "taxon 632"],
        ["SciPy", "1.18.1", "statistics"],
        ["R", "4.5.3", "figure rendering: ggplot2 4.0.3, patchwork 1.3.2, ragg 1.5.2, ggrepel 0.9.8"],
        ["R", "4.4.3", "remaining R analysis"],
        ["MUMmer4/dnadiff", "-", "AlignedBases"],
        ["CheckM", "v1.2.4 (reported by NCBI)", "completeness, contamination"]]
wr("S1_software.tsv", ["tool", "version", "notes"], rows, "pipeline/config.yaml; env/ENV_VERSIONS.txt; env/R_SESSION_INFO.txt")

# S2 recurrent blocks + accessions
blocks = rd("results/qc/uncalled_blocks_all.tsv")
raw = [r for r in blocks if int(r["length"]) >= 10000]
coords = {}
for r in raw:
    coords.setdefault((r["replicon"], r["start"], r["length"]), 0)
    coords[(r["replicon"], r["start"], r["length"])] += 1
rows = [[k[0], k[1], k[2], v] for k, v in sorted(coords.items(), key=lambda x: -x[1])]
wr("S2_recurrent_blocks.tsv", ["replicon", "start_5kb_bin", "length_bp", "n_genomes"], rows,
   "results/qc/uncalled_blocks_all.tsv")
acc = [r["accession"] for r in rd("results/qc/tierB_candidates.tsv")]
wr("S2_accessions.tsv", ["accession"], [[a] for a in acc], "results/qc/tierB_candidates.tsv")

# S3 threshold sweep, Tier A (archived) and Tier B (recomputed from the matrix)
sa = rd("results/figures/fig4d_sweep.csv", ",")
rows = [[r["threshold"], r["retained"], r["max_diameter"], r["max_size"], r["dup_pairs_merged"]] for r in sa]
wr("S3_sweep_tierA.tsv", ["t", "retained", "max_diameter", "max_size", "same_strain_pairs_merged"], rows,
   "results/figures/fig4d_sweep.csv")
try:
    import numpy as np
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    D = np.load(os.path.join(ROOT, "results/dedup/snp_distance_raw_tierB.npy")).astype(np.int64)
    np.fill_diagonal(D, 0)
    ss = max(int(D.max()), 1)
    Z = linkage(squareform(D.astype(float), checks=False), method="complete")
    rows = []
    for t in [0, 3, 5, 10, 15, 18, 19, 20, 25, 30, 50, 75, 100]:
        lab = fcluster(Z, t=t, criterion="distance")
        sizes = np.bincount(lab)[1:]
        rows.append([t, int(len(sizes)), int(sizes.max()), "grid" if t in (0, 3, 5, 10, 15, 20, 25, 30, 50, 75, 100) else "off-grid"])
    wr("S3_sweep_tierB.tsv", ["t", "retained", "max_size", "grid_point"], rows,
       "results/dedup/snp_distance_raw_tierB.npy (complete linkage)")
except Exception as e:                                        # noqa: BLE001
    print("  (S3 tierB skipped: %r)" % e)

# S4 plasmid three-state for the closed genomes
prof = {r["accession"]: r for r in rd("results/qc/deletion_profile_all.tsv")}
meta = {}
for r in rd("results/qc/assembly_metadata.tsv"):
    meta[r["accession"]] = r
    if r.get("paired_accession"):
        meta.setdefault(r["paired_accession"], r)
ts = rd("results/qc/replicon_three_state_seqreport.tsv")
rows = []
for r in ts:
    a = r["accession"]
    m = meta.get(a) or meta.get(a.split(".")[0])
    if not m:
        continue
    if not (m.get("assembly_level") or "").lower().startswith(("complete", "chromosome")):
        continue
    rows.append([a, r["n_plasmid_seq"], r["pCD1_state"], r["pPCP1_state"], r["pMT1_state"]])
wr("S4_plasmid_three_state.tsv", ["accession", "n_plasmid_sequences", "pCD1", "pPCP1", "pMT1"], rows,
   "results/qc/replicon_three_state_seqreport.tsv + assembly_metadata.tsv")

# S5 baseline layer per genome
lab = {r["accession"] for r in rd("results/qc/lab_strain_inventory.tsv")}
rows = []
for a, r in sorted(prof.items()):
    m = meta.get(a, {})
    if not (m.get("assembly_level") or "").lower().startswith(("complete", "chromosome")):
        continue
    rows.append([a, m.get("assembly_level", ""), "yes" if a in lab else "no",
                 r.get("replicons_absent", ""), r.get("big_blocks_5kb", ""),
                 r.get("n_big", ""), m.get("strain", "")])
wr("S5_baseline_layer.tsv", ["accession", "assembly_level", "laboratory_named",
                             "replicons_absent", "big_blocks_5kb", "n_big_blocks", "strain"], rows,
   "results/qc/deletion_profile_all.tsv + assembly_metadata.tsv + lab_strain_inventory.tsv")

# S6 pseudogene entanglement
pe = rd("results/qc/pseudogene_entanglement.tsv")
cols = list(pe[0].keys())
wr("S6_pseudogene_entanglement.tsv", cols, [[r[c] for c in cols] for r in pe],
   "results/qc/pseudogene_entanglement.tsv")

# S7 interpretation guidance (authored text)
rows = [
 ["Laboratory-derived vs naturally circulating",
  "SNP distance to a reference; a small distance is read as the same as a known strain",
  "SNP distance is invalid for this question; a 102.6 kb deletion contributes zero. Check the deletion profile independently: >=10 kb uncalled runs plus the replicon three-state"],
 ["Including an archive genome in a phylogeny or PAV analysis",
  "Dereplicate by strain field, or by an ANI or SNP threshold",
  "Dereplication must carry a deletion-profile concordance constraint, and that constraint must distinguish absent from not submitted"],
 ["Reporting 'X % of strains in a region lack F1'",
  "Compute a frequency from public genomes",
  "That is an archive frequency: an upper bound, reportable only as an interval and only for closed genomes"],
 ["Deciding whether an absence is biological",
  "Reference-relative absence is scored as a loss",
  "Determine polarity first; the dataset's most frequent absence is AMBIGUOUS and 2 of 13 recurrent regions are gains of the reference lineage"],
]
wr("S7_interpretation.tsv", ["setting", "current practice", "the rule this paper gives"], rows,
   "Section 4.7")

# S8 second reference
def latest(path, key):
    v = None
    for line in open(os.path.join(ROOT, path)):
        line = line.strip()
        if line:
            d = json.loads(line)["numbers"]
            v = d.get(key, v)
    return v
c = latest("results/ref2/numbers.json", "constraints")
e1 = latest("results/ref2/numbers.json", "endpoint_E1_co92_pair")
rows = [
 ["Reference accession", "GCF_000009065.1", "GCF_000013825.1"],
 ["Biovar", "Orientalis", "Antiqua"],
 ["Assembly length (bp)", "4,829,855", "4,879,836"],
 ["CheckM completeness / contamination (%)", "98.69 / 0.86", "97.7 / 0.96"],
 ["Union variant sites (Tier A)", "12,926", "12,835"],
 ["Same-strain pairs merging at t = 15", "3 of 4", "3 of 4"],
 ["Same-strain pairs merging at t = 20", "4 of 4", "4 of 4"],
 ["Largest cluster at t = 20", "12", "12"],
 ["CO92 pair, raw distance", "13", str(e1.get("raw_distance")) if e1 else ""],
 ["CO92 pair, distance alone", "merged", "merged"],
 ["CO92 pair, conservative constraint", "separated", "separated"],
 ["Multi-member clusters, no constraint", "36", str(c["off"]["multi_member"]) if c else ""],
 ["Discordant clusters, resolved", "24 (21 / 3)", "25 (23 / 2)"],
 ["Discordant proportion", "66.7 %", "73.5 %"],
 ["Genomes retained, off -> conservative", "175 -> 209", "175 -> 208"],
]
wr("S8_second_reference.tsv", ["quantity", "CO92 reference", "Antiqua reference"], rows,
   "results/ref2/numbers.json and verify.json; results/qc/variant_freq_strata.tsv")
print("done")
