#!/usr/bin/env python3
"""
Merge Tier A assembly PAV, Tier C coarse PAV and pseudogene states into the
character matrix PastML reconstructs on.

DECISION THAT MUST NOT BE LOST:
Three-state output (present / absent / unknown). Tier C unknowns must stay
unknown, not be imputed to absent.

WHY THE PANAROO CSV IS AN INPUT
-------------------------------
The Rtab gives presence/absence per orthogroup but carries no locus tags, and
a pseudogene call is a locus tag in one strain. gene_presence_absence.csv is
where Panaroo records which locus tags belong to which orthogroup, so it is
the only bridge from the pseudogene table to the orthogroup axis. Without it
pseudogene calls cannot be attributed to any column and would be dropped
silently.

STATE CODING
------------
    1  present and functional
    0  lost -- either absent from the assembly, or present but a pseudogene
    ?  unknown -- no usable observation (Tier C without breadth, or a strain
       absent from a source table)

Collapsing pseudogene into 0 is deliberate: this study measures loss of
function, and a pseudogene is a loss of function even though the sequence is
still there. Deletions and pseudogenisation stay distinguishable in
pseudo_calls_long.tsv for anyone who needs the split.

Tier A and Tier C are merged into one matrix because ancestral
reconstruction needs a single character matrix over one tree. Tier A rows are
dense (assembly PAV); Tier C rows come from read breadth and are largely '?'.
That is the honest representation -- imputing Tier C absences from absence of
evidence would manufacture losses along the deepest branches of the tree.
"""
import csv
import sys
from pathlib import Path

snakemake = globals().get("snakemake")  # noqa: F821

PRESENT, ABSENT, UNKNOWN = "1", "0", "?"


def read_rtab(path):
    """Panaroo Rtab: first column = orthogroup, remaining columns = samples."""
    p = Path(path)
    if not p.exists():
        raise SystemExit(f"FATAL: missing {p}")
    with p.open(encoding="utf-8", errors="replace") as fh:
        rows = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]
    if not rows:
        raise SystemExit(f"FATAL: {p} is empty")
    header = rows[0]
    samples = header[1:]
    table = {}
    for r in rows[1:]:
        if len(r) < 2:
            continue
        table[r[0]] = dict(zip(samples, r[1:]))
    return samples, table


def read_panaroo_locus_map(path):
    """gene_presence_absence.csv -> {(sample, locus_tag): orthogroup}."""
    mapping = {}
    p = Path(path)
    if not p.exists():
        return mapping, set()
    with p.open(encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.reader(fh)
        try:
            header = next(reader)
        except StopIteration:
            return mapping, set()
        try:
            i_group = header.index("Gene") if "Gene" in header else 0
        except ValueError:
            i_group = 0
        sample_cols = {}
        for i, name in enumerate(header):
            if name in ("Gene", "Non-unique Gene name", "Annotation",
                        "No. isolates", "No. sequences", "Avg sequences per isolate",
                        "Genome Fragment", "Order within Fragment", "Accessory Fragment",
                        "Accessory Order with Fragment", "QC", "Min group size nuc",
                        "Max group size nuc", "Avg group size nuc"):
                continue
            sample_cols[i] = name
        for row in reader:
            if len(row) <= i_group:
                continue
            group = row[i_group]
            for i, sample in sample_cols.items():
                if i >= len(row):
                    continue
                for tag in row[i].split(";"):
                    tag = tag.strip()
                    if tag and tag != "":
                        mapping[(sample, tag)] = group
    return mapping, set(sample_cols.values())


def read_pseudo_calls(path):
    """pseudo_calls_long.tsv -> {(accession, locus_tag)} for agreed calls."""
    calls = set()
    p = Path(path)
    if not p.exists():
        return calls
    with p.open(encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if str(r.get("agreement", "")).strip() in ("1", "True", "true"):
                acc = r.get("accession", "")
                tag = (r.get("locus_tag") or "").strip()
                if acc and tag:
                    calls.add((acc, tag))
    return calls


def read_coarse_pav(path):
    """Tier C breadth-based states -> {taxon: {group: state}}.

    Expected columns: accession, group, state (present/absent/unknown).
    Absent file -> Tier C simply does not appear.
    """
    out = {}
    p = Path(path)
    if not p.exists():
        return out
    with p.open(encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            tax = r.get("accession") or r.get("taxon") or r.get("strain")
            grp = r.get("group") or r.get("orthogroup") or r.get("gene")
            st = (r.get("state") or "").strip().lower()
            if not tax or not grp:
                continue
            if st in ("present", "1", "p"):
                v = PRESENT
            elif st in ("absent", "0", "a"):
                v = ABSENT
            else:
                v = UNKNOWN
            out.setdefault(tax, {})[grp] = v
    return out


def main():
    samples, rtab = read_rtab(snakemake.input.rtab)
    # gene_presence_absence.csv is where Panaroo records which locus tags belong
    # to which orthogroup. It is the only bridge from a locus-level pseudogene
    # call to the orthogroup axis; without it pseudogene calls would be silently
    # dropped. Declared as an explicit input rather than derived from the Rtab
    # path -- a string substitution that silently reads the wrong file is
    # exactly the failure mode this pipeline is built to avoid.
    locus_map, csv_samples = read_panaroo_locus_map(snakemake.input.csv)
    pseudo = read_pseudo_calls(snakemake.input.pseudo)
    coarse = read_coarse_pav(snakemake.input.coarse)

    groups = sorted(rtab.keys())
    group_idx = {g: i for i, g in enumerate(groups)}

    # ---- Tier A rows -----------------------------------------------------
    states = {s: [UNKNOWN] * len(groups) for s in samples}
    for s in samples:
        row = states[s]
        for g in groups:
            v = rtab[g].get(s, "")
            row[group_idx[g]] = PRESENT if str(v).strip() in ("1", "1.0") else ABSENT

    # pseudogene collapse: present locus that Panaroo grouped -> set to 0
    n_pseudo_mapped = 0
    n_pseudo_unmapped = 0
    for (acc, tag) in pseudo:
        if acc not in states:
            n_pseudo_unmapped += 1
            continue
        grp = locus_map.get((acc, tag))
        if grp is None or grp not in group_idx:
            n_pseudo_unmapped += 1
            continue
        states[acc][group_idx[grp]] = ABSENT
        n_pseudo_mapped += 1

    # ---- Tier C rows -----------------------------------------------------
    tier_c_taxa = sorted(coarse.keys())
    for tax in tier_c_taxa:
        row = [UNKNOWN] * len(groups)
        for grp, v in coarse[tax].items():
            if grp in group_idx:
                row[group_idx[grp]] = v
        states[tax] = row

    # ---- emit ------------------------------------------------------------
    out = Path(snakemake.output.tsv)
    out.parent.mkdir(parents=True, exist_ok=True)
    taxa = sorted(states.keys())
    nA = len(samples)
    nC = len(tier_c_taxa)
    # Columns carry the Panaroo orthogroup name itself, so downstream
    # enrichment can go straight to KEGG/COG without an index indirection.
    # A sanitised alias is written alongside for tools that dislike the name.
    aliases = {g: f"g{i}" for i, g in enumerate(groups)}
    with out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["taxon"] + groups)
        for t in taxa:
            w.writerow([t] + states[t])

    # group name lookup for downstream, so any gene id can be resolved
    with out.with_name("combined_group_index.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["orthogroup", "alias"])
        for g in groups:
            w.writerow([g, aliases[g]])

    total = sum(len(v) for v in states.values())
    n_unknown = sum(1 for v in states.values() for x in v if x == UNKNOWN)
    n_absent = sum(1 for v in states.values() for x in v if x == ABSENT)
    n_present = total - n_unknown - n_absent

    # Two distinct causes produce a low join rate and they need different
    # fixes, so the message must separate them. --remove-invalid-genes in
    # panaroo's extra args strips internal-stop genes (pseudogenes), which
    # leaves the call present in pseudo_calls_long but absent from the
    # orthogroup table -- that is NOT a namespace problem even though the
    # symptom looks the same.
    if n_pseudo_mapped + n_pseudo_unmapped:
        rate = n_pseudo_mapped / (n_pseudo_mapped + n_pseudo_unmapped)
        if rate < 0.95:
            unmapped_tags = sorted({t for _a, t in pseudo
                                    if _a in states and (locus_map.get((_a, t)) is None)})
            n_absent_from_panaroo = len(unmapped_tags)
            n_contig_failed = n_pseudo_unmapped - n_absent_from_panaroo
            raise SystemExit(
                f"FATAL: pseudogene->orthogroup join rate {rate:.1%} "
                f"({n_pseudo_mapped}/{n_pseudo_mapped + n_pseudo_unmapped}).\n"
                f"  calls whose locus_tag is absent from the panaroo csv: "
                f"{n_absent_from_panaroo}\n"
                f"    -> suspect --remove-invalid-genes stripping internal-stop "
                f"genes (config.pangenome.extra)\n"
                f"  calls unattributable for other reasons: {n_contig_failed}\n"
                f"    -> unmapped tag examples: {unmapped_tags[:3]}\n"
                f"    -> csv samples: {sorted(csv_samples)[:3]}"
            )

    msg = "\n".join([
        f"[combine_pav] orthogroups={len(groups)} taxa={len(taxa)} "
        f"(TierA={nA}, TierC={nC})",
        f"[combine_pav] cells={total} present={n_present} absent={n_absent} "
        f"unknown={n_unknown} ({100 * n_unknown / max(total, 1):.1f}%)",
        f"[combine_pav] pseudogene calls mapped={n_pseudo_mapped} "
        f"unmapped={n_pseudo_unmapped}",
        f"[combine_pav] wrote {out}",
    ])
    print(msg, file=sys.stderr)
    Path(snakemake.log[0]).write_text(msg + "\n")


if __name__ == "__main__":
    main()
