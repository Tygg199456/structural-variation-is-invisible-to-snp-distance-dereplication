#!/usr/bin/env python3
"""
Which biosynthetic capacities are conditionally essential, and under which
nutrient limitation?

WHY THIS REPLACES THE MAMMAL-VERSUS-FLEA AXIS
The measured Jaccard between the mammalian-blood and flea-midgut essential
gene sets is 0.974 -- they differ by 3 genes. The cause is structural, not a
tuning problem: comparing the two biomass objectives shows 47 of 62
metabolites identical and all 15 differences confined to membrane lipids
(clpn/pe/pg across chain lengths, plus LPS). Amino acids, nucleotides,
cofactors and ions are identical. Since FBA essentiality is determined by
whether a deletion blocks production of a required biomass component, and
both objectives require the same amino acids, nucleotides and cofactors in
the same amounts, the essential sets are forced to coincide. Three genes is
the structural ceiling.

Defining a genuinely different flea-midgut medium WOULD change essentiality
-- if the medium lacks a nutrient the blood medium supplies, the
biosynthetic pathway for it becomes conditionally essential. But the
quantitative composition of the flea midgut after a blood meal is not
characterised, and inventing exchange bounds would make the resulting
"difference" an artefact of invented numbers rather than a finding.

So the question changes:

    not   "how do lineages differ between mammal and flea"
    but   "which biosynthetic capacities are conditionally essential, and
           under which nutrient limitation would losing one matter"

That is a sensitivity analysis, not a claim to have modelled the vector.
Which nutrient limitations are plausible for the flea midgut belongs in the
Discussion, as discussion.

The output joins directly onto the parallel-loss gene list from Q2: for each
gene that was lost, this says whether the loss has a predicted consequence
and under what restriction.

Usage:
    python3 scripts/nutrient_dropout_scan.py \
        --sbml data/models/iPC815.xml \
        --media results/fba/media.json \
        --base mammalian_blood \
        --out results/fba/nutrient_dropout.tsv

Expected magnitude: one scan per nutrient class (roughly 6-10 classes) x 815
genes. Classes whose dropout makes the wild type unable to grow are reported
and skipped -- a lethal dropout carries no information, since every knockout
would also be lethal.
"""
import argparse
import csv
import json
import re
import sys
from pathlib import Path

LETHAL = 0.01
IMPAIRED = 0.95

# Nutrient classes, matched on metabolite names rather than reaction ids,
# because exchange ids differ between model builds.
CLASSES = {
    "amino_acids": re.compile(
        r"alanine|arginine|asparagin|aspart|cystein|glutam|glycine|histidin|"
        r"isoleucin|leucin|lysin|methionin|phenylalanin|prolin|serin|"
        r"threonin|tryptophan|tyrosin|valin", re.I),
    "nucleotide_precursors": re.compile(
        r"adenin|guanin|cytosin|thymin|uracil|inosine|xanthin|orotate|"
        r"nucleosid", re.I),
    "cofactors_vitamins": re.compile(
        r"thiamin|riboflavin|niacin|pantothen|pyridox|biotin|folate|"
        r"cobalamin|lipoate|quinone|heme|siroheme", re.I),
    "iron": re.compile(r"\bfe\s?\(?[23]\+?\)?|ferric|ferrous|\biron\b", re.I),
    "sulfur": re.compile(r"sulfate|sulfite|sulfide|thiosulfate", re.I),
    "nitrogen": re.compile(r"ammonium|ammonia|nitrate|nitrite|urea", re.I),
    "phosphate": re.compile(r"phosphate|\bpi\b", re.I),
    "carbon_sources": re.compile(
        r"glucose|glycerol|succinate|fumarate|malate|pyruvate|acetate|"
        r"citrate|galactose|mannose|maltose|ribose|xylose|gluconate|"
        r"lactate|arabinose|fructose|trehalose|sorbitol|mannitol", re.I),
}

# Components the scan must never close: the model cannot grow without them
# at all, so a dropout including them is uninformative.
NEVER_DROP = re.compile(r"\bo2\b|oxygen|\bh2o\b|water|\bh\+|proton", re.I)


def classify(frac):
    if frac is None:
        return "error"
    if frac <= LETHAL:
        return "essential"
    if frac < IMPAIRED:
        return "impaired"
    return "silent"


def exchange_label(model, rid):
    r = model.reactions.get_by_id(rid)
    for m in r.metabolites:
        return m.name or m.id
    return rid


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sbml", required=True)
    ap.add_argument("--media", required=True)
    ap.add_argument("--base", default="mammalian_blood")
    ap.add_argument("--out", default="results/fba/nutrient_dropout.tsv")
    a = ap.parse_args()

    try:
        import cobra
    except ImportError:
        sys.exit("cobra not importable; see envs/cobra.yaml")

    model = cobra.io.read_sbml_model(a.sbml)
    spec = json.loads(Path(a.media).read_text())
    if a.base not in spec["media"]:
        sys.exit(f"FATAL: '{a.base}' not in {list(spec['media'])}")
    base_spec = spec["media"][a.base]
    base_medium = {k: v for k, v in base_spec["medium"].items() if v > 0}
    objective = base_spec.get("objective")

    labels = {rid: exchange_label(model, rid) for rid in base_medium}

    print("=" * 84)
    print("Nutrient-dropout scan -- conditional essentiality")
    print("=" * 84)
    print(f"model      : {model.id}")
    print(f"base medium: {a.base}  ({len(base_medium)} components)")
    print(f"objective  : {objective}")

    # assign each base component to a class
    membership = {}
    for rid, lab in labels.items():
        if NEVER_DROP.search(f"{lab} {rid}"):
            membership[rid] = "never_drop"
            continue
        hit = None
        for cls, pat in CLASSES.items():
            if pat.search(f"{lab} {rid}"):
                hit = cls
                break
        membership[rid] = hit or "unclassified"

    print()
    print("base components by class")
    for cls in list(CLASSES) + ["never_drop", "unclassified"]:
        members = [labels[r] for r, c in membership.items() if c == cls]
        if members:
            print(f"  {cls:<24}{len(members):>3}  {members[:5]}")

    with model:
        if objective:
            model.objective = objective
        model.medium = base_medium
        wt_base = float(model.optimize().objective_value or 0.0)
    print()
    print(f"wild-type growth on the base medium: {wt_base:.6f}")
    if wt_base <= 1e-9:
        sys.exit("FATAL: the base medium does not support growth")

    gene_ids = [g.id for g in model.genes]
    scans = {}
    for cls in CLASSES:
        drop = [r for r, c in membership.items() if c == cls]
        if not drop:
            continue
        med = {k: v for k, v in base_medium.items() if k not in drop}
        with model:
            if objective:
                model.objective = objective
            try:
                model.medium = med
                wt = float(model.optimize().objective_value or 0.0)
            except Exception as exc:                        # noqa: BLE001
                print(f"  {cls}: medium rejected ({exc}) -- skipped")
                continue
            if wt <= 1e-9:
                print(f"  {cls:<24}wild type cannot grow without this class "
                      f"-> UNINFORMATIVE (every knockout would be lethal too)")
                scans[cls] = {"wt": 0.0, "genes": None,
                              "dropped": [labels[r] for r in drop]}
                continue
            per_gene = {}
            for gid in gene_ids:
                try:
                    with model:
                        model.genes.get_by_id(gid).knock_out()
                        g = float(model.optimize().objective_value or 0.0)
                except Exception:                            # noqa: BLE001
                    g = None
                per_gene[gid] = (g, None if g is None else g / wt)
            scans[cls] = {"wt": wt, "genes": per_gene,
                          "dropped": [labels[r] for r in drop]}
            n_ess = sum(1 for v in per_gene.values()
                        if classify(v[1]) == "essential")
            print(f"  {cls:<24}wt {wt:.6f} ({100 * wt / wt_base:.0f}% of base)"
                  f"   essential genes {n_ess}")

    usable = {c: s for c, s in scans.items() if s["genes"]}
    if not usable:
        sys.exit("FATAL: no dropout left the wild type viable")

    # baseline essentiality on the full base medium
    with model:
        if objective:
            model.objective = objective
        model.medium = base_medium
        base_ess = set()
        for gid in gene_ids:
            try:
                with model:
                    model.genes.get_by_id(gid).knock_out()
                    g = float(model.optimize().objective_value or 0.0)
            except Exception:                                # noqa: BLE001
                continue
            if g / wt_base <= LETHAL:
                base_ess.add(gid)

    print()
    print(f"essential on the full base medium: {len(base_ess)}")
    print()
    print("CONDITIONALLY ESSENTIAL -- essential only once a class is withheld")
    print(f"  {'dropout class':<24}{'newly essential':>16}")
    rows = []
    for cls, s in usable.items():
        newly = [gid for gid, v in s["genes"].items()
                 if classify(v[1]) == "essential" and gid not in base_ess]
        print(f"  {cls:<24}{len(newly):>16}")
        for gid in newly:
            gene = model.genes.get_by_id(gid)
            rows.append({
                "gene": gid, "name": getattr(gene, "name", "") or "",
                "dropout_class": cls,
                "withheld": ";".join(s["dropped"][:6]),
                "growth_without_gene": round(s["genes"][gid][0] or 0.0, 6),
                "frac_of_dropout_wt": round(s["genes"][gid][1] or 0.0, 4),
                "essential_on_full_base": "no",
            })
    for gid in sorted(base_ess):
        gene = model.genes.get_by_id(gid)
        rows.append({
            "gene": gid, "name": getattr(gene, "name", "") or "",
            "dropout_class": "unconditional",
            "withheld": "", "growth_without_gene": 0.0,
            "frac_of_dropout_wt": 0.0,
            "essential_on_full_base": "yes",
        })

    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    cond = len({r["gene"] for r in rows
                if r["essential_on_full_base"] == "no"})
    print()
    print("=" * 84)
    print(f"unconditionally essential : {len(base_ess)}")
    print(f"conditionally essential   : {cond}")
    print(f"total with a predicted phenotype: {len(base_ess) + cond} / "
          f"{len(gene_ids)} model genes "
          f"= {100 * (len(base_ess) + cond) / 4000:.1f}% of the ~4,000-gene "
          f"genome")
    print()
    print("FOR Fig.5")
    print("  Replace the mammal-versus-flea columns (Jaccard 0.974, a")
    print("  structural artefact of two biomass objectives that differ only")
    print("  in membrane lipids) with:")
    print("    axis 1  carbon-source breadth per lineage (wild type 45/76)")
    print("    axis 2  conditional essentiality x lost genes, from this table")
    print("  and report nutrient-rich versus oligotrophic as the contrast,")
    print("  which does discriminate (Jaccard 0.715-0.734).")
    print()
    print("FOR METHODS")
    print("  State the resolution honestly: the model has predictive")
    print("  resolution over this fraction of the genome, not the 20% that")
    print("  815/4000 implies.")
    print(f"written: {a.out}")
    print("=" * 84)


if __name__ == "__main__":
    main()
