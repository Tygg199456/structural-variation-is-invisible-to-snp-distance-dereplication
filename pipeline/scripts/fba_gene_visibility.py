#!/usr/bin/env python3
"""
Which of iPC815's 815 genes can FBA actually see, and do the three media
discriminate at all?

WHY RUN THIS NOW, BEFORE ANY STRAIN DATA EXISTS
strain_fba.py is blocked on the PAV matrix, which is blocked on panaroo,
which is blocked on the annotation. But this analysis needs no strain data:
it is a property of the model alone, and it answers three things that decide
how Fig.5 can be written.

1. IT TURNS THE 20% COVERAGE CAVEAT INTO A NUMBER
   "iPC815 covers 815 of ~4,000 genes" is the stated limitation. The real
   resolution of Fig.5 is smaller still: only the genes whose deletion
   changes predicted growth are visible at all. Everything else is a gene
   the model carries but cannot report on. That count is Fig.5's actual
   discriminating power and belongs in Methods.

2. IT PRE-SORTS THE PAV MATRIX
   When the matrix arrives, every lost gene can be split immediately into
   "FBA can speak to this" and "FBA is silent here" -- rather than running
   the whole analysis and then discovering which conclusions were empty.

3. THE ONE THAT COULD INVALIDATE Fig.5's MAIN AXIS
   If the essential-gene sets of mammalian blood, flea midgut and
   oligotrophic medium are nearly identical, the three columns of Fig.5
   carry almost no independent information and the mammal-versus-flea
   contrast -- the figure's whole point -- does not exist in this model.
   That is knowable today, with no strain data, and it should be known
   before 20 hours of annotation finishes rather than after the figure is
   drafted.

Method: single-gene deletion (GPR rules propagate; reactions are not
knocked out directly), per medium, with the temperature-matched biomass
objective bound for each. Genes are classed by the fraction of wild-type
growth retained.

Usage:
    python3 scripts/fba_gene_visibility.py \
        --sbml data/models/iPC815.xml \
        --media results/fba/media.json \
        --out results/fba/gene_visibility.tsv

Expected magnitude: 815 genes x 3 media. Expect a majority of genes to be
silent (no growth change) in every medium -- that is the point. If ZERO
genes are essential in some medium, that medium is not constraining and
needs rebuilding.
"""
import argparse
import csv
import json
import sys
from pathlib import Path

LETHAL = 0.01      # fraction of wild-type growth below which we call lethal
IMPAIRED = 0.95    # below this, growth is measurably reduced


def classify(frac):
    if frac is None:
        return "error"
    if frac <= LETHAL:
        return "essential"
    if frac < IMPAIRED:
        return "impaired"
    return "silent"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sbml", required=True)
    ap.add_argument("--media", required=True,
                    help="media.json from define_media.py")
    ap.add_argument("--out", default="results/fba/gene_visibility.tsv")
    a = ap.parse_args()

    try:
        import cobra
    except ImportError:
        sys.exit("cobra not importable; see envs/cobra.yaml")

    model = cobra.io.read_sbml_model(a.sbml)
    spec = json.loads(Path(a.media).read_text())
    media = spec["media"]

    print("=" * 86)
    print("iPC815 gene visibility -- what can FBA see, and do the media differ?")
    print("=" * 86)
    print(f"model : {model.id}   genes {len(model.genes)}   "
          f"reactions {len(model.reactions)}")
    print(f"media : {list(media)}")
    print()

    gene_ids = [g.id for g in model.genes]
    results = {}        # medium -> gene -> (growth, frac, class)
    wt = {}

    for name, m in media.items():
        obj = m.get("objective")
        medium = {k: v for k, v in m["medium"].items() if v > 0}
        # the oligotrophic base is carbon-free by design; open one source so
        # the knockout scan has a non-zero baseline to measure against
        if name == "oligotrophic":
            trials = m.get("single_carbon_trials") or []
            usable = [t for t in trials if t.get("growth", 0) > 1e-9]
            if usable:
                medium[usable[0]["exchange"]] = 10.0
                print(f"  {name}: opened {usable[0]['exchange']} "
                      f"({usable[0]['metabolite']}) for the scan baseline")

        with model:
            if obj:
                model.objective = obj
            try:
                model.medium = medium
            except Exception as exc:                        # noqa: BLE001
                print(f"  {name}: medium rejected ({exc}) -- skipped")
                continue
            base = float(model.optimize().objective_value or 0.0)
            wt[name] = base
            if base <= 1e-9:
                print(f"  {name}: wild-type growth is zero -- every knockout "
                      f"would read lethal. Skipped; fix the medium.")
                continue

            per_gene = {}
            for gid in gene_ids:
                gene = model.genes.get_by_id(gid)
                try:
                    with model:
                        gene.knock_out()
                        g = float(model.optimize().objective_value or 0.0)
                except Exception:                           # noqa: BLE001
                    g = None
                frac = None if g is None else (g / base)
                per_gene[gid] = (g, frac, classify(frac))
            results[name] = per_gene

    if not results:
        sys.exit("FATAL: no medium produced a usable baseline")

    print()
    print("wild-type growth per medium")
    for name in results:
        print(f"  {name:<18}{wt[name]:.6f}   objective "
              f"{media[name].get('objective')}")

    # ---- per-medium class counts ----------------------------------------
    print()
    print("gene classes per medium")
    print(f"  {'medium':<18}{'essential':>11}{'impaired':>10}{'silent':>9}"
          f"{'visible':>9}{'visible%':>10}")
    for name, per in results.items():
        ess = sum(1 for v in per.values() if v[2] == "essential")
        imp = sum(1 for v in per.values() if v[2] == "impaired")
        sil = sum(1 for v in per.values() if v[2] == "silent")
        vis = ess + imp
        print(f"  {name:<18}{ess:>11}{imp:>10}{sil:>9}{vis:>9}"
              f"{100 * vis / len(gene_ids):>9.1f}%")

    # ---- union of visibility: the real resolution of Fig.5 --------------
    visible_any = {gid for gid in gene_ids
                   if any(results[n][gid][2] in ("essential", "impaired")
                          for n in results)}
    print()
    print(f"genes visible in AT LEAST ONE medium : {len(visible_any)} / "
          f"{len(gene_ids)} ({100 * len(visible_any) / len(gene_ids):.1f}% "
          f"of the model)")
    print(f"as a fraction of the ~4,000-gene genome: "
          f"{100 * len(visible_any) / 4000:.1f}%")
    print("  -> THIS is Fig.5's actual discriminating power. Quote it in")
    print("     Methods instead of the 815/4000 figure, which overstates it.")

    # ---- do the media discriminate? -------------------------------------
    print()
    print("DO THE MEDIA DISCRIMINATE? pairwise Jaccard of essential sets")
    ess_sets = {n: {gid for gid, v in per.items() if v[2] == "essential"}
                for n, per in results.items()}
    names = list(ess_sets)
    worrying = []
    for i, n1 in enumerate(names):
        for n2 in names[i + 1:]:
            a_, b_ = ess_sets[n1], ess_sets[n2]
            union = a_ | b_
            j = len(a_ & b_) / len(union) if union else 1.0
            only1, only2 = len(a_ - b_), len(b_ - a_)
            print(f"  {n1} vs {n2}: Jaccard {j:.3f}   "
                  f"only {n1}: {only1}   only {n2}: {only2}")
            if j > 0.95:
                worrying.append((n1, n2, j))

    print()
    if worrying:
        print("WARNING -- Fig.5's main axis may not exist in this model:")
        for n1, n2, j in worrying:
            print(f"  {n1} and {n2} share {j:.1%} of their essential genes.")
        print("  Nearly identical essential sets mean those columns of Fig.5")
        print("  carry little independent information. Before drafting the")
        print("  figure, either differentiate the media further (tighter")
        print("  exchange bounds, the temperature-matched BOF actually")
        print("  bound) or drop the contrast and present a single condition")
        print("  with the carbon-breadth scan as the environmental axis.")
    else:
        print("The media give distinct essential-gene sets, so the")
        print("mammal-versus-flea contrast has a basis in this model.")

    # ---- write ----------------------------------------------------------
    rows = []
    for gid in gene_ids:
        gene = model.genes.get_by_id(gid)
        row = {"gene": gid, "name": getattr(gene, "name", "") or "",
               "n_reactions": len(gene.reactions)}
        for name in results:
            g, frac, cls = results[name][gid]
            row[f"{name}_growth"] = "" if g is None else round(g, 6)
            row[f"{name}_frac_wt"] = "" if frac is None else round(frac, 4)
            row[f"{name}_class"] = cls
        row["visible_any"] = "yes" if gid in visible_any else "no"
        rows.append(row)

    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print()
    print(f"written: {a.out}")
    print("When the PAV matrix arrives, join on this file first: every lost")
    print("gene splits into 'FBA can speak to this' and 'FBA is silent "
          "here'.")
    print("=" * 86)


if __name__ == "__main__":
    main()
