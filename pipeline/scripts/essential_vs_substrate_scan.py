#!/usr/bin/env python3
"""
Which medium components are auxotrophic REQUIREMENTS, and which are
substrates the strain could get elsewhere?

WHY A GROWTH CURVE CANNOT ANSWER THIS
Your diagnostic showed:

    base + phe            -> 0
    base + phe + cit      -> 0.0885
    base + phe + gly      -> 0.0456
    base + phe + cit + gly-> 0.2833   (full)

Every one of those components "makes it grow when added", so growth alone
cannot separate "fills a biosynthetic gap" from "is a usable substrate".
Ask a different question instead:

    CAN ANY OTHER SINGLE COMPONENT SUBSTITUTE FOR IT?

    nothing substitutes  -> auxotrophic REQUIREMENT. Must stay open in the
                            oligotrophic medium, or the wild type cannot
                            grow and every knockout reads as lethal.
    something substitutes-> SUBSTRATE. Belongs in the carbon-source scan as
                            a variable.

This matters because iPC815 covers 815 of ~4,000 genes and the published
model explicitly identifies gaps in several biosynthetic pathways and
offers them as explanations for this strain's known nutritional
requirements. A gap in the model is indistinguishable from a real
auxotrophy by growth rate alone, and the distinction decides what the
oligotrophic medium even means.

WHAT THE OLIGOTROPHIC MEDIUM IS FOR
Not a single growth number. The readout is the BREADTH of usable carbon
sources per lineage -- how many of N sources support growth. That is what
"when was environmental persistence lost" actually asks, and a lineage that
uses 20 sources is more versatile than one that uses 5. A medium on which
the wild type scores zero carries no information at all, because every
knockout scores zero too.

Usage:
    python3 scripts/essential_vs_substrate_scan.py \
        --sbml data/models/iPC815.xml \
        --out results/fba/component_roles.tsv
"""
import argparse
import csv
import sys
from pathlib import Path

GROWTH_EPS = 1e-9


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sbml", required=True)
    ap.add_argument("--objective", default=None,
                    help="biomass reaction; default = the model's own")
    ap.add_argument("--max-substitutes", type=int, default=40,
                    help="how many candidate substitutes to try per component")
    ap.add_argument("--out", default="results/fba/component_roles.tsv")
    a = ap.parse_args()

    try:
        import cobra
    except ImportError:
        sys.exit("cobra not importable; see envs/cobra.yaml")

    model = cobra.io.read_sbml_model(a.sbml)
    if a.objective:
        model.objective = a.objective

    default_medium = dict(model.medium)
    full = float(model.optimize().objective_value or 0.0)
    if full <= GROWTH_EPS:
        sys.exit("FATAL: the model does not grow on its own default medium")

    all_exchanges = [r.id for r in model.exchanges]
    met_name = {}
    for r in model.exchanges:
        for m in r.metabolites:
            met_name[r.id] = m.name or m.id
            break

    print("=" * 84)
    print("Essential requirement versus substitutable substrate")
    print("=" * 84)
    print(f"model            : {model.id}")
    print(f"default medium   : {len(default_medium)} components")
    print(f"full growth      : {full:.6f}")
    print()

    # ---- step 1: which default components are individually necessary? ----
    necessary = []
    print("[1] drop one component at a time from the default medium")
    for rid in sorted(default_medium):
        with model:
            med = dict(default_medium)
            med.pop(rid, None)
            try:
                model.medium = med
                g = float(model.optimize().objective_value or 0.0)
            except Exception:                                 # noqa: BLE001
                g = 0.0
        frac = g / full
        tag = "NECESSARY" if g <= GROWTH_EPS else (
            "limiting" if frac < 0.5 else "dispensable")
        if g <= GROWTH_EPS:
            necessary.append(rid)
        print(f"    {rid:<18}{met_name.get(rid, ''):<34}"
              f"{g:>10.6f}  {frac:>6.1%}  {tag}")

    print(f"\n    individually necessary: {len(necessary)} "
          f"-> {[met_name.get(r, r) for r in necessary]}")

    # ---- step 2: for each necessary one, can anything substitute? --------
    print("\n[2] substitution scan -- can any single other exchange replace it?")
    rows = []
    candidates = [r for r in all_exchanges if r not in default_medium]
    for rid in necessary:
        base = dict(default_medium)
        base.pop(rid, None)
        subs = []
        for cand in candidates[: a.max_substitutes]:
            with model:
                med = dict(base)
                med[cand] = 10.0
                try:
                    model.medium = med
                    g = float(model.optimize().objective_value or 0.0)
                except Exception:                             # noqa: BLE001
                    g = 0.0
            if g > GROWTH_EPS:
                subs.append((cand, met_name.get(cand, cand), g))

        role = "AUXOTROPHIC_REQUIREMENT" if not subs else "SUBSTRATE"
        rows.append({
            "exchange": rid,
            "metabolite": met_name.get(rid, rid),
            "role": role,
            "n_substitutes": len(subs),
            "substitutes": ",".join(f"{m}({g:.3f})" for _, m, g in subs[:5]),
        })
        print(f"    {rid:<18}{met_name.get(rid, ''):<30}{role:<24}"
              f"{len(subs)} substitute(s)")
        for _, m, g in subs[:3]:
            print(f"        <- {m}  growth {g:.6f}")

    # dispensable defaults are substrates by definition
    for rid in sorted(default_medium):
        if rid in necessary:
            continue
        rows.append({"exchange": rid, "metabolite": met_name.get(rid, rid),
                     "role": "SUBSTRATE", "n_substitutes": "",
                     "substitutes": "dispensable in the default medium"})

    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    req = [r for r in rows if r["role"] == "AUXOTROPHIC_REQUIREMENT"]
    sub = [r for r in rows if r["role"] == "SUBSTRATE"]

    print()
    print("=" * 84)
    print("OLIGOTROPHIC MEDIUM DEFINITION")
    print("=" * 84)
    print("always open (auxotrophic requirements; without these the wild type")
    print("cannot grow and every knockout reads as lethal):")
    for r in req:
        print(f"    {r['exchange']:<18}{r['metabolite']}")
    print()
    print(f"carbon-source variables (closed in the base, reopened one at a "
          f"time): {len(sub)}")
    print()
    print("Expected wild-type growth on base + one carbon source: positive")
    print("and sub-maximal. The measured 0.0885 (31% of full) is a good")
    print("baseline -- headroom in both directions. Zero is useless, because")
    print("every knockout would also be zero.")
    print()
    print("The Fig.5 oligotrophic readout is the COUNT of carbon sources a")
    print("lineage can use, not a single growth rate.")
    print(f"written: {a.out}")


if __name__ == "__main__":
    main()
