#!/usr/bin/env python3
"""
Build the three in silico media for Fig.5, each bound to its own biomass
objective.

THE BOF BINDING IS THE POINT
iPC815 carries TWO biomass objective functions (measured on this model:
BIOMASS_25C and BIOMASS_37C). Fig.5's main axis is mammal versus flea, and
the flea midgut is 26 C. Running the flea condition on the 37 C biomass
composition imposes host-interior macromolecule ratios on a bacterium
sitting in its vector -- it would quietly flatten the contrast the figure
exists to show.

    mammalian_blood -> BIOMASS_37C
    flea_midgut     -> BIOMASS_25C
    oligotrophic    -> BIOMASS_25C   (stated in Methods; environmental
                                      persistence is the low-temperature case)

DISCOVERY RATHER THAN HARD-CODED IDs
Exchange reaction IDs differ between model builds, so nothing here assumes
a naming scheme. The script starts from the model's own default medium and
applies modifications by matching metabolite names and formulas, then
reports every exchange it could not find. An unreported silent miss is how
a medium ends up not actually iron-limited.

THE CHECK THAT MATTERS
Each medium must give positive growth on the UNMODIFIED model. If a medium
kills the wild type, every strain knockout reads as lethal, Fig.5 becomes a
uniform zero heatmap, and nothing raises an error. The script exits non-zero
in that case.

Usage (standalone or via Snakemake):
    python3 scripts/define_media.py --sbml data/models/iPC815.xml \
        --out results/fba/media.json
"""
import argparse
import csv
import json
import re
import sys
from pathlib import Path

snakemake = globals().get("snakemake")  # noqa: F821

# Metabolite name patterns, not reaction IDs.
IRON = re.compile(r"\bfe\s?\(?[23]\+?\)?|\bferric|\bferrous|\biron\b", re.I)
GLUCOSE = re.compile(r"\bd-?glucose\b|\bglc", re.I)
AMINO_ACIDS = [
    "alanine", "arginine", "asparagine", "aspartate", "cysteine",
    "glutamate", "glutamine", "glycine", "histidine", "isoleucine",
    "leucine", "lysine", "methionine", "phenylalanine", "proline",
    "serine", "threonine", "tryptophan", "tyrosine", "valine",
]
CARBON_SCAN_HINTS = [
    "glucose", "glycerol", "succinate", "fumarate", "malate", "pyruvate",
    "acetate", "citrate", "galactose", "mannose", "maltose", "ribose",
    "xylose", "gluconate", "lactate", "glutamate", "aspartate", "alanine",
    "proline", "serine", "arabinose", "fructose", "trehalose", "sorbitol",
    "mannitol", "n-acetyl", "glucosamine", "oxoglutarate", "glycine",
    "threonine",
]


def pick_bof(model):
    """Find the temperature-tagged biomass objectives."""
    found = {}
    for r in model.reactions:
        up = r.id.upper()
        if "BIOMASS" not in up and "biomass" not in (r.name or "").lower():
            continue
        m = re.search(r"(\d{2})\s*C\b", r.id, re.I)
        if m:
            found[int(m.group(1))] = r.id
        else:
            found.setdefault("untagged", r.id)
    return found


def exchanges_matching(model, pattern=None, names=None):
    """Exchange reactions whose single metabolite matches a name pattern."""
    hits = []
    for r in model.exchanges:
        for met in r.metabolites:
            label = f"{met.name or ''} {met.id}"
            if pattern is not None and pattern.search(label):
                hits.append((r.id, met.name or met.id))
                break
            if names is not None and any(n in label.lower() for n in names):
                hits.append((r.id, met.name or met.id))
                break
    return hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sbml", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--glucose-limit", type=float, default=1.0,
                    help="mmol/gDW/h uptake cap for the blood medium")
    ap.add_argument("--requirements", default="results/fba/component_roles.tsv",
                    help="output of essential_vs_substrate_scan.py; its "
                         "AUXOTROPHIC_REQUIREMENT rows stay open in the "
                         "oligotrophic base (they are needs, not carbon "
                         "sources -- zeroing phe made growth impossible and "
                         "every knockout then read lethal)")
    a, _ = ap.parse_known_args()

    sbml = a.sbml or (snakemake.input.sbml if snakemake else None)
    out = a.out or (snakemake.output.json if snakemake else None)
    req_path = a.requirements
    if not sbml or not out:
        sys.exit("FATAL: --sbml and --out are required outside Snakemake")

    try:
        import cobra
    except ImportError:
        sys.exit("cobra not importable; see envs/cobra.yaml")

    model = cobra.io.read_sbml_model(sbml)
    base_growth = float(model.optimize().objective_value or 0.0)

    bofs = pick_bof(model)
    print("=" * 76)
    print("iPC815 media definition")
    print("=" * 76)
    print(f"model           : {model.id}")
    print(f"reactions       : {len(model.reactions)}  "
          f"exchanges: {len(model.exchanges)}")
    print(f"default growth  : {base_growth:.6f}")
    print(f"biomass objectives found: {bofs}")

    bof_37 = bofs.get(37)
    bof_25 = bofs.get(25) or bofs.get(26)
    if not bof_37 or not bof_25:
        print()
        print("WARNING: expected both a 37 C and a 25/26 C biomass objective.")
        print("         Found: " + str(bofs))
        print("         Falling back to a single objective for every medium.")
        print("         Record in Limitations that biomass composition was")
        print("         NOT adjusted for temperature -- that weakens the")
        print("         mammal-versus-flea contrast and must be stated.")
        bof_37 = bof_37 or bofs.get("untagged")
        bof_25 = bof_25 or bof_37

    iron = exchanges_matching(model, pattern=IRON)
    glucose = exchanges_matching(model, pattern=GLUCOSE)
    aas = exchanges_matching(model, names=AMINO_ACIDS)
    carbon = exchanges_matching(model, names=CARBON_SCAN_HINTS)

    print()
    print(f"iron exchanges found    : {len(iron)}  "
          f"{[i for i, _ in iron][:6]}")
    print(f"glucose exchanges found : {len(glucose)} "
          f"{[i for i, _ in glucose][:4]}")
    print(f"amino-acid exchanges    : {len(aas)}")
    print(f"carbon-source candidates: {len(carbon)}")

    missing = []
    if not iron:
        missing.append("iron (blood medium cannot be made iron-limited, "
                       "so ybt/HPI necessity will not be reproduced)")
    if not glucose:
        missing.append("glucose (blood medium cannot be glucose-limited)")
    if len(aas) < 10:
        missing.append(f"amino acids (only {len(aas)} found; the blood "
                       f"medium will not be amino-acid rich)")

    default_medium = dict(model.medium)

    media = {}

    # ---- mammalian blood: AA-rich, glucose-limited, iron withheld --------
    blood = dict(default_medium)
    for rid, _ in aas:
        blood[rid] = 10.0
    for rid, _ in glucose:
        blood[rid] = a.glucose_limit
    for rid, _ in iron:
        blood[rid] = 0.0
    media["mammalian_blood"] = {
        "objective": bof_37,
        "temperature_c": 37,
        "medium": blood,
        "rationale": ("amino-acid rich, glucose limited, iron strictly "
                      "withheld to represent host iron sequestration; "
                      "explains the necessity of ybt/HPI"),
    }

    # ---- flea midgut: blood meal at 26 C --------------------------------
    flea = dict(default_medium)
    for rid, _ in aas:
        flea[rid] = 5.0
    for rid, _ in glucose:
        flea[rid] = 0.5
    media["flea_midgut"] = {
        "objective": bof_25,
        "temperature_c": 26,
        "medium": flea,
        "rationale": ("blood-meal components at vector temperature; the "
                      "25 C biomass objective is what makes the "
                      "mammal-versus-flea contrast meaningful"),
    }

    # ---- oligotrophic: minimal, single carbon source scan ---------------
    # The base keeps the AUXOTROPHIC REQUIREMENTS open. Measured 2026-10-10:
    # iPC815 cannot synthesise phenylalanine (and needs methionine and O2), so
    # zeroing every carbon-source candidate also zeroed phe and made growth
    # impossible for ANY of the 76 candidates. A medium that kills the wild
    # type makes every knockout read lethal and Fig.5 uniformly zero, with no
    # error raised -- so the requirements are part of the base, not part of
    # the variable.
    requirements = {}
    if req_path and Path(req_path).exists():
        with open(req_path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("role") == "AUXOTROPHIC_REQUIREMENT":
                    requirements[r["exchange"]] = r.get("metabolite", "")
    else:
        print(f"WARNING: {req_path} not found -- run "
              f"scripts/essential_vs_substrate_scan.py first. Falling back "
              f"to an empty requirement set, which reproduces the "
              f"zero-growth defect.", file=sys.stderr)

    oligo = dict(default_medium)
    for rid, _ in carbon:
        if rid in requirements:
            continue
        oligo[rid] = 0.0
    media["oligotrophic"] = {
        "objective": bof_25,
        "temperature_c": 26,
        "medium": oligo,
        "carbon_source_scan": [rid for rid, _ in carbon],
        "auxotrophic_requirements": requirements,
        "rationale": ("minimal medium with every candidate carbon source "
                      "closed but the auxotrophic requirements open; "
                      "strain_fba.py reopens one carbon source at a time, "
                      "and the readout is the COUNT of usable sources, not "
                      "a single growth rate"),
    }

    # ---- verify each medium supports wild-type growth -------------------
    print()
    print("wild-type growth check")
    failures = []
    for name, spec in media.items():
        with model:
            model.objective = spec["objective"]
            try:
                model.medium = {k: v for k, v in spec["medium"].items()
                                if v > 0}
            except Exception as exc:                      # noqa: BLE001
                print(f"  {name:<18} medium rejected: {exc}")
                failures.append(name)
                continue
            g = float(model.optimize().objective_value or 0.0)
        spec["wt_growth"] = g
        flag = "OK" if g > 1e-9 else "ZERO"
        print(f"  {name:<18}{spec['objective']:<16}{g:>12.6f}  [{flag}]")
        if g <= 1e-9 and name != "oligotrophic":
            failures.append(name)

    # oligotrophic carries no carbon source BY DESIGN, so a zero on the base is
    # expected and says nothing. What must not be zero is the number of carbon
    # sources that can restore growth -- that is the guard against a medium
    # that kills the wild type (every knockout would read lethal and Fig.5
    # would be a uniform zero heatmap with no error raised).
    if "oligotrophic" in media and carbon:
        tried = []
        for rid, mname in carbon:
            with model:
                model.objective = media["oligotrophic"]["objective"]
                med = {k: v for k, v in media["oligotrophic"]["medium"].items()
                       if v > 0}
                med[rid] = 10.0
                try:
                    model.medium = med
                    g = float(model.optimize().objective_value or 0.0)
                except Exception:                          # noqa: BLE001
                    g = 0.0
            tried.append({"exchange": rid, "metabolite": mname, "growth": g})
        usable = [t for t in tried if t["growth"] > 1e-9]
        for t in usable:
            print(f"  {'oligotrophic+C':<18}{t['exchange']:<16}"
                  f"{t['growth']:>12.6f}  [OK]  {t['metabolite']}")
        print(f"  -> {len(usable)}/{len(tried)} single carbon sources "
              f"support growth")
        print(f"     requirements held open: "
              f"{', '.join(sorted(requirements)) or '(none)'}")
        media["oligotrophic"]["single_carbon_trials"] = tried
        media["oligotrophic"]["n_usable_carbon_sources"] = len(usable)
        media["oligotrophic"]["usable_carbon_sources"] = [
            {"exchange": t["exchange"], "metabolite": t["metabolite"],
             "growth": t["growth"]} for t in usable]
        if not usable:
            failures.append(
                f"oligotrophic (none of the {len(tried)} carbon-source "
                "candidates supports growth with the auxotrophic "
                "requirements open -- the medium, not the strain, is wrong)")

    Path(out).parent.mkdir(parents=True, exist_ok=True)
    Path(out).write_text(json.dumps(
        {"model": model.id, "source_sbml": sbml,
         "default_growth": base_growth, "biomass_objectives": bofs,
         "exchanges_not_found": missing, "media": media}, indent=2))

    print()
    if missing:
        print("EXCHANGES NOT FOUND -- each one weakens a specific claim:")
        for m in missing:
            print(f"  - {m}")
        print("  Inspect the model's exchange list and extend the name")
        print("  patterns at the top of this script rather than letting the")
        print("  medium silently lack the constraint.")
    if failures:
        print()
        print("FATAL: these media do not support wild-type growth: "
              f"{failures}")
        print("Every strain knockout would read as lethal and Fig.5 would be")
        print("a uniform zero heatmap with no error raised. Fix the medium")
        print("before running strain_fba.py.")
        print(f"written anyway for inspection: {out}")
        sys.exit(1)

    print(f"written: {out}")
    print("=" * 76)


if __name__ == "__main__":
    main()
