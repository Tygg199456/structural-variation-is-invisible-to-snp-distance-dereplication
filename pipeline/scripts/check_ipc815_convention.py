#!/usr/bin/env python3
"""
Decide whether the BiGG iPC815 download differs from the publication in
SUBSTANCE or only in COUNTING CONVENTION.

Measured on the BiGG file:
    genes       815   vs published  815   -- exact match
    species    1552   vs published  963   -- +589
    reactions  1961   vs published 1678   -- +283

The exact gene match is the strong evidence. A different reconstruction
would not land on 815 genes; a curation expansion would normally ADD genes
while filling gaps, not hold them constant while adding several hundred
reactions. That points at counting convention, and there are two standard
ones to test:

  METABOLITES -- BiGG counts species PER COMPARTMENT. Glucose in cytosol,
  periplasm and extracellular space are three species records
  (M_glc__D_c / _p / _e). The paper reports 963 UNIQUE metabolites, i.e.
  distinct chemical entities with compartments collapsed. With 3
  compartments, 1552 records collapsing to ~963 is exactly what this
  convention produces.

  REACTIONS -- BiGG models always carry exchange reactions (R_EX_*), usually
  demand and sink reactions too. A published internal-reaction count would
  exclude them. The gap is 283, which is the right order for the exchange
  set of a three-compartment model.

VERDICT RULES
  both hit   -> convention difference. The model IS the published iPC815.
                Fig.5 needs no downgrade on this account.
  one hits   -> partial explanation; record the residual in Methods.
  neither    -> genuinely different model. Redesign Fig.5 around a de novo
                reconstruction (gapseq / CarveMe) before investing further.

The ~20% gene coverage caveat is unaffected either way -- that was always
the real ceiling on Fig.5, and it is unchanged by this result.

Usage:  python3 scripts/check_ipc815_convention.py data/models/iPC815.xml
"""
import collections
import re
import sys

PUB = {"genes": 815, "metabolites": 963, "reactions": 1678}
TOL = 0.03          # 3% tolerance on the "hit" judgement


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "data/models/iPC815.xml"
    text = open(path).read()

    species = re.findall(r'<species\b[^>]*\bid="([^"]+)"', text)
    reactions = re.findall(r'<reaction\b[^>]*\bid="([^"]+)"', text)
    genes = re.findall(r'<fbc:geneProduct\b[^>]*\bfbc:id="([^"]+)"', text)
    comps = re.findall(r'<compartment\b[^>]*\bid="([^"]+)"', text)

    # BiGG: M_<base>_<compartment letter>
    def strip_compartment(s):
        s = s[2:] if s.startswith("M_") else s
        return re.sub(r"_[a-z0-9]{1,2}$", "", s)

    unique_mets = {strip_compartment(s) for s in species}
    boundary = [r for r in reactions
                if r.startswith(("R_EX_", "R_DM_", "R_SK_", "EX_", "DM_", "SK_"))]
    internal = len(reactions) - len(boundary)

    suffixes = collections.Counter(s.rsplit("_", 1)[-1] for s in species)

    def hit(obs, exp):
        return abs(obs - exp) <= TOL * exp

    met_hit = hit(len(unique_mets), PUB["metabolites"])
    rxn_hit = hit(internal, PUB["reactions"])
    gene_hit = len(genes) == PUB["genes"]

    print("=" * 70)
    print(f"iPC815 counting-convention check  --  {path}")
    print("=" * 70)
    print(f"compartments declared     : {len(comps)}  {comps}")
    print(f"species compartment tags  : {dict(suffixes.most_common(8))}")
    print()
    print(f"{'genes':<28}{len(genes):>6}   published {PUB['genes']:>5}"
          f"   [{'OK' if gene_hit else 'MISMATCH'}]")
    print(f"{'species records':<28}{len(species):>6}   "
          f"(compartment-resolved)")
    print(f"{'  -> unique metabolites':<28}{len(unique_mets):>6}   "
          f"published {PUB['metabolites']:>5}   "
          f"[{'HIT' if met_hit else 'MISS'}]")
    print(f"{'reactions total':<28}{len(reactions):>6}")
    print(f"{'  -> exchange/demand/sink':<28}{len(boundary):>6}")
    print(f"{'  -> internal only':<28}{internal:>6}   "
          f"published {PUB['reactions']:>5}   "
          f"[{'HIT' if rxn_hit else 'MISS'}]")
    print()

    if gene_hit and met_hit and rxn_hit:
        verdict = (
            "CONVENTION DIFFERENCE. This is the published iPC815.\n"
            "  Methods: state that counts are reported compartment-resolved\n"
            "  and including boundary reactions, and give both figures.\n"
            "  Fig.5 evidence strength: UNCHANGED.")
        code = 0
    elif gene_hit and (met_hit or rxn_hit):
        verdict = (
            "PARTIAL. Gene set matches and one count is explained by\n"
            "  convention; the other is not. Report the residual explicitly\n"
            "  in Methods and check the BiGG model's version history.\n"
            "  Fig.5 evidence strength: slightly reduced, still usable.")
        code = 0
    elif gene_hit:
        verdict = (
            "UNEXPLAINED. Gene set matches but neither count is explained by\n"
            "  convention, so this is likely a curated expansion. Usable, but\n"
            "  cite the BiGG version, not the 2011 counts.")
        code = 0
    else:
        verdict = (
            "DIFFERENT MODEL. The gene set does not match. Do NOT present\n"
            "  this as iPC815. Redesign Fig.5 around a de novo reconstruction\n"
            "  (gapseq / CarveMe) and decide that NOW, not after the\n"
            "  phylogenetics is finished.")
        code = 1

    print("VERDICT")
    print("  " + verdict)
    print()
    print("UNCHANGED REGARDLESS: 815 genes is ~20% of the ~4,000-gene genome.")
    print("Lost genes with no GPR are invisible to FBA. Scope every Fig.5")
    print("conclusion to metabolic MODULES and have strain_fba.py count the")
    print("lost genes it could not represent.")
    print("=" * 70)
    sys.exit(code)


if __name__ == "__main__":
    main()
