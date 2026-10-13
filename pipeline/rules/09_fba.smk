# ------------------------------------------------------------------ 09_fba ---
# Metabolic consequence of gene decay.
#
# iPC815 (Charusanti 2011): 815 genes, 678 proteins, 963 metabolites,
# 1678 reactions -- roughly 20% of the ~4000-gene genome. Many lost genes
# therefore have no GPR in the model and their effect is invisible to FBA.
# Scope every conclusion to metabolic MODULES, state the coverage in the
# Methods, and cross-check with de novo reconstructions.
#
# Cheapest step in the whole pipeline: < 30 min for 210 strains x 3 media.

rule fetch_ipc815:
    output:
        sbml = f"{DATA}/models/iPC815.xml"
    log: f"{LOGS}/09_fetch_model.log"
    conda: "../envs/py.yaml"
    script: "../scripts/fetch_ipc815.py"


rule define_media:
    """
    Three in silico media, each with an explicit biological rationale that
    belongs verbatim in the Methods:
      mammalian blood -- amino-acid rich, glucose limited, iron withheld
                         (host iron sequestration; explains ybt/HPI necessity)
      flea midgut     -- blood-meal components at 26 C, high osmolarity
                         (biofilm/hms-associated metabolism)
      oligotrophic    -- minimal medium, 30 single carbon sources scanned
                         (when environmental persistence was lost)
    """
    input:
        sbml = f"{DATA}/models/iPC815.xml"
    output:
        json = f"{RES}/fba/media.json"
    params:
        media = config["fba"]["media"]
    log: f"{LOGS}/09_define_media.log"
    conda: "../envs/cobra.yaml"
    script: "../scripts/define_media.py"


rule strain_knockout_fba:
    """
    For each Tier A strain and each reconstructed ancestral node: apply the
    lost/pseudogenised gene set, let the GPR rules propagate, then pFBA.
    """
    input:
        sbml   = f"{DATA}/models/iPC815.xml",
        media  = f"{RES}/fba/media.json",
        states = f"{RES}/asr/combined_states.tsv",
        events = f"{RES}/asr/pastml_events.tsv"
    output:
        growth = f"{RES}/fba/growth_matrix.tsv",
        auxo   = f"{RES}/fba/auxotrophy_matrix.tsv"
    threads: 16
    params:
        solver = config["fba"]["solver"]
    log: f"{LOGS}/09_knockout_fba.log"
    conda: "../envs/cobra.yaml"
    script: "../scripts/strain_fba.py"


rule model_cross_validation:
    """
    De novo reconstructions for 10 representatives. If gapseq/CarveMe models
    reproduce the iPC815 trend, the 20% coverage caveat is survivable; if they
    diverge, the FBA panel must be demoted to supplementary.
    """
    input:
        faas = lambda w: expand(f"{RES}/annotate/{{acc}}/{{acc}}.faa",
                                acc=fba_representatives()),
        ipc  = f"{RES}/fba/growth_matrix.tsv"
    output:
        tsv = f"{RES}/fba/cross_validation.tsv"
    threads: 16
    log: f"{LOGS}/09_cross_validation.log"
    conda: "../envs/gapseq.yaml"
    script: "../scripts/model_crossval.py"


def fba_representatives():
    """One representative per major lineage, anchors preferred."""
    p = f"{RES}/fba/representatives.txt"
    if not os.path.exists(p):
        return []
    return [ln.strip() for ln in open(p) if ln.strip()]
