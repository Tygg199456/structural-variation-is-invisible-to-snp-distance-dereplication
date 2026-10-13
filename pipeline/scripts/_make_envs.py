#!/usr/bin/env python3
"""Write the conda env files for the pipeline. Run once: python scripts/_make_envs.py"""
import pathlib

CHANNELS = ["conda-forge", "bioconda", "defaults"]

ENVS = {
    "ncbi":        ["ncbi-datasets-cli=18.38.0", "unzip", "jq"],
    "py":          ["python=3.11", "numpy=1.26", "pandas=2.2", "scipy=1.13",
                    "matplotlib=3.8", "seaborn=0.13", "biopython=1.83",
                    "dendropy=4.6", "ete3=3.1", "pysam=0.22", "scikit-learn=1.5",
                    "statsmodels=0.14", "pycirclize", "requests"],
    "r":           ["r-base=4.3", "bioconductor-ggtree", "r-ggplot2", "r-dplyr",
                    "bioconductor-clusterprofiler", "r-ape", "r-phytools"],
    "mummer":      ["mummer4=4.0"],
    "checkm2":     ["checkm2=1.0.2"],
    "snippy":      ["snippy=4.6.0"],
    "bedtools":    ["bedtools=2.31"],
    "gubbins":     ["gubbins=3.3"],
    "snpsites":    ["snp-sites=2.5.1", "snp-dists=0.8.2"],
    "iqtree":      ["iqtree=2.3.6"],
    "veryfasttree":["veryfasttree=4.0"],
    "bakta":       ["bakta=1.9.4"],
    "panaroo":     ["panaroo=1.5.0"],
    "ppanggolin":  ["ppanggolin=2.1.0"],
    "pseudofinder":["pseudofinder=1.1", "diamond=2.1"],
    "isescan":     ["isescan=1.7.2.3"],
    "pastml":      ["pastml=1.9.34"],
    "beast":       ["beast2=2.7.7", "beagle-lib", "tracer"],
    "cobra":       ["python=3.11", "cobra=0.29", "glpk", "python-libsbml",
                    "pandas", "numpy"],
    "gapseq":      ["gapseq", "carveme"],
    "amdirt":      ["amdirt"],
    "nextflow":    ["nextflow=23.10", "openjdk=17"],
}

NOTES = {
    "beast": "# BEAST2 packages ORC, RLC and bModelTest must be installed\n"
             "# separately:  packagemanager -add ORC -add RLC -add bModelTest\n",
    "bakta": "# Set $BAKTA_DB and record the exact DB release in config.yaml\n"
             "# (annotation.db_version). One DB version for every genome.\n",
    "pseudofinder": "# Build the DIAMOND reference db from Y. pseudotuberculosis\n"
                    "# plus related Enterobacterales before first use.\n",
}

out = pathlib.Path(__file__).resolve().parent.parent / "envs"
out.mkdir(exist_ok=True)

for name, deps in ENVS.items():
    lines = [NOTES.get(name, "").rstrip("\n")] if name in NOTES else []
    lines += [f"name: {name}", "channels:"]
    lines += [f"  - {c}" for c in CHANNELS]
    lines += ["dependencies:"]
    lines += [f"  - {d}" for d in deps]
    (out / f"{name}.yaml").write_text("\n".join(x for x in lines if x) + "\n")

print(f"wrote {len(ENVS)} env files to {out}")
