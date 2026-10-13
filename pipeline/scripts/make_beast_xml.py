#!/usr/bin/env python3
"""
Build the BEAST2 XML, and refuse to build it wrong.

THE SILENT KILLER THIS SCRIPT EXISTS TO PREVENT
-----------------------------------------------
BEAST2 treats the supplied alignment as the complete set of sites. Hand it
~25,000 variable sites with no constant-site weights and every branch length
is inflated by roughly

    4,650,000 constant+variable sites / 25,000 variable sites  ~=  186x

The run still converges. Tracer still shows healthy ESS. Every clock rate and
every node age is wrong by more than two orders of magnitude, and nothing in
the output says so. The only defence is to supply constantSiteWeights (the
AscertainedAlignment machinery) and to assert they are present before the XML
is ever written.

The same four counts come from `snp-sites -C`, which also feeds iqtree
-fconst. Note that -fconst and +ASC are NOT interchangeable: +ASC conditions
the likelihood on constant sites being unobservable, which is the wrong model
when they have been observed and counted.
"""

import csv
import sys
from pathlib import Path
from xml.sax.saxutils import escape

snakemake = globals().get("snakemake")  # noqa: F821

MIN_PLAUSIBLE_CONSTANT_SITES = 1_000_000   # Y. pestis chromosome is ~4.65 Mb


def read_fconst(path: Path) -> list[int]:
    raw = path.read_text().strip()
    parts = [p for p in raw.replace("\n", ",").split(",") if p.strip()]
    counts = [int(float(p)) for p in parts]
    if len(counts) != 4:
        raise SystemExit(
            f"FATAL: expected 4 constant-site counts (A,C,G,T), got {counts} "
            f"from {path}. Regenerate with `snp-sites -C`."
        )
    return counts


def read_fasta(path: Path) -> dict[str, str]:
    seqs, name, buf = {}, None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def main():
    fconst = read_fconst(Path(snakemake.input.fconst))
    total_const = sum(fconst)

    # ---- assertion 1: the weights exist and are biologically plausible ----
    if total_const < MIN_PLAUSIBLE_CONSTANT_SITES:
        raise SystemExit(
            f"FATAL: constant-site total {total_const:,} is implausibly low "
            f"for a ~4.65 Mb genome. Branch lengths would be inflated and the "
            f"clock rate meaningless. Check that snippy-core produced a "
            f"FULL-LENGTH alignment before snp-sites -C."
        )

    seqs = read_fasta(Path(snakemake.input.var))
    keep = {r["accession"] for r in
            csv.DictReader(open(snakemake.input.taxa), delimiter="\t")}
    seqs = {k: v for k, v in seqs.items() if k in keep}
    if not seqs:
        raise SystemExit("FATAL: no taxa left after filtering to dating set.")

    n_var = len(next(iter(seqs.values())))
    inflation = (total_const + n_var) / max(n_var, 1)

    dates = {r["taxon"]: r["date"] for r in
             csv.DictReader(open(snakemake.input.dates), delimiter="\t")
             if r["taxon"] in seqs}
    missing = sorted(set(seqs) - set(dates))
    if missing:
        raise SystemExit(
            f"FATAL: {len(missing)} taxa have no tip date, e.g. {missing[:5]}. "
            f"Drop them from the dating set or supply dates."
        )

    # ------------------------------------------------------------------ XML
    seq_xml = "\n".join(
        f'        <sequence id="seq_{escape(t)}" taxon="{escape(t)}" '
        f'totalcount="4" value="{s}"/>' for t, s in sorted(seqs.items()))

    trait = ",".join(f"{t}={d}" for t, d in sorted(dates.items()))

    xml = f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<beast beautitemplate='Standard' beautistatus='' namespace="beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood" required="" version="2.7">

    <!--
      Constant-site weights: A,C,G,T = {fconst[0]},{fconst[1]},{fconst[2]},{fconst[3]}
      total constant = {total_const:,}
      variable sites = {n_var:,}
      Without these weights every branch length would be inflated ~{inflation:.0f}x.
    -->

    <data id="alignment" spec="FilteredAlignment" filter="-">
      <data id="rawAlignment" spec="Alignment" dataType="nucleotide">
{seq_xml}
      </data>
      <constantSiteWeights id="constantSiteWeights" spec="parameter.IntegerParameter"
                           dimension="4" estimate="false">{fconst[0]} {fconst[1]} {fconst[2]} {fconst[3]}</constantSiteWeights>
    </data>

    <trait id="dateTrait" spec="beast.evolution.tree.TraitSet"
           traitname="date" value="{trait}">
      <taxa id="TaxonSet" spec="TaxonSet" alignment="@alignment"/>
    </trait>

    <!-- clock: {snakemake.params.clock}  site: {snakemake.params.site}  prior: {snakemake.params.prior} -->
    <!-- NOTE: install the BEAST2 packages first:
         packagemanager -add ORC -add RLC -add bModelTest -->

    <run id="mcmc" spec="MCMC" chainLength="{snakemake.params.length}"
         storeEvery="{snakemake.params.store}">
      <!-- state, prior, operators and loggers: emit from a BEAUti template
           keyed on clock={snakemake.params.clock}, site={snakemake.params.site},
           treePrior={snakemake.params.prior}.
           The data and trait blocks above are already complete and correct. -->
    </run>
</beast>
"""

    Path(snakemake.output.xml).write_text(xml)

    msg = "\n".join([
        "=" * 72,
        "BEAST2 XML written",
        "=" * 72,
        f"taxa                 : {len(seqs)}",
        f"variable sites       : {n_var:,}",
        f"constant sites ACGT  : {fconst} (total {total_const:,})",
        f"branch-length inflation avoided: ~{inflation:.0f}x",
        "constantSiteWeights  : PRESENT (asserted)",
        "=" * 72,
    ])
    print(msg, file=sys.stderr)
    Path(snakemake.log[0]).write_text(msg + "\n")


if __name__ == "__main__":
    main()
