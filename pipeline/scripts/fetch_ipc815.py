#!/usr/bin/env python3
"""
Fetch and validate iPC815, the Y. pestis CO92 genome-scale metabolic model.

RUN THIS OUT OF DAG ORDER, EARLY.
FBA sits at the end of the pipeline and costs under 30 minutes, but carries
the second of the study's two novelty claims. If iPC815 is unusable, Fig.5
needs redesigning -- a decision worth three weeks early rather than three
weeks late.

v3: curl/wget first. On this machine Python's TLS stack fails
(`SSL: UNEXPECTED_EOF_WHILE_READING`, urlopen timeouts) while curl on the
same host works, so urllib is the fallback, not the default. Any future
remote fetch in this pipeline should follow the same order.

Published composition (Charusanti et al. 2011, BMC Syst Biol 5:163):
    815 genes | 678 proteins | 963 unique metabolites | 1678 reactions

Counts from a repository download may legitimately differ by COUNTING
CONVENTION rather than substance -- BiGG counts metabolites per compartment
and includes exchange reactions. `check_ipc815_convention.py` tests that
directly; this script only reports the raw comparison and the baseline
growth check.
"""
import json
import shutil
import subprocess
import sys
from pathlib import Path

snakemake = globals().get("snakemake")  # noqa: F821

EXPECTED = {"genes": 815, "reactions": 1678, "metabolites": 963}

SOURCES = [
    ("BiGG", "http://bigg.ucsd.edu/static/models/iPC815.xml"),
    ("BioModels", "https://www.ebi.ac.uk/biomodels/model/download/"
                  "MODEL1507180045?filename=MODEL1507180045_url.xml"),
]

MANUAL_NOTE = """
Every candidate source failed.

Retrieve iPC815 manually, place the SBML at the output path, and re-run --
the script short-circuits on an existing file and goes straight to
validation.

  1. BiGG Models   http://bigg.ucsd.edu/models        search 'iPC815'
  2. BioModels     https://www.ebi.ac.uk/biomodels/   search 'Yersinia pestis'
  3. Supplementary material of Charusanti et al. 2011,
     BMC Systems Biology 5:163 (open access)

If none yields a loadable SBML, Fig.5 must be redesigned around a de novo
reconstruction (gapseq / CarveMe). Decide that now.
"""

MIN_BYTES = 10_000


def fetch(url: str, dest: Path) -> str | None:
    """curl, then wget, then urllib. Returns the tool that worked."""
    tmp = dest.with_suffix(dest.suffix + ".part")

    if shutil.which("curl"):
        r = subprocess.run(
            ["curl", "-fsSL", "--retry", "3", "--max-time", "300",
             "-o", str(tmp), url],
            capture_output=True, text=True)
        if r.returncode == 0 and tmp.exists() and tmp.stat().st_size > MIN_BYTES:
            tmp.replace(dest)
            return "curl"
        print(f"    curl failed: rc={r.returncode} {r.stderr.strip()[:200]}",
              file=sys.stderr)

    if shutil.which("wget"):
        r = subprocess.run(
            ["wget", "-q", "-T", "300", "-t", "3", "-O", str(tmp), url],
            capture_output=True, text=True)
        if r.returncode == 0 and tmp.exists() and tmp.stat().st_size > MIN_BYTES:
            tmp.replace(dest)
            return "wget"
        print(f"    wget failed: rc={r.returncode}", file=sys.stderr)

    try:
        import urllib.request
        with urllib.request.urlopen(url, timeout=120) as resp:
            data = resp.read()
        if len(data) > MIN_BYTES:
            dest.write_bytes(data)
            return "urllib"
    except Exception as exc:                      # noqa: BLE001
        print(f"    urllib failed: {type(exc).__name__}: {exc}",
              file=sys.stderr)

    if tmp.exists():
        tmp.unlink()
    return None


def main():
    out = Path(snakemake.output.sbml)
    out.parent.mkdir(parents=True, exist_ok=True)

    via = None
    if out.exists() and out.stat().st_size > MIN_BYTES:
        print(f"[fetch_ipc815] using existing {out} "
              f"({out.stat().st_size:,} bytes)", file=sys.stderr)
        via = "pre-existing"
    else:
        for name, url in SOURCES:
            print(f"[fetch_ipc815] trying {name}: {url}", file=sys.stderr)
            tool = fetch(url, out)
            if tool:
                via = f"{name} (via {tool})"
                break

    if via is None:
        raise SystemExit(MANUAL_NOTE)

    try:
        import cobra
    except ImportError:
        raise SystemExit(
            "cobra not importable. Install it before this step:\n"
            "  mamba create -p ~/envs/cobra -c conda-forge cobra python-libsbml\n"
            "The baseline growth check below is not optional: if the "
            "unmodified model does not grow, every strain knockout reads as "
            "lethal and Fig.5 is a uniform zero heatmap with no error raised."
        )

    model = cobra.io.read_sbml_model(str(out))
    observed = {"genes": len(model.genes),
                "reactions": len(model.reactions),
                "metabolites": len(model.metabolites)}

    growth = float(model.optimize().objective_value or 0.0)

    lines = ["=" * 72, "iPC815 validation", "=" * 72,
             f"source                 : {via}",
             f"model id               : {model.id}"]
    for k in ("genes", "reactions", "metabolites"):
        flag = "OK" if observed[k] == EXPECTED[k] else "DIFFERS"
        lines.append(f"{k:<23}: {observed[k]} "
                     f"(published {EXPECTED[k]}) [{flag}]")
    lines.append(f"exchange reactions     : {len(model.exchanges)}")
    lines.append(f"baseline growth rate   : {growth:.6f}")

    if observed != EXPECTED:
        lines += ["",
                  "Counts differ from the publication. Before treating this",
                  "as a different model, run:",
                  f"    python3 scripts/check_ipc815_convention.py {out}",
                  "BiGG counts metabolites per compartment and includes",
                  "exchange reactions; both are counting conventions, not",
                  "substance. An exact gene match is the signal that the",
                  "underlying reconstruction is unchanged."]

    if growth <= 1e-9:
        lines += ["",
                  "FATAL: the unmodified model does not grow on its default",
                  "medium. Strain knockouts would all read as lethal and",
                  "Fig.5 would be meaningless. Fix the medium before",
                  "proceeding to define_media.py."]

    lines += ["",
              f"gene coverage vs ~4000-gene genome: "
              f"{100 * observed['genes'] / 4000:.1f}%",
              "Scope every FBA conclusion to metabolic MODULES. Lost genes",
              "with no GPR are invisible here; strain_fba.py must count them",
              "explicitly rather than letting them read as 'no effect'.",
              "=" * 72]

    msg = "\n".join(lines)
    print(msg, file=sys.stderr)
    Path(snakemake.log[0]).write_text(msg + "\n")
    Path(str(out) + ".validation.json").write_text(json.dumps(
        {"source": via, "model_id": model.id, "expected": EXPECTED,
         "observed": observed, "matches_publication": observed == EXPECTED,
         "baseline_growth_rate": growth,
         "n_exchange_reactions": len(model.exchanges)}, indent=2))

    if growth <= 1e-9:
        sys.exit(1)


if __name__ == "__main__":
    main()
