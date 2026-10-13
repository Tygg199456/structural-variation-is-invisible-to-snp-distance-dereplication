#!/usr/bin/env python3
"""
Does a gene get lost more often than drift alone would produce?

The reductive-evolution literature asserts convergent loss without ever
testing it against a null. This is that test, and it is what makes the
parallelism claim in this paper falsifiable.

Procedure, per gene:
  1. Count independent loss events on the reconstructed tree (from PastML).
  2. Estimate that gene's per-unit-branch-length loss rate by ML.
  3. Simulate loss along the SAME tree, at that rate, N times.
  4. Empirical p = P(simulated independent losses >= observed).
  5. Benjamini-Hochberg across genes.

Step 2 matters: simulating every gene at one shared rate would flag any
frequently-lost gene as parallel. Conditioning on the gene's own rate asks the
sharper question -- given how often this gene is lost at all, is it lost in
more INDEPENDENT places than chance allows?
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

snakemake = globals().get("snakemake")  # noqa: F821

try:
    import dendropy
except ImportError:  # pragma: no cover
    dendropy = None


def load_tree(path):
    if dendropy is None:
        raise SystemExit("dendropy required; see envs/py.yaml")
    return dendropy.Tree.get(path=path, schema="newick",
                             preserve_underscores=True)


def simulate_losses(branch_lengths, rate, rng, n_sim):
    """
    Independent-loss count per simulation under a Poisson loss process.
    A branch either carries a loss or it does not; the number of distinct
    branches carrying a loss is the independent-loss count.
    """
    p = 1.0 - np.exp(-rate * branch_lengths)          # (n_branches,)
    draws = rng.random((n_sim, p.size)) < p
    return draws.sum(axis=1)


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def main():
    n_sim = int(snakemake.params.sims)
    alpha = float(snakemake.params.alpha)
    rng = np.random.default_rng(20261009)

    tree = load_tree(snakemake.input.tree)
    bl = np.array([e.length for e in tree.edges()
                   if e.length is not None and e.length > 0], dtype=float)
    total_bl = bl.sum()
    if total_bl <= 0:
        raise SystemExit("FATAL: tree has no positive branch lengths.")

    observed = defaultdict(int)
    n_ambiguous_skipped = 0
    n_rows = 0
    for r in csv.DictReader(open(snakemake.input.events), delimiter="\t"):
        n_rows += 1
        # Ambiguous transitions carry event='ambiguous' and must NOT enter the
        # loss count: this null estimates each gene's rate FROM the observed
        # count, so counting unresolved transitions inflates numerator and
        # rate together and the test loses power.
        conf = (r.get("confidence") or "unambiguous").strip().lower()
        if conf != "unambiguous":
            n_ambiguous_skipped += 1
            continue
        if r.get("event", "").lower() in ("loss", "pseudogenisation",
                                          "pseudogenization"):
            observed[r["gene"]] += 1

    print(f"[parallel_loss_null] event rows={n_rows} "
          f"ambiguous_skipped={n_ambiguous_skipped} "
          f"({100 * n_ambiguous_skipped / max(n_rows, 1):.1f}%) "
          f"genes_with_losses={len(observed)}", file=sys.stderr)

    if not observed:
        raise SystemExit("FATAL: no loss events parsed from PastML output.")

    rows = []
    for gene, obs in sorted(observed.items()):
        # ML rate under a Poisson loss process over total tree length
        rate = obs / total_bl
        sims = simulate_losses(bl, rate, rng, n_sim)
        p = float((sims >= obs).sum() + 1) / (n_sim + 1)   # add-one, never 0
        rows.append({
            "gene": gene,
            "observed_independent_losses": obs,
            "ml_loss_rate_per_bl": f"{rate:.6g}",
            "sim_mean": f"{sims.mean():.3f}",
            "sim_p95": f"{np.percentile(sims, 95):.1f}",
            "p_empirical": f"{p:.6g}",
        })

    q = bh_fdr([float(r["p_empirical"]) for r in rows])
    for r, qq in zip(rows, q):
        r["q_BH"] = f"{qq:.6g}"
        r["significant"] = "yes" if qq < alpha else "no"

    rows.sort(key=lambda r: (float(r["q_BH"]),
                             -int(r["observed_independent_losses"])))

    with open(snakemake.output.tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    n_sig = sum(1 for r in rows if r["significant"] == "yes")
    msg = (f"[parallel_loss_null] genes={len(rows)}  sims={n_sim}  "
           f"total_branch_length={total_bl:.4f}  "
           f"significant_at_FDR{alpha}={n_sig}")
    print(msg, file=sys.stderr)
    Path(snakemake.log[0]).write_text(msg + "\n")


if __name__ == "__main__":
    main()
