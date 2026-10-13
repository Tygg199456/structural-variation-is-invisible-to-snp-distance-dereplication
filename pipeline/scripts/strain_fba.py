#!/usr/bin/env python3
"""
Per-strain and per-ancestral-node knockout + pFBA on iPC815.

DECISION THAT MUST NOT BE LOST:
Let GPR rules propagate; do not knock out reactions directly. Record which
lost genes have NO GPR -- that is the 20% coverage caveat made
quantitative.

Inputs / outputs are supplied by Snakemake via the `snakemake` object.
"""

snakemake = globals().get("snakemake")  # noqa: F821


def main():
    raise NotImplementedError(
        "strain_fba.py: not implemented. See the module docstring for the contract."
    )


if __name__ == "__main__":
    main()
