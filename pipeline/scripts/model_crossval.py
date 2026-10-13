#!/usr/bin/env python3
"""
gapseq / CarveMe reconstructions for 10 representatives.

DECISION THAT MUST NOT BE LOST:
If de novo models contradict the iPC815 trend, demote the FBA panel to
supplementary rather than arguing around it.

Inputs / outputs are supplied by Snakemake via the `snakemake` object.
"""

snakemake = globals().get("snakemake")  # noqa: F821


def main():
    raise NotImplementedError(
        "model_crossval.py: not implemented. See the module docstring for the contract."
    )


if __name__ == "__main__":
    main()
