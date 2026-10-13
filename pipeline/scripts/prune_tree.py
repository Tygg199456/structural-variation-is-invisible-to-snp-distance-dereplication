#!/usr/bin/env python3
"""
Prune the dated MCC tree to Tier A + Tier C for ASR.

DECISION THAT MUST NOT BE LOST:
Keep dated branch lengths: loss events must land on a real time axis even
though only Tier A carries trustworthy PAV.

Inputs / outputs are supplied by Snakemake via the `snakemake` object.
"""

snakemake = globals().get("snakemake")  # noqa: F821


def main():
    raise NotImplementedError(
        "prune_tree.py: not implemented. See the module docstring for the contract."
    )


if __name__ == "__main__":
    main()
