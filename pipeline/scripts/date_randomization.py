#!/usr/bin/env python3
"""
Date-randomization test, >=20 year-shuffled replicates.

DECISION THAT MUST NOT BE LOST:
Pass = the real 95% HPD does not overlap the shuffled distribution.

Inputs / outputs are supplied by Snakemake via the `snakemake` object.
"""

snakemake = globals().get("snakemake")  # noqa: F821


def main():
    raise NotImplementedError(
        "date_randomization.py: not implemented. See the module docstring for the contract."
    )


if __name__ == "__main__":
    main()
