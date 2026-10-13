#!/usr/bin/env python3
"""
Fig.5 -- lineage x medium growth heatmap, auxotrophy onset.

DECISION THAT MUST NOT BE LOST:
Mark lineages lacking a calibrated clock; their onset times are relative,
not calendar dates.

Inputs / outputs are supplied by Snakemake via the `snakemake` object.
"""

snakemake = globals().get("snakemake")  # noqa: F821


def main():
    raise NotImplementedError(
        "fig5_fba.py: not implemented. See the module docstring for the contract."
    )


if __name__ == "__main__":
    main()
