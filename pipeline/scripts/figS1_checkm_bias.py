#!/usr/bin/env python3
"""
Fig.S1 -- CheckM completeness vs pseudogene count, with regression.

DECISION THAT MUST NOT BE LOST:
A significant negative slope is the published justification for not
filtering on completeness. Report slope, R2 and p in the table.

Inputs / outputs are supplied by Snakemake via the `snakemake` object.
"""

snakemake = globals().get("snakemake")  # noqa: F821


def main():
    raise NotImplementedError(
        "figS1_checkm_bias.py: not implemented. See the module docstring for the contract."
    )


if __name__ == "__main__":
    main()
