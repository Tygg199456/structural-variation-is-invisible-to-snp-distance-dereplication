#!/usr/bin/env python3
"""
Breadth-based gene presence for ancient genomes.

DECISION THAT MUST NOT BE LOST:
Presence = breadth >= 90% of gene length at >= 3x. Emit
presence/absence/UNKNOWN; never emit a frameshift or pseudogene call.

Inputs / outputs are supplied by Snakemake via the `snakemake` object.
"""

snakemake = globals().get("snakemake")  # noqa: F821


def main():
    raise NotImplementedError(
        "coarse_pav.py: not implemented. See the module docstring for the contract."
    )


if __name__ == "__main__":
    main()
