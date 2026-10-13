#!/usr/bin/env python3
"""
Fig.S4 -- all-3453 VeryFastTree vs deduplicated backbone.

DECISION THAT MUST NOT BE LOST:
Report Robinson-Foulds distance and the support-value distributions side
by side. This answers 'why discard 1800 genomes'.

Inputs / outputs are supplied by Snakemake via the `snakemake` object.
"""

snakemake = globals().get("snakemake")  # noqa: F821


def main():
    raise NotImplementedError(
        "figS4_topology_compare.py: not implemented. See the module docstring for the contract."
    )


if __name__ == "__main__":
    main()
