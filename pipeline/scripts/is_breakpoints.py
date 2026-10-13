#!/usr/bin/env python3
"""
Co-localisation of IS copies with loss/pseudogenisation breakpoints.

DECISION THAT MUST NOT BE LOST:
Test against a position-randomised null; raw overlap counts mean nothing
when IS elements occupy a large share of the genome.

Inputs / outputs are supplied by Snakemake via the `snakemake` object.
"""

snakemake = globals().get("snakemake")  # noqa: F821


def main():
    raise NotImplementedError(
        "is_breakpoints.py: not implemented. See the module docstring for the contract."
    )


if __name__ == "__main__":
    main()
