#!/usr/bin/env python3
"""
Per-population root-to-tip regression (clock gate 3).

DECISION THAT MUST NOT BE LOST:
Run per population, never on the whole species -- the global Y. pestis
clock is unstable and a pooled regression hides it.

Inputs / outputs are supplied by Snakemake via the `snakemake` object.
"""

snakemake = globals().get("snakemake")  # noqa: F821


def main():
    raise NotImplementedError(
        "root_to_tip.py: not implemented. See the module docstring for the contract."
    )


if __name__ == "__main__":
    main()
