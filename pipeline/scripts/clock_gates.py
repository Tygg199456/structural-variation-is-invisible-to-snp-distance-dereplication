#!/usr/bin/env python3
"""
Apply the four clock gates and decide BEAST vs LSD2 per population.

DECISION THAT MUST NOT BE LOST:
All four must pass: n>=30 dated tips, span>=50 y, R2>=0.30 with positive
slope, and a passing date-randomization test.

Inputs / outputs are supplied by Snakemake via the `snakemake` object.
"""

snakemake = globals().get("snakemake")  # noqa: F821


def main():
    raise NotImplementedError(
        "clock_gates.py: not implemented. See the module docstring for the contract."
    )


if __name__ == "__main__":
    main()
