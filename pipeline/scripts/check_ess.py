#!/usr/bin/env python3
"""
Parse BEAST logs, assert ESS >= threshold for every parameter.

DECISION THAT MUST NOT BE LOST:
Exit non-zero on failure; a converged-looking run with low ESS on one
parameter is still unusable.

Inputs / outputs are supplied by Snakemake via the `snakemake` object.
"""

snakemake = globals().get("snakemake")  # noqa: F821


def main():
    raise NotImplementedError(
        "check_ess.py: not implemented. See the module docstring for the contract."
    )


if __name__ == "__main__":
    main()
