#!/usr/bin/env python3
"""
Fig.6 -- loss events per branch vs node age (step vs continuous).

DECISION THAT MUST NOT BE LOST:
Plot gate-failing populations on relative time and flag them visually; do
not mix calendar and relative axes silently.

Inputs / outputs are supplied by Snakemake via the `snakemake` object.
"""

snakemake = globals().get("snakemake")  # noqa: F821


def main():
    raise NotImplementedError(
        "fig6_decay_rate.py: not implemented. See the module docstring for the contract."
    )


if __name__ == "__main__":
    main()
