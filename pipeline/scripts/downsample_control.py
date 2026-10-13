#!/usr/bin/env python3
"""
Fig.S2 -- false-loss rate when Tier A is downsampled to Tier C depth.

DECISION THAT MUST NOT BE LOST:
Without this control the ancient losses in Fig.3 are indefensible. Match
the empirical Tier C depth distribution, not a single value.

Inputs / outputs are supplied by Snakemake via the `snakemake` object.
"""

snakemake = globals().get("snakemake")  # noqa: F821


def main():
    raise NotImplementedError(
        "downsample_control.py: not implemented. See the module docstring for the contract."
    )


if __name__ == "__main__":
    main()
