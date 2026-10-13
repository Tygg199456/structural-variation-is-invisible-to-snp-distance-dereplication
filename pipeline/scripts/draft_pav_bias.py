#!/usr/bin/env python3
"""
Fig.S3 -- spurious losses introduced by draft assemblies.

DECISION THAT MUST NOT BE LOST:
Use strains with BOTH complete and draft assemblies in NCBI; run each
through identical panaroo settings. This is the evidence for barring Tier
B from PAV, not an assertion.

Inputs / outputs are supplied by Snakemake via the `snakemake` object.
"""

snakemake = globals().get("snakemake")  # noqa: F821


def main():
    raise NotImplementedError(
        "draft_pav_bias.py: not implemented. See the module docstring for the contract."
    )


if __name__ == "__main__":
    main()
