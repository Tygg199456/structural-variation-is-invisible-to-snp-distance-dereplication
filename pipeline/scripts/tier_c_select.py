#!/usr/bin/env python3
"""
Select ancient libraries passing depth + damage authentication.

DECISION THAT MUST NOT BE LOST:
UDG-treated libraries carry little terminal deamination and will fail the
damage gate; read library prep from AncientMetagenomeDir and exempt them
rather than silently dropping real data.

Inputs / outputs are supplied by Snakemake via the `snakemake` object.
"""

snakemake = globals().get("snakemake")  # noqa: F821


def main():
    raise NotImplementedError(
        "tier_c_select.py: not implemented. See the module docstring for the contract."
    )


if __name__ == "__main__":
    main()
