#!/usr/bin/env python3
"""
Assemble the ~350-taxon dating set.

DECISION THAT MUST NOT BE LOST:
Backbone supplements are chosen to MAXIMISE sampling-year span, not
assembly quality. Dating needs temporal spread; PAV needs Tier A.

WHY ANCHORS ARE *NOT* EXEMPT HERE (deliberate, not an oversight)
----------------------------------------------------------------
Whitelisted anchors bypass every other reduction point in the pipeline.
This one is the exception on purpose: a tip with no collection_date cannot
be placed on a time axis, and forcing an anchor in would mean fabricating a
date. An anchor that also has no date is therefore dropped like any other
undated tip. Outgroups are excluded separately (config.outgroup.*,
exclude_from_dating) -- they violate the molecular clock.

Inputs / outputs are supplied by Snakemake via the `snakemake` object.
"""

snakemake = globals().get("snakemake")  # noqa: F821


def main():
    raise NotImplementedError(
        "select_dating_taxa.py: not implemented. See the module docstring for the contract."
    )


if __name__ == "__main__":
    main()
