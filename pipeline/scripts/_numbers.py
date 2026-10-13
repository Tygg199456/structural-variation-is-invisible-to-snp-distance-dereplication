#!/usr/bin/env python3
"""_numbers.py — shared producer for results/numbers.json.

External review 2026-10-13: every analysis step should append the numbers it
computed to results/numbers.json, so that a later audit is a diff instead of a
re-run, and so that a number in the manuscript with no producing record is
visible as a gap.  One JSON object per line (JSON Lines):

    {"ts": "...Z", "git": "abc1234", "script": "s3_specificity.py",
     "numbers": {"pgm_arm1_labnamed": 9, "pgm_arm1_n": 41, ...}}

Deliberately append-only and never fatal: recording a number must not change
the result of the analysis that produced it.
"""
import datetime
import json
import os
import subprocess

DEFAULT = "results/numbers.json"


def _git():
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return "no-git"


def record(script, numbers, path=DEFAULT):
    """Append {ts, git, script, numbers} to *path* as one JSON line.

    `numbers` is any JSON-serialisable mapping.  Returns the path written, or
    None if the write failed (a failure to record must not abort a run).
    """
    rec = {
        "ts": datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"),
        "git": _git(),
        "script": os.path.basename(script),
        "numbers": numbers,
    }
    try:
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        with open(path, "a") as fh:
            fh.write(json.dumps(rec, sort_keys=True) + "\n")
        return path
    except Exception as exc:                      # never fail the analysis
        print("[_numbers] WARNING: could not record numbers: %r" % exc)
        return None
