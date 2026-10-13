#!/usr/bin/env python3
"""C.4 single-source decomposition for SY1011 (handoff §2.2 / §2.3).

Everything derives from ONE source: the M-to-M coordinate list of the dnadiff
run, results/qc/dnadiff/<acc>.mcoords.  The complement of the union of aligned
reference intervals is the unaligned set; it is partitioned into

    del_replicon : unaligned bases falling on a reference replicon the assembly
                   does NOT carry (union coverage < CARRIED_FRAC of that replicon)
    del_block    : unaligned contiguous runs >= BLOCK bp on carried replicons
    residual     : every other unaligned base

The identity  del_replicon + del_block + residual == unaligned  is asserted per
genome (it is a partition, not a sum of independently-derived terms).

The unaligned total is asserted against the dnadiff .report UnalignedBases --
the number behind the published 93.00% figure -- so the two axes of the
"x+y=5" figure come from the same tool.

WHY M-to-M AND NOT 1-to-1 (this corrects the handoff's §2.3 prescription):
the .1coords (delta-filter -1) complement is NOT the report's unaligned figure.
1-to-1 filtering drops repeat alignments, so on C-792 the .1coords complement
is 414,093 bp while the report says 338,039 bp (off by 76,054 bp = 1.57 pp).
The .mcoords (M-to-M = the full alignment set) complement reproduces the report
to within a few bp.  Using .1coords would have introduced a bigger error than
every quantity C.4 was arguing about.

Outputs
  results/qc/c4_decomposition.tsv     per-genome rows
  stdout                              aggregate + the reported 5 low-REF genomes
"""
import csv, glob, os, re, sys, collections
import os

ROOT = os.environ.get("PLAGUE_DECAY_ROOT",
                       os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.chdir(ROOT)
DD = "results/qc/dnadiff"
BLOCK = 10000          # bp, "block" threshold from handoff §2.2
CARRIED_FRAC = 0.50    # replicon counted as carried if >= half its length aligns
TOL = 10               # bp tolerance vs the dnadiff report


def merge(iv):
    iv = sorted(iv)
    m = []
    for s, e in iv:
        if m and s <= m[-1][1] + 1:
            m[-1] = (m[-1][0], max(m[-1][1], e))
        else:
            m.append((s, e))
    return m


def complement(m, L):
    out = []
    pos = 1
    for s, e in m:
        if s > pos:
            out.append((pos, s - 1))
        pos = max(pos, e + 1)
    if pos <= L:
        out.append((pos, L))
    return out


def read_mcoords(path):
    by = collections.defaultdict(list)
    reflen = {}
    nrows = 0
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 13:
            continue
        try:
            s1, e1 = int(f[0]), int(f[1])
            reflen[f[11]] = int(f[7])
        except ValueError:
            raise SystemExit("bad mcoords row in %s: %r" % (path, line[:80]))
        by[f[11]].append((min(s1, e1), max(s1, e1)))
        nrows += 1
    assert nrows > 0, "no alignment rows in %s" % path
    return by, reflen


def report_unaligned(path):
    """First UnalignedBases number (REF column) from a dnadiff .report."""
    val = None
    for line in open(path):
        if line.startswith("UnalignedBases"):
            toks = line.split()
            # "UnalignedBases  338039(6.9989%)  3325(0.0742%)"
            m = re.match(r"(\d+)", toks[1])
            assert m, "cannot parse %r" % line
            val = int(m.group(1))
            break
    assert val is not None, "no UnalignedBases in %s" % path
    return val


rows = []
mismatch = []
for mc in sorted(glob.glob(DD + "/*.mcoords")):
    acc = os.path.basename(mc)[: -len(".mcoords")]
    rep = DD + "/" + acc + ".report"
    if not os.path.exists(rep):
        mismatch.append((acc, "NO_REPORT", 0, 0))
        continue
    by, reflen = read_mcoords(mc)
    unaligned = 0
    del_rep = 0
    del_blk = 0
    resid = 0
    absent_tags = []
    carried_tags = []
    for tag in sorted(reflen):
        L = reflen[tag]
        cov = merge(by.get(tag, []))
        covlen = sum(e - s + 1 for s, e in cov)
        uncov = complement(cov, L)
        u = sum(e - s + 1 for s, e in uncov)
        unaligned += u
        if L > 0 and covlen >= CARRIED_FRAC * L:
            carried_tags.append(tag)
            for s, e in uncov:
                if e - s + 1 >= BLOCK:
                    del_blk += e - s + 1
                else:
                    resid += e - s + 1
        else:
            absent_tags.append(tag)
            del_rep += u
    assert del_rep + del_blk + resid == unaligned, (
        "partition broken for %s: %d+%d+%d != %d"
        % (acc, del_rep, del_blk, resid, unaligned)
    )
    rep_un = report_unaligned(rep)
    diff = unaligned - rep_un
    if abs(diff) > TOL:
        mismatch.append((acc, "MISMATCH", unaligned, rep_un))
    rows.append((acc, unaligned, rep_un, diff, del_rep, del_blk, resid,
                 ";".join(absent_tags), len(carried_tags)))

with open("results/qc/c4_decomposition.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["accession", "unaligned_mcoords", "unaligned_report", "diff",
                "del_replicon", "del_block", "residual", "absent_replicons",
                "n_carried_replicons"])
    for r in rows:
        w.writerow(r)

print("genomes processed            : %d" % len(rows))
print("identity del_rep+blk+residual: OK for all %d" % len(rows))
print("report mismatches > %d bp     : %d" % (TOL, len(mismatch)))
for m in mismatch[:20]:
    print("   MISMATCH %s %s %s %s" % m)
print()
diffs = [abs(r[3]) for r in rows]
if diffs:
    print("|mcoords - report| max=%d  mean=%.2f" % (max(diffs), sum(diffs) / len(diffs)))
print()
print("=== the 5 lowest REF-side genomes (report dnadiff%) ===")
sc = sorted(rows, key=lambda r: -r[2])
# need per-genome total; use the fixed reference total
REF_TOTAL = 4829855
for r in sc[:8]:
    print("  %-18s report_unaligned=%7d  (%.2f%% aligned)  mcoords=%7d  diff=%+d"
          % (r[0], r[2], 100.0 * (REF_TOTAL - r[2]) / REF_TOTAL, r[1], r[3]))
print()
print("=== C-792 (GCF_015337085.2) ===")
for r in rows:
    if r[0] == "GCF_015337085.2":
        print("  unaligned=%d (report %d, diff %+d)" % (r[1], r[2], r[3]))
        print("  del_replicon=%d  del_block=%d  residual=%d  sum=%d"
              % (r[4], r[5], r[6], r[4] + r[5] + r[6]))
        print("  absent replicons: %s" % r[7])
print()
print("WROTE results/qc/c4_decomposition.tsv")
