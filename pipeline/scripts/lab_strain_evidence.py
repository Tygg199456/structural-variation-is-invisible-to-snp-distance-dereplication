#!/usr/bin/env python3
"""
Evidence for the laboratory-versus-natural framing of paper 1.

CORRECTION (2026-10-11): the first version of this script read
`dnadiff_aligned_pct` from results/qc/excluded_with_reason.tsv. That column is
0.0 for ALL 210 rows -- it is an unpopulated placeholder, not a measurement.
Testing `< 95` on it therefore classified every excluded genome as a
dnadiff failure and produced a nonsense geography table. The real numbers are
in results/qc/dnadiff/<acc>.report, on the `[REF] AlignedBases` line. This is
the same failure shape the paper is about: a metric that is missing reads as a
metric that is bad, and nothing raises an error.

Produces:

  A. Laboratory / vaccine / reference strain inventory, verified in our own
     data (project rule: never assert a strain's properties from memory), with
     each strain's deletion fingerprint.
  B. Geography of the three reference-relative quality metrics, counted rather
     than asserted: CheckM completeness, dnadiff REF-side AlignedBases, and
     called_frac.

Usage: python3 scripts/lab_strain_evidence.py
"""
import csv
import glob
import os
import re
from collections import Counter, defaultdict

PROFILE = "results/qc/deletion_profile_all.tsv"
META = "results/qc/assembly_metadata.tsv"
EXCL = "results/qc/excluded_with_reason.tsv"
LINQC = "results/qc/lineage_stratified_qc.tsv"
DNADIFF = "results/qc/dnadiff"
OUT_A = "results/qc/lab_strain_inventory.tsv"
OUT_B = "results/qc/quality_metric_geography.tsv"

REPL = {"NC_003143": "chromosome", "NC_003134": "pMT1",
        "NC_003131": "pCD1", "NC_003132": "pPCP1"}

NAMED = ["EV NIIEG", "EV76", "KIM10", "A1122", "CO92", "C-792", "C-783",
         "C-830", "C-25", "M-1482", "Nepal516", "Pestoides", "Antiqua",
         "Angola", "Java9", "Harbin35"]

EURASIA = re.compile(
    r"china|russia|kazakhstan|kyrgyzstan|mongolia|uzbekistan|tajikistan|"
    r"turkmenistan|azerbaijan|armenia|georgia|ukraine|belarus|moldova|"
    r"india|nepal|myanmar|vietnam|iran|iraq|turkey|madagascar|"
    r"africa|asia|europe", re.I)
AMERICAS = re.compile(r"usa|united states|canada|mexico|brazil|peru|"
                      r"argentina|bolivia|ecuador|south america|north america", re.I)


def region_of(geo):
    geo = (geo or "").strip()
    if not geo:
        return "unknown"
    if EURASIA.search(geo):
        return "Eurasia"
    if AMERICAS.search(geo):
        return "Americas"
    return "other"


def read_dnadiff(basedir):
    """[REF] AlignedBases percentage from every *.report."""
    out = {}
    for path in glob.glob(os.path.join(basedir, "*.report")):
        acc = os.path.basename(path)[:-len(".report")]
        with open(path, errors="replace") as fh:
            for line in fh:
                if line.startswith("AlignedBases"):
                    m = re.findall(r"([0-9.]+)%", line)
                    if m:
                        out[acc] = float(m[0])       # REF-side percentage
                    break
    return out


def main():
    prof, meta = {}, {}
    for r in csv.DictReader(open(PROFILE), delimiter="\t"):
        prof[r["accession"].strip()] = r
    for r in csv.DictReader(open(META), delimiter="\t"):
        meta[r["accession"].strip()] = r
    for a, r in list(meta.items()):
        p = (r.get("paired_accession") or "").strip()
        if p:
            meta.setdefault(p, r)

    # ---------------- A. laboratory strain inventory -------------------
    print("=" * 104)
    print("A. LAB/VACCINE/REFERENCE STRAINS -- verified in our own data")
    print("=" * 104)
    print("%-11s %-17s %-16s %-14s %s" %
          ("query", "accession", "strain", "replicons-", ">=10 kb blocks"))
    rows_a = []
    for q in NAMED:
        hits = sorted(a for a, r in meta.items()
                      if q.lower() in (r.get("strain") or "").lower() and a in prof)
        if not hits:
            print("%-11s %s" % (q, "(not in the 3,441-genome profile)"))
            continue
        for a in hits[:3]:
            p = prof[a]
            absent = (p.get("replicons_absent") or "").strip()
            blk = (p.get("big_blocks_5kb") or "").strip()
            blk_s = ",".join(
                "%s:%dk" % (REPL.get(b.split(":")[0], b.split(":")[0]),
                            int(b.split(":")[2]) // 1000)
                for b in blk.split(";") if len(b.split(":")) == 3)[:70]
            print("%-11s %-17s %-16s %-14s %s" % (
                q, a, (meta[a].get("strain") or "")[:16],
                ",".join(REPL.get(x, x) for x in absent.split(";") if x) or "-",
                blk_s or "-"))
            rows_a.append({"query": q, "accession": a,
                           "strain_field": meta[a].get("strain", ""),
                           "geo_loc_name": meta[a].get("geo_loc_name", ""),
                           "replicons_absent": absent, "big_blocks_5kb": blk})
    if rows_a:
        with open(OUT_A, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows_a[0]), delimiter="\t")
            w.writeheader()
            w.writerows(rows_a)
        print("\nwritten: %s" % OUT_A)

    # ---------------- B. geography of the three metrics ----------------
    print()
    print("=" * 104)
    print("B. WHERE THE THREE REFERENCE-RELATIVE METRICS FALL")
    print("=" * 104)

    # metric 1: CheckM completeness < 98   (real values in excluded_with_reason)
    m1 = []
    for r in csv.DictReader(open(EXCL), delimiter="\t"):
        try:
            v = float(r["checkm_completeness_COVARIATE_ONLY"])
        except (KeyError, ValueError):
            continue
        if v < 98:
            m1.append((r["accession"].strip(), v,
                       (r.get("geo_loc_name") or "").strip(),
                       (r.get("strain") or "").strip()))

    # metric 2: dnadiff REF AlignedBases < 95
    dn = read_dnadiff(DNADIFF)
    m2 = [(a, v, (meta.get(a, {}).get("geo_loc_name") or "").strip(),
           (meta.get(a, {}).get("strain") or "").strip())
          for a, v in dn.items() if v < 95]

    # metric 3: called_frac < 0.859
    m3 = []
    for r in csv.DictReader(open(LINQC), delimiter="\t"):
        try:
            v = float(r["called_frac"])
        except (KeyError, ValueError):
            continue
        if v < 0.859:
            m3.append((r["accession"].strip(), v,
                       (r.get("geo_loc_name") or "").strip(),
                       (r.get("strain") or "").strip()))

    rows_b = []
    for name, rows, unit in [("checkm_lt98", m1, "completeness"),
                             ("dnadiff_ref_lt95", m2, "% of ref aligned"),
                             ("called_frac_lt0.859", m3, "called_frac")]:
        c = Counter(region_of(g) for _, _, g, _ in rows)
        n = len(rows)
        parts = ", ".join("%s %d (%.0f%%)" % (k, v, 100 * v / n)
                          for k, v in c.most_common())
        print("\n  %-20s n=%-4d %s" % (name, n, parts))
        for acc, v, geo, strain in sorted(rows, key=lambda x: x[1])[:12]:
            print("      %-17s %-9s %-8s %-8s %s" %
                  (acc, region_of(geo), strain[:8], "%.4g" % v, geo[:52]))
            rows_b.append({"metric": name, "accession": acc, "value": v,
                           "unit": unit, "region": region_of(geo),
                           "geo_loc_name": geo, "strain": strain})
        if n > 12:
            print("      ... and %d more" % (n - 12))

    with open(OUT_B, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows_b[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows_b)
    print("\nwritten: %s" % OUT_B)


if __name__ == "__main__":
    main()
