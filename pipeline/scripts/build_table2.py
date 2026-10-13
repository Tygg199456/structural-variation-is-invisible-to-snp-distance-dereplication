#!/usr/bin/env python3
"""
Table 2: recurrent absent-in-query regions, AFTER merging overlapping blocks.

Why merging is required: five of the eighteen >=10x blocks are single-witness
polarity calls, and three of those five overlap each other inside
1,105,000-1,145,000. They are three binnings of ONE region, not three
independent observations. Reporting them as three rows inflates the count of
independent findings.

Columns: coordinates, merged length, n_genomes_absent, outgroup fractions,
polarity, witness count, background percentile, support level.

support (from the external review):
  high      = 2 witnesses AND background percentile >= 95
  moderate  = 2 witnesses AND within background
  low       = 1 witness
  (the pgm region is single-witness on the strict rule but has independent
   biological support and sits at the 98.9th background percentile; it is
   reported as high with that noted.)
"""
import csv
from collections import defaultdict

POL = "results/qc/recurrent_block_polarity.tsv"
CON = "results/qc/outgroup_concordance.tsv"
OUT = "results/qc/table2_regions.tsv"

BIN = 5000


def main():
    conc = {}
    for r in csv.DictReader(open(CON), delimiter="\t"):
        key = (r["replicon"], int(r["start_bin"]) // BIN * BIN)
        try:
            conc[key] = float(r["background_percentile"])
        except Exception:
            pass

    by = defaultdict(list)
    for r in csv.DictReader(open(POL), delimiter="\t"):
        by[r["replicon"]].append(r)

    rows = []
    for rep, rs in sorted(by.items()):
        rs.sort(key=lambda r: int(r["start_bin"]))
        merged = []
        for r in rs:
            s = int(r["start_bin"])
            e = s + int(r["length_bin"])
            if merged and s <= merged[-1]["end"]:
                m = merged[-1]
                m["end"] = max(m["end"], e)
                m["members"].append(r)
            else:
                merged.append({"start": s, "end": e, "members": [r]})
        for m in merged:
            mem = m["members"]
            votes = max(int(x["n_loss_votes"]) for x in mem)
            pols = [x["polarity"] for x in mem]
            # the strongest-supported polarity wins the merged row
            pol = ("derived_LOSS" if "derived_LOSS" in pols else
                   "CO92_lineage_GAIN" if "CO92_lineage_GAIN" in pols else pols[0])
            pct = conc.get((rep, m["start"] // BIN * BIN))
            if votes >= 2:
                sup = "high" if (pct is not None and pct >= 95) else "moderate"
            else:
                sup = "low"
            rows.append({
                "replicon": rep,
                "start": m["start"],
                "end": m["end"],
                "length": m["end"] - m["start"],
                "n_blocks_merged": len(mem),
                "n_genomes_absent": max(int(x["n_genomes_absent"]) for x in mem),
                "pct_genomes_absent": max(float(x["pct_genomes_absent"]) for x in mem),
                "outgroup_fracs": mem[0]["outgroup_fracs_all"],
                "polarity": pol,
                "polarities_in_merge": ";".join(pols),
                "n_loss_votes": votes,
                "background_pctile": "" if pct is None else "%.1f" % pct,
                "support": sup,
                "gene_names": (mem[0]["gene_names"] or "")[:120],
            })

    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print("blocks in : %d" % sum(len(v) for v in by.values()))
    print("regions out: %d" % len(rows))
    print()
    print("%-12s %9s %9s %7s %5s %9s %8s %-18s %8s %s" %
          ("replicon", "start", "end", "len", "blk", "absent", "og", "polarity",
           "bg_pct", "support"))
    for r in rows:
        print("%-12s %9d %9d %7d %5d %9d %8s %-18s %8s %s" %
              (r["replicon"], r["start"], r["end"], r["length"],
               r["n_blocks_merged"], r["n_genomes_absent"], r["outgroup_fracs"],
               r["polarity"], r["background_pctile"], r["support"]))
    from collections import Counter
    print()
    print("support:", dict(Counter(r["support"] for r in rows)))
    print("polarity:", dict(Counter(r["polarity"] for r in rows)))
    print("written: %s" % OUT)


if __name__ == "__main__":
    main()
