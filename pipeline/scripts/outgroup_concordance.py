#!/usr/bin/env python3
"""
Is the mid-band meaningful, or is it just how much two outgroups differ
everywhere?

THE QUESTION
Measured on the recurrent blocks, the two outgroups disagree substantially:

    pgm 105 kb        IP32953 0.951   IP31758 0.627   delta -0.32
    1,105,000 +40 kb          0.751           0.659         -0.09
    795,000 +15 kb            0.593           0.555         -0.04

Their GLOBAL callable fractions are nearly identical (0.881 vs 0.879), so
this is not an overall quality difference. But before reading those deltas
as biology, we need the background: how much do two Y. pseudotuberculosis
strains differ, per window, across the whole chromosome?

    block deltas no larger than background  -> the mid band is interspecies
                                               /intraspecies divergence
                                               noise. Every polarity call is
                                               weak and the limitation must
                                               be stated more strongly.
    block deltas clearly larger             -> those regions really are
                                               special and the calls hold.

This also bears on the asymmetric consensus rule. It accepts a single
witness for ancestral presence, which is correct when the other outgroup
lost the region independently -- but if the two outgroups disagree broadly
and randomly, a "single witness" is as likely to be noise as signal.

Usage:
    python3 scripts/outgroup_concordance.py \
        results/snippy/GCA_000047365.1/snps.aligned.fa \
        results/snippy/GCF_000016945.1/snps.aligned.fa \
        --blocks results/qc/recurrent_block_polarity.tsv \
        --window 25000 --n-windows 1000

Expected magnitude: 1,000 background windows sampled; the block deltas (18
of them) compared against that distribution. A block at the 95th percentile
or beyond of background is "special".
"""
import argparse
import csv
import random
import sys
from pathlib import Path

MISSING = set("Nn-.*")


def norm_rep(rep):
    if "." in rep and rep.rsplit(".", 1)[1].isdigit():
        return rep.rsplit(".", 1)[0]
    return rep


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = norm_rep(line[1:].split()[0]), []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def called_frac(seq, start, length):
    s = max(start - 1, 0)
    e = min(s + length, len(seq))
    win = seq[s:e]
    if not win:
        return None
    return sum(1 for c in win if c not in MISSING) / len(win)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("og1")
    ap.add_argument("og2")
    ap.add_argument("--blocks", help="recurrent_block_polarity.tsv")
    ap.add_argument("--window", type=int, default=25000)
    ap.add_argument("--n-windows", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=20261011)
    ap.add_argument("--out", default="results/qc/outgroup_concordance.tsv")
    a = ap.parse_args()

    s1, s2 = read_fasta(a.og1), read_fasta(a.og2)
    shared = sorted(set(s1) & set(s2))
    if not shared:
        sys.exit(f"FATAL: no shared replicons.\n  {sorted(s1)[:4]}\n"
                 f"  {sorted(s2)[:4]}")

    rng = random.Random(a.seed)

    print("=" * 80)
    print("Outgroup concordance -- background versus the recurrent blocks")
    print("=" * 80)
    for rep in shared:
        f1 = called_frac(s1[rep], 1, len(s1[rep]))
        f2 = called_frac(s2[rep], 1, len(s2[rep]))
        print(f"  {rep:<16}og1 {f1:.3f}   og2 {f2:.3f}   "
              f"len {len(s1[rep]):,}")

    # ---- background: random windows on the longest shared replicon -------
    main_rep = max(shared, key=lambda r: len(s1[r]))
    L = min(len(s1[main_rep]), len(s2[main_rep]))
    deltas = []
    for _ in range(a.n_windows):
        st = rng.randrange(1, max(L - a.window, 2))
        c1 = called_frac(s1[main_rep], st, a.window)
        c2 = called_frac(s2[main_rep], st, a.window)
        if c1 is None or c2 is None:
            continue
        deltas.append(abs(c1 - c2))
    if not deltas:
        sys.exit("FATAL: no background windows sampled")

    deltas.sort()

    def pct(p):
        return deltas[min(int(p / 100 * len(deltas)), len(deltas) - 1)]

    print()
    print(f"background |og1 - og2| over {len(deltas)} windows of "
          f"{a.window:,} bp on {main_rep}")
    for p in (50, 75, 90, 95, 99):
        print(f"  p{p:<3} {pct(p):.4f}")
    print(f"  max  {deltas[-1]:.4f}")

    # ---- the blocks ------------------------------------------------------
    rows = []
    if a.blocks and Path(a.blocks).exists():
        print()
        print("recurrent blocks against that background")
        print(f"  {'replicon':<14}{'start':>11}{'len':>9}{'og1':>7}{'og2':>7}"
              f"{'delta':>8}{'bg_pctile':>11}  verdict")
        for r in csv.DictReader(open(a.blocks), delimiter="\t"):
            fr = r.get("outgroup_fracs_all") or ""
            vals = [float(x) for x in fr.split(";") if x.strip()]
            if len(vals) < 2:
                continue
            d = abs(vals[0] - vals[1])
            rank = sum(1 for x in deltas if x <= d) / len(deltas) * 100
            verdict = ("SPECIAL" if rank >= 95 else
                       "elevated" if rank >= 75 else "within background")
            rows.append({
                "replicon": r["replicon"], "start_bin": r["start_bin"],
                "length_bin": r["length_bin"],
                "og1": vals[0], "og2": vals[1], "delta": round(d, 4),
                "background_percentile": round(rank, 1),
                "polarity": r.get("polarity", ""),
                "n_loss_votes": r.get("n_loss_votes", ""),
                "verdict": verdict,
            })
            print(f"  {r['replicon']:<14}{int(r['start_bin']):>11,}"
                  f"{int(r['length_bin']):>9,}{vals[0]:>7.3f}{vals[1]:>7.3f}"
                  f"{d:>8.3f}{rank:>10.1f}%  {verdict}  "
                  f"[{r.get('polarity', '')}]")

        Path(a.out).parent.mkdir(parents=True, exist_ok=True)
        if rows:
            with open(a.out, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(rows[0]),
                                   delimiter="\t")
                w.writeheader()
                w.writerows(rows)

        n_special = sum(1 for r in rows if r["verdict"] == "SPECIAL")
        n_bg = sum(1 for r in rows if r["verdict"] == "within background")
        print()
        print(f"SPECIAL (>= p95 of background) : {n_special}")
        print(f"within background              : {n_bg}")
        if n_bg > len(rows) / 2:
            print()
            print("More than half the blocks sit within the background range.")
            print("The mid-band is then largely outgroup-to-outgroup variation,")
            print("not a property of these regions. Polarity calls resting on")
            print("a single witness (n_loss_votes = 1) should carry a")
            print("confidence flag in Table 2, and the limitation needs to")
            print("say that outgroup discordance is genome-wide rather than")
            print("local.")
        print(f"written: {a.out}")
    print("=" * 80)


if __name__ == "__main__":
    main()
