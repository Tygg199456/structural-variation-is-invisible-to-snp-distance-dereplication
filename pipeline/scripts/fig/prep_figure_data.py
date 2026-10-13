#!/usr/bin/env python3
"""
Tidy every number paper 1's figures need.  R scripts read only the CSVs written
here; no .npy/.npz/GFF access from R.

Design rules (external review 2026-10-11):
  * a failed step RAISES -- it is never printed and skipped, because R would
    then read a stale CSV and draw a wrong figure without an error;
  * every CSV is written to .tmp and os.replace'd, so a failure leaves no
    half file;
  * every CSV carries a provenance header line (UTC time, source files, git).

Outputs (results/figures/):  see the CSV list in RESULTS.md / the docstring of
each function.

Run:  python3 scripts/fig/prep_figure_data.py
"""
import csv, os, re, subprocess, sys, datetime, random
from collections import Counter, defaultdict
import numpy as np

OUT = "results/figures"
SENT = 10 ** 6
REF_TOTAL = 4829855
os.makedirs(OUT, exist_ok=True)
REPL = {"NC_003143": "chromosome", "NC_003134": "pMT1",
        "NC_003131": "pCD1", "NC_003132": "pPCP1"}
PROV = []


def _git():
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"],
                                       stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return "no-git"


GIT = _git()
NOW = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
CO92_ANNOT = "results/annotate/GCF_000009065.1/GCF_000009065.1.gff3"


def w(name, header, rows, sources=()):
    path = os.path.join(OUT, name)
    tmp = path + ".tmp"
    with open(tmp, "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(header)
        wr.writerows(rows)
    os.replace(tmp, path)
    with open(path + ".prov", "w") as ph:
        ph.write("generated %s UTC | git %s | sources: %s\n"
                 % (NOW, GIT, "; ".join(sources)))
    PROV.append((name, len(rows)))
    print("  wrote %-30s %d rows" % (name, len(rows)))


def rd(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def norm(a):
    return re.sub(r"\.\d+$", "", (a or "").strip())


def reads_fasta_run(path):
    """Maximal [Nn-.*] runs on a snps.aligned.fa, per replicon."""
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    out = []
    for r, s in seqs.items():
        for m in re.finditer(r"[Nn\-.]+", s):
            out.append((r, m.start() + 1, m.end(), m.end() - m.start()))
    return out


META = {}
for r in rd("results/qc/assembly_metadata.tsv"):
    META[r["accession"]] = r
    META.setdefault(norm(r["accession"]), r)
    p = (r.get("paired_accession") or "").strip()
    if p:
        META.setdefault(p, r)
        META.setdefault(norm(p), r)


def meta(a):
    return META.get(a, META.get(norm(a), {}))


# ============================================================ Fig.1
def fig1a():
    """Called fraction per genome. NOTE: the denominator is the union variant
    site set, NOT the chromosome length -- so the 0.955 chromosome-length
    self-comparison ceiling is a different quantity and is deliberately not
    drawn on this axis."""
    rows = []
    for tier in ("tierA", "tierB"):
        d = np.load("results/dedup/genotypes_%s.npz" % tier, allow_pickle=True)
        for s, v in zip(d["samples"], d["called_frac"]):
            rows.append([tier, str(s), "%.6f" % float(v)])
    w("fig1a_called_frac.csv", ["tier", "sample", "called_frac"], rows,
      ["results/dedup/genotypes_tier{A,B}.npz"])


def fig1b():
    D = np.load("results/dedup/snp_distance_raw_tierB.npy").astype(float)
    iu = np.triu_indices(D.shape[0], k=1)
    v = D[iu]
    v = v[v < SENT]
    edges = np.logspace(np.log10(0.5), np.log10(v.max() + 1), 61)
    cnt, _ = np.histogram(np.maximum(v, 0.5), bins=edges)
    rows = [["%.3f" % edges[i], "%.3f" % edges[i + 1], int(cnt[i])]
            for i in range(len(cnt))]
    w("fig1b_distance_hist.csv", ["bin_lo", "bin_hi", "n_pairs"], rows,
      ["results/dedup/snp_distance_raw_tierB.npy"])
    w("fig1b_scalars.csv", ["statistic", "value"], [
        ["n_pairs", int(v.size)], ["n_zero", int((v == 0).sum())],
        ["n_le_25", int((v <= 25).sum())],
        ["median", int(np.median(v))], ["p95", int(np.percentile(v, 95))],
        ["max", int(v.max())]],
      ["results/dedup/snp_distance_raw_tierB.npy"])
    pairs = [("Angola", "GCF_000018805.1", "GCF_000834845.1"),
             ("Pestoides F", "GCF_000016445.1", "GCF_000834315.1"),
             ("Antiqua", "GCF_000013825.1", "GCF_000834825.1"),
             ("CO92", "GCF_000009065.1", "GCF_001293415.1"),
             ("A1122", "GCF_000222975.1", "GCF_000834755.1")]
    mats = {}
    for tier in ("tierA", "tierB"):
        mats[tier] = (np.load("results/dedup/snp_distance_raw_%s.npy" % tier).astype(float),
                      {a: i for i, a in enumerate(
                          l.strip() for l in open(
                              "results/dedup/snp_distance_raw_%s.order.txt" % tier)
                          if l.strip())})
    srows = []
    for name, a, b in pairs:
        r = [name, a, b]
        for tier in ("tierA", "tierB"):
            Dt, it = mats[tier]
            r.append("%.0f" % Dt[it[a], it[b]] if a in it and b in it else "NA")
        srows.append(r)
    w("fig1b_same_strain.csv",
      ["strain", "acc_a", "acc_b", "raw_distance_tierA", "raw_distance_tierB"],
      srows, ["results/dedup/snp_distance_raw_tier{A,B}.npy"])


def fig1c():
    runs = reads_fasta_run("results/snippy/GCF_000009065.1/snps.aligned.fa")
    rows = [[i, r[3], r[0], r[1], r[2]] for i, r in enumerate(sorted(runs, key=lambda x: -x[3]))]
    w("fig1c_self_blocks.csv", ["block_index", "length_bp", "replicon", "start", "end"],
      rows, ["results/snippy/GCF_000009065.1/snps.aligned.fa"])
    big = [r for r in rows if r[1] >= 500]
    print("      self-comparison: %d runs, %d >=500 bp, max %d, %d >=10 kb"
          % (len(rows), len(big), max(r[1] for r in rows), sum(1 for r in rows if r[1] >= 10000)))


def fig1d():
    D = np.load("results/dedup/snp_distance_raw_tierB.npy").astype(float)
    order = [l.strip() for l in
             open("results/dedup/snp_distance_raw_tierB.order.txt") if l.strip()]
    idx = {a: i for i, a in enumerate(order)}
    target = None
    for r in rd("results/dedup/tierB_clusters.tsv"):
        if int(r.get("size", 0)) >= 400:
            target = r
            break
    if target is None:
        raise RuntimeError("no >=400 cluster found in tierB_clusters.tsv")
    mem = [m for m in target["members"].split(",") if m in idx]
    sub = D[np.ix_([idx[a] for a in mem], [idx[a] for a in mem])]
    d = sub[np.triu_indices(len(mem), 1)]
    d = d[d < SENT]
    rows = [["%d" % k, int((d == k).sum())] for k in range(0, int(d.max()) + 1)]
    w("fig1d_cluster_internal.csv", ["distance_snp", "n_pairs"], rows,
      ["results/dedup/snp_distance_raw_tierB.npy", "results/dedup/tierB_clusters.tsv"])
    yrs = []
    for a in mem:
        y = (meta(a).get("collection_date") or "")[:4]
        if y.isdigit():
            yrs.append(int(y))
    w("fig1d_cluster_years.csv", ["statistic", "value"],
      [["n_members", len(mem)], ["year_min", min(yrs)], ["year_max", max(yrs)],
       ["n_dated", len(yrs)], ["p50", int(np.percentile(d, 50))],
       ["p95", int(np.percentile(d, 95))], ["max", int(d.max())]],
      ["results/dedup/tierB_clusters.tsv", "results/qc/assembly_metadata.tsv"])


# ============================================================ Fig.2
def fig2a():
    rb = {norm(r["accession"]): r for r in rd("results/qc/reference_bias.tsv")}
    c4 = {norm(r["accession"]): 100.0 * (1 - int(r["unaligned_report"]) / REF_TOTAL)
          for r in rd("results/qc/c4_decomposition.tsv")}
    ctry = {}
    for a in rb:
        ctry[a] = (meta(a).get("geo_loc_name") or "").split(":")[0].strip() or "unknown"
    keep = {c for c, n in Counter(ctry.values()).items() if n >= 8}
    rows = []
    for a, r in rb.items():
        reg = ctry[a] if ctry[a] in keep else "other"
        try:
            ck = float(r["checkm_completeness"])
        except ValueError:
            ck = float("nan")
        try:
            cf = float(r["called_frac"])
        except ValueError:
            cf = float("nan")
        dd = c4.get(a, float("nan"))
        for m, v, thr in (("checkm", ck, 98.0), ("called_frac", cf, 0.859),
                          ("dnadiff_ref", dd, 95.0)):
            flag = "" if v != v else int(v < thr)
            rows.append([a, reg, m, "%s" % v, flag])
    w("fig2a_metric_flags.csv",
      ["accession", "region", "metric", "value", "flagged"], rows,
      ["results/qc/reference_bias.tsv", "results/qc/assembly_metadata.tsv"])


def fig2b():
    rows = []
    for r in rd("results/qc/reference_bias.tsv"):
        a = norm(r["accession"])
        lab = (r.get("strain") or "").strip() if float(r.get("residual_z") or 0) < -4 else ""
        anchor = r["accession"] if r["accession"] in (
            "GCF_000013805.1", "GCA_000006645.1", "GCF_000009065.1",
            "GCF_000016445.1", "GCF_000013825.1", "GCF_000018805.1",
            "GCF_000222975.1", "GCF_000007885.1") else ""
        rows.append([r["accession"], r["dist_to_reference"], r["called_frac"],
                     r["residual_z"], r["flag"], anchor])
    w("fig2b_regression.csv",
      ["accession", "distance_to_CO92", "called_frac", "resid_z", "flag", "anchor_label"],
      rows, ["results/qc/reference_bias.tsv"])


def _dnadiff_pct(acc):
    p = "results/qc/dnadiff/%s.report" % acc
    if not os.path.exists(p):
        return (float("nan"), float("nan"))
    ref = qry = float("nan")
    for line in open(p):
        if line.startswith("AlignedBases"):
            f = line.split()
            ref = float(f[1][f[1].find("(") + 1:f[1].find("%")])
            qry = float(f[2][f[2].find("(") + 1:f[2].find("%")])
            break
    return (ref, qry)


def fig2c():
    ts = {norm(r["accession"]): r for r in rd("results/qc/replicon_three_state_seqreport.tsv")}
    dmap = {}
    for fn in os.listdir("results/qc/dnadiff"):
        if fn.endswith(".report"):
            dmap[norm(fn[:-len(".report")])] = fn[:-len(".report")]
    lab5 = {"GCF_015337085": "C-792", "GCF_015337445": "KM 567",
            "GCF_015338205": "C-830", "GCF_015337285": "C-783",
            "GCF_000590535": "EV NIIEG"}
    rows = []
    for a in sorted(ts):
        v = dmap.get(a)
        if v is None:
            continue
        ref, qry = _dnadiff_pct2(v)
        if ref != ref:
            continue
        npl = int(ts[a]["n_plasmid_seq"])
        rows.append([v, "%.2f" % ref, "%.2f" % qry, "yes" if npl > 0 else "no",
                     lab5.get(a, "")])
    w("fig2c_dnadiff.csv",
      ["accession", "ref_pct", "qry_pct", "plasmid_submitted", "strain_label"],
      rows, ["results/qc/dnadiff/*.report", "results/qc/replicon_three_state_seqreport.tsv"])


def _dnadiff_pct2(acc):
    p = "results/qc/dnadiff/%s.report" % acc
    ref = qry = float("nan")
    if not os.path.exists(p):
        return (ref, qry)
    for line in open(p):
        if line.startswith("AlignedBases"):
            f = line.split()
            ref = float(f[1][f[1].find("(") + 1:f[1].find("%")])
            qry = float(f[2][f[2].find("(") + 1:f[2].find("%")])
            break
    return (ref, qry)


def fig2d():
    anch = {"GCF_000013825.1": "Antiqua", "GCF_000018805.1": "Angola",
            "GCF_000016445.1": "Pestoides F"}
    rows = []
    for r in rd("results/qc/reference_bias.tsv"):
        rows.append([r["accession"], r["checkm_completeness"],
                     anch.get(r["accession"], "")])
    w("fig2d_checkm.csv", ["accession", "completeness", "anchor_label"], rows,
      ["results/qc/reference_bias.tsv"])


# ============================================================ Fig.3
def fig3a():
    rows = []
    for r in rd("results/qc/uncalled_blocks.tsv"):
        rows.append([r["accession"], r["replicon"], int(r["start"]), int(r["length"])])
    w("fig3a_block_lengths.csv", ["accession", "replicon", "start", "length_bp"],
      rows, ["results/qc/uncalled_blocks.tsv"])


def fig3b():
    cnt = Counter()
    for r in rd("results/qc/uncalled_blocks_all.tsv"):
        if int(r["length"]) >= 10000:
            cnt[(r["replicon"], int(r["start"]), int(r["length"]))] += 1
    rows = [[k[0], k[1], k[2], v] for k, v in sorted(cnt.items(), key=lambda x: -x[1])]
    w("fig3b_recurrence.csv", ["replicon", "start", "length_bp", "n_genomes"],
      rows, ["results/qc/uncalled_blocks_all.tsv"])
    prof = {r["accession"]: r for r in rd("results/qc/deletion_profile_all.tsv")}
    lv = Counter()
    tot = Counter()
    for a, r in prof.items():
        level = (meta(a).get("assembly_level") or "unknown").lower()
        level = level.replace(" genome", "").replace("complete genome", "complete")
        key = ("complete" if level.startswith("complete") else
               "chromosome" if level.startswith("chromosome") else
               "scaffold" if level.startswith("scaffold") else "contig")
        tot[key] += 1
        if any(int(b.split(":")[2]) >= 10000 for b in (r.get("big_blocks_5kb") or "").split(";")
               if len(b.split(":")) == 3):
            lv[key] += 1
    rows = [[k, tot[k], lv[k], "%.4f" % (lv[k] / tot[k])]
            for k in ("complete", "chromosome", "scaffold", "contig") if tot[k]]
    w("fig3b_level_frac.csv", ["level", "n_total", "n_with_block", "frac"], rows,
      ["results/qc/deletion_profile_all.tsv", "results/qc/assembly_metadata.tsv"])


def fig3c():
    prof = {r["accession"]: r for r in rd("results/qc/deletion_profile_all.tsv")}
    rows = []
    for r in rd("results/qc/uncalled_blocks.tsv"):
        if r["replicon"] != "NC_003143":
            continue
        s, e = int(r["start"]), int(r["end"])
        if 2130000 <= s <= 2145000 and (e - s) >= 100000 and e >= 2230000:
            rows.append([meta(r["accession"]).get("strain") or r["accession"],
                         r["accession"], s, e, e - s])
    rows.sort(key=lambda x: x[2])
    w("fig3c_pgm_align.csv", ["strain", "accession", "start", "end", "length_bp"],
      rows, ["results/qc/uncalled_blocks.tsv"])
    print("      pgm exact runs: %d assemblies" % len(rows))
    grp = [("ybt", re.compile(r"ybt[A-Z]|yersiniabactin", re.I)),
           ("irp", re.compile(r"\birp[12]\b|HMWP", re.I)),
           ("fyuA", re.compile(r"fyuA|psn|pesticin.siderophore", re.I)),
           ("pga", re.compile(r"pga[ABCD]|pigment|Hms", re.I))]
    grows = []
    for line in open(CO92_ANNOT):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS" or not f[0].startswith("NC_003143"):
            continue
        s, e = int(f[3]), int(f[4])
        if e < 2100000 or s > 2300000:
            continue
        attrs = f[8]
        gn = ""
        m = re.search(r"(?:gene|Name)=([^;]+)", attrs)
        if m:
            gn = m.group(1)
        for tag, pat in grp:
            if pat.search(gn) or pat.search(attrs):
                grows.append([gn or "(%s)" % tag, s, e, tag])
                break
    w("fig3c_pgm_genes.csv", ["gene", "start", "end", "group"], grows, [CO92_ANNOT])


# ============================================================ Fig.4
def fig4a():
    w("fig4a_mechanism.csv", ["site_class", "handling", "affects_distance"], [
        ["jointly called, matches reference", "enters C", "no"],
        ["jointly called, differs", "enters diff", "yes"],
        ["not called in query (deletion)", "enters neither", "no"],
        ["not called in query (coverage gap)", "enters neither", "no"]], [])


def fig4b():
    rows = []
    for acc, lab in (("GCF_000009065.1", "CO92 (reference)"),
                     ("GCF_001293415.1", "CO92 (LANL 2015, pgm-/pPCP1-)")):
        for rep, st, en, ln in reads_fasta_run("results/snippy/%s/snps.aligned.fa" % acc):
            if ln >= 10000 and rep.startswith("NC_003143"):
                rows.append([lab, acc, rep, st, en, ln])
    w("fig4b_co92_tracks.csv",
      ["label", "accession", "replicon", "start", "end", "length_bp"], rows,
      ["results/snippy/GCF_000009065.1/snps.aligned.fa",
       "results/snippy/GCF_001293415.1/snps.aligned.fa"])


def fig4c():
    rows = []
    for mode, label in (("off", "off"), ("full", "full (naive two-state)"),
                        ("threestate", "threestate"), ("conservative", "conservative (used)")):
        k = sum(1 for _ in open("results/dedup/tierA_keep.%s.tsv" % mode)) - 1
        cl = rd("results/dedup/tierA_clusters.%s.tsv" % mode)
        multi = [r for r in cl if int(r["size"]) > 1]
        mixed = sum(1 for r in multi
                    if len({(next((x for x in rd("results/qc/deletion_profile_all.tsv")
                                   if x["accession"] == m), {}) or {}).get("replicons_absent", "")
                            for m in r["members"].split(",")}) > 1)
        rows.append([label, k, len(cl), len(multi)])
    w("fig4c_constraints.csv",
      ["constraint", "retained", "clusters", "multi_member_clusters"], rows,
      ["results/dedup/tierA_keep.{off,full,threestate,conservative}.tsv"])


def fig4d():
    src = os.path.expanduser("~/work/SY/logs/sweep_256.txt")
    if not os.path.exists(src):
        src = "logs/sweep_256.txt"
    if not os.path.exists(src):
        raise RuntimeError("sweep file not found")
    rows = []
    for line in open(src):
        f = line.split()
        if len(f) >= 7 and f[0].isdigit():
            rows.append([int(f[0]), int(f[2]), int(f[4]), int(f[5]),
                         f[6].split("/")[0]])
    w("fig4d_sweep.csv",
      ["threshold", "retained", "max_diameter", "max_size", "dup_pairs_merged"],
      rows, [src])


# ============================================================ Fig.5
def _merge_regions(blocks):
    """18 raw blocks -> 13 regions, merging overlaps only when polarity agrees."""
    out = []
    for rep in sorted({b["replicon"] for b in blocks}):
        bl = sorted([b for b in blocks if b["replicon"] == rep],
                    key=lambda x: int(x["start_bin"]))
        cur = None
        for b in bl:
            s, ln = int(b["start_bin"]), int(b["length_bin"])
            if cur and b["polarity"] == cur["polarity"] and s <= cur["end"] + 1:
                cur["end"] = max(cur["end"], s + ln)
                cur["members"].append(b)
            else:
                if cur:
                    out.append(cur)
                cur = {"replicon": rep, "start": s, "end": s + ln,
                       "polarity": b["polarity"], "members": [b]}
        if cur:
            out.append(cur)
    return out


def fig5a():
    blocks = rd("results/qc/outgroup_concordance.tsv")
    pol = {norm(r["replicon"]): r for r in rd("results/qc/recurrent_block_polarity.tsv")}
    # n_genomes from the binned recurrence
    cnt = Counter()
    for r in rd("results/qc/uncalled_blocks_all.tsv"):
        if int(r["length"]) >= 10000:
            cnt[(r["replicon"], int(r["start"]), int(r["length"]))] += 1
    regs = _merge_regions(blocks)
    rows = []
    for i, g in enumerate(sorted(regs, key=lambda x: (x["replicon"], x["start"])), 1):
        og1 = float(g["members"][0]["og1"]); og2 = float(g["members"][0]["og2"])
        n = max(cnt.get((m["replicon"], int(m["start_bin"]), int(m["length_bin"])), 0)
                for m in g["members"])
        hl = "primary" if (g["start"] <= 2365000 <= g["end"] or g["start"] <= 2370000 <= g["end"]) else (
             "secondary" if g["start"] <= 795000 <= g["end"] else "")
        rows.append(["region_%d" % i, g["replicon"], g["start"], g["end"],
                     g["end"] - g["start"], n, "%.3f" % og1, "%.3f" % og2,
                     g["polarity"], pol.get(norm(g["replicon"]), {}).get("n_loss_votes", ""),
                     hl])
    w("fig5a_region_polarity.csv",
      ["region_id", "replicon", "start", "end", "length_bp", "n_genomes",
       "og1_callable", "og2_callable", "polarity", "n_loss_votes", "highlight"],
      rows, ["results/qc/outgroup_concordance.tsv", "results/qc/uncalled_blocks_all.tsv"])


def fig5b():
    MISSING = set("Nn-.*")  # noqa

    def read_fa(p):
        seqs, name, buf = {}, None, []
        for line in open(p):
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                rep = line[1:].split()[0]
                name = rep.rsplit(".", 1)[0] if rep.rsplit(".", 1)[-1].isdigit() else rep
                buf = []
            else:
                buf.append(line.strip())
        if name is not None:
            seqs[name] = "".join(buf)
        return seqs

    def cf(seq, s, ln):
        win = seq[max(s - 1, 0):max(s - 1, 0) + ln]
        return sum(1 for c in win if c not in MISSING) / len(win) if win else None

    s1 = read_fa("results/snippy/GCA_000047365.1/snps.aligned.fa")
    s2 = read_fa("results/snippy/GCF_000016945.1/snps.aligned.fa")
    shared = sorted(set(s1) & set(s2))
    main = max(shared, key=lambda r: len(s1[r]))
    L = min(len(s1[main]), len(s2[main]))
    rng = random.Random(20261011)
    deltas = []
    for _ in range(1000):
        st = rng.randrange(1, max(L - 25000, 2))
        a, b = cf(s1[main], st, 25000), cf(s2[main], st, 25000)
        if a is not None and b is not None:
            deltas.append(abs(a - b))
    if not deltas:
        raise RuntimeError("no background windows")
    w("fig5b_background.csv", ["window_id", "delta"],
      [[i, "%.4f" % d] for i, d in enumerate(deltas, 1)], ["outgroups snps.aligned.fa"])
    deltas.sort()

    def pctl(d):
        return sum(1 for x in deltas if x <= d) / len(deltas) * 100
    rows = []
    for r in rd("results/qc/outgroup_concordance.tsv"):
        d = abs(float(r["og1"]) - float(r["og2"]))
        bpct = float(r["background_percentile"])
        rows.append(["%s:%s:%s" % (r["replicon"], r["start_bin"], r["length_bin"]),
                     "%.4f" % d, "%.1f" % bpct, int(bpct < 75)])
    w("fig5b_block_delta.csv",
      ["region_id", "delta", "percentile", "within_background"], rows,
      ["results/qc/outgroup_concordance.tsv"])



# ============================================================ Fig.1 (new schematic)
M1_NODES = """node_id;xmin;xmax;ymin;ymax;label;sublabel;style
assembly;20;80;99;107;assembly (one of 3,441);;plain
vcf;0;44;74;87;snps.vcf;differences only;plain
aligned;56;100;74;87;snps.aligned.fa;reference length, three states;blue
pass1;0;44;52;65;pass 1;union of variant sites;plain
pass2;56;100;52;65;pass 2: read union sites;ACGT -> 1-4, N / - / masked -> 0;blue
matrix;14;86;28;41;genotype matrix (int8);Tier A 256 x 12926\\nTier B 3441 x 67711;plain
raw;0;46;2;17;raw distance;over jointly called sites\\nthreshold applied here;orange
scaled;54;100;2;17;scaled distance;raw x n_sites / C\\ndiagnostic only;plain
"""
M1_EDGES = """edge_id;seq;x;y;arrow_end;dashed
stem;1;50;99;0;0
stem;2;50;95.5;0;0
split_l;1;50;95.5;1;0
split_l;2;50;91;1;0
split_l;3;22;91;1;0
split_l;4;22;87;1;0
split_r;1;50;95.5;1;0
split_r;2;50;91;1;0
split_r;3;78;91;1;0
split_r;4;78;87;1;0
v1;1;22;74;1;0
v1;2;22;65;1;0
v2;1;78;74;1;0
v2;2;78;65;1;0
dash;1;45;58.5;1;1
dash;2;55;58.5;1;1
conv_l;1;22;52;1;0
conv_l;2;22;46;1;0
conv_l;3;50;46;1;0
conv_l;4;50;41;1;0
conv_r;1;78;52;0;0
conv_r;2;78;46;0;0
conv_r;3;50;46;0;0
out_l;1;28;28;1;0
out_l;2;28;17;1;0
out_r;1;72;28;1;0
out_r;2;72;17;1;0
"""
M1_ELAB = """x;y;label;hjust
50;93;snippy 4.6.0 --ctgs, reference GCF_000009065.1;0.5
"""


def fig_m1():
    for name, txt in (("fig_m1_nodes.csv", M1_NODES),
                      ("fig_m1_edges.csv", M1_EDGES),
                      ("fig_m1_edge_labels.csv", M1_ELAB)):
        path = os.path.join(OUT, name); tmp = path + ".tmp"
        open(tmp, "w").write(txt)
        os.replace(tmp, path)
        open(path + ".prov", "w").write(
            "schematics geometry, fixed; written %s UTC\n" % NOW)
        PROV.append((name, txt.count("\n") - 1))
        print("  wrote %-30s %d rows" % (name, txt.count("\n") - 1))
    rows = []
    st_a = ["called_match", "called_match", "called_diff", "called_match", "called_match",
            "called_match", "called_match", "called_match", "called_match", "called_match",
            "called_match", "called_match", "called_diff", "called_match", "called_match", "called_match"]
    st_b = ["called_match", "called_match", "called_diff", "called_match", "called_match",
            "called_match", "uncalled", "uncalled", "uncalled", "uncalled", "uncalled",
            "called_match", "called_diff", "called_match", "called_match", "called_match"]
    for i, st in enumerate(st_a, 1):
        rows.append([i, "A", st])
    for i, st in enumerate(st_b, 1):
        rows.append([i, "B", st])
    w("fig_m1_sites.csv", ["site_index", "genome", "state"], rows,
      ["schematic; C/d counted in R"])


def main():
    print("preparing figure data -> %s/  (git %s, %s)" % (OUT, GIT, NOW))
    for fn in (fig1a, fig1b, fig1c, fig1d, fig2a, fig2b, fig2c, fig2d,
               fig3a, fig3b, fig3c, fig4a, fig4b, fig4c, fig4d, fig5a, fig5b, fig_m1):
        fn()
    print("\n%d CSVs written" % len(PROV))


if __name__ == "__main__":
    main()
