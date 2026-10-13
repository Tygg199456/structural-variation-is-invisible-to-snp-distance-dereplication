#!/usr/bin/env python3
"""Replicon three-state inventory for SY1011 (external review 2026-10-11).

Revised: the plasmid identity is read from the plasmid's own DESCRIPTION (the
text after the word "plasmid"), NOT from the whole header.  Two failure modes
were found and fixed:

  * a strain NAME may contain a plasmid label -- "strain 6/69_pPCP1- plasmid
    pMT1" made the pMT1 sequence match the pPCP1 pattern (the strain is a
    pPCP1-cured derivative whose name records the phenotype);
  * an unnamed plasmid could be assigned by nearest length -- a 38 kb unnamed
    plasmid was matched to the 9.6 kb pPCP1.

Length fallback is now only used for an unnamed plasmid and only when within
25% of a reference plasmid's length; otherwise the sequence is recorded as
"other" and does not create a CARRIED call.

The submitted sequence list is cached so re-classification does not re-read
17 GB of FASTA.
"""
import os, glob, csv, re, collections

ROOT = os.environ.get("PLAGUE_DECAY_ROOT",
                       os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.chdir(ROOT)
GN = "data/genomes"
CACHE = "results/qc/_replinv_seqs.tsv"
OUT = "results/qc/replicon_three_state.tsv"
REFP = {"pCD1": 70305, "pPCP1": 9612, "pMT1": 96210}
NAMEPAT = {"pCD1": re.compile(r"pCD1", re.I),
           "pPCP1": re.compile(r"pPCP1", re.I),
           "pMT1": re.compile(r"pMT1", re.I)}
PLASMID = re.compile(r"plasmid", re.I)
DESC = re.compile(r"plasmid\s+([^,;]+)", re.I)
TAG2P = {"NC_003131": "pCD1", "NC_003132": "pPCP1", "NC_003134": "pMT1"}
TOL = 0.25
norm = lambda a: re.sub(r"\.\d+$", "", (a or "").strip())


def scan(path):
    out = []
    hdr = None
    n = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if hdr is not None:
                    out.append((hdr, n))
                hdr = line[1:].rstrip("\n")
                n = 0
            else:
                n += len(line.strip())
    if hdr is not None:
        out.append((hdr, n))
    return out


def load_cache():
    cache = collections.defaultdict(list)
    if os.path.exists(CACHE):
        for r in csv.DictReader(open(CACHE), delimiter="\t"):
            cache[r["accession"]].append((r["header"], int(r["length"])))
    return cache


def get_seqs(d, cache, write):
    if d in cache and cache[d]:
        return cache[d]
    fnas = sorted(glob.glob(os.path.join(GN, d, "*.fna")))
    if not fnas:
        return []
    ss = scan(fnas[0])
    for h, l in ss:
        write.writerow([d, h.replace("\t", " "), l])
    cache[d] = ss
    return ss


def classify(plasmid):
    state, ev = {}, {}
    if not plasmid:
        for p in REFP:
            state[p] = "UNRESOLVED"
            ev[p] = "no plasmid sequence submitted"
        return state, ev
    assigned = collections.defaultdict(list)
    for h, l in plasmid:
        m = DESC.search(h)
        desc = m.group(1) if m else h
        hit = None
        for p, pat in NAMEPAT.items():
            if pat.search(desc):
                hit = p
                break
        if hit is None:
            cands = [(abs(l - REFP[p]), p) for p in REFP if abs(l - REFP[p]) <= TOL * REFP[p]]
            if cands:
                hit = min(cands)[1]
        if hit:
            assigned[hit].append((h, l))
    for p in REFP:
        if assigned.get(p):
            state[p] = "CARRIED"
            ev[p] = ";".join("%s(%d)" % (h, l) for h, l in assigned[p])
        else:
            state[p] = "LOST"
            ev[p] = "plasmids submitted, not " + p
    return state, ev


old = {norm(r["accession"]): (r.get("replicons_absent") or "").strip()
       for r in csv.DictReader(open("results/qc/deletion_profile_all.tsv"), delimiter="\t")}

cache = load_cache()
mode = "a" if cache else "w"
with open(CACHE, mode, newline="") as cachefh:
    cw = csv.writer(cachefh, delimiter="\t")
    if mode == "w":
        cw.writerow(["accession", "header", "length"])
    rows = []
    for d in sorted(os.listdir(GN)):
        dd = os.path.join(GN, d)
        if not os.path.isdir(dd):
            continue
        ss = get_seqs(d, cache, cw)
        if not ss:
            continue
        plasmid = [(h, l) for h, l in ss if PLASMID.search(h)]
        state, ev = classify(plasmid)
        rows.append((d, len(ss), len(plasmid), sum(l for _, l in ss), state, ev, plasmid))

with open(OUT, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["accession", "n_seq", "n_plasmid_seq", "total_bp",
                "pCD1_state", "pPCP1_state", "pMT1_state",
                "old_replicons_absent", "evidence"])
    for d, n_seq, npl, tot, state, ev, plasmid in rows:
        allev = "; ".join("%s|%d" % (h, l) for h, l in plasmid) or "(none)"
        w.writerow([d, n_seq, npl, tot, state["pCD1"], state["pPCP1"], state["pMT1"],
                    old.get(norm(d), ""), allev])

c = collections.Counter()
for d, n_seq, npl, tot, state, ev, plasmid in rows:
    for p in REFP:
        c[(p, state[p])] += 1
print("genomes scanned:", len(rows))
for p in REFP:
    print("  %-6s CARRIED=%-4d LOST=%-4d UNRESOLVED=%-4d"
          % (p, c[(p, "CARRIED")], c[(p, "LOST")], c[(p, "UNRESOLVED")]))
print()
print("=== old claim vs new three-state (all-profile) ===")
cf = collections.Counter()
for d, n_seq, npl, tot, state, ev, plasmid in rows:
    if norm(d) not in old:
        continue
    oa = set(old.get(norm(d), "").split(";")) if old.get(norm(d), "") else set()
    for tag, p in TAG2P.items():
        cf[("OLD_absent" if tag in oa else "OLD_present", state[p])] += 1
for k in sorted(cf):
    print("  %-12s -> %-11s %5d" % (k[0], k[1], cf[k]))
print()
print("WROTE " + OUT)
