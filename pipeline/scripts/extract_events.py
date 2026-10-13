#!/usr/bin/env python3
"""
PastML work dir -> tidy gene/branch/event table.

DECISION THAT MUST NOT BE LOST:
Count a loss once per branch; repeated marginal states on one branch are
not independent events.

A branch is one edge of the reconstructed tree. PastML reports the marginal
ancestral state at each NODE; an event lives on the EDGE between a node and
its parent. Reading states off nodes and calling each one an event double
counts: one loss produces a "lost" state at every descendant node, and the
null model downstream would then see dozens of independent losses where there
was one.

So the transition is computed per edge, parent state -> child state, and at
most one event is emitted per (gene, edge).

Called as:  python scripts/extract_events.py <pastml_work_dir> > events.tsv
Columns: gene, parent, child, branch_length, event
Files in the work dir are discovered rather than assumed -- PastML's naming
changed across versions.
"""

import csv
import os
import re
import sys
from pathlib import Path

LOSS_TOKENS = {"0", "absent", "loss", "lost", "no"}
PRESENT_TOKENS = {"1", "present", "gain", "yes"}
PSEUDO_TOKENS = {"2", "pseudogene", "pseudogenisation", "pseudogenization", "pseudo"}


def find_first(root, patterns):
    for pat in patterns:
        hits = sorted(Path(root).rglob(pat))
        if hits:
            return hits[0]
    return None


def read_tree_edges(path):
    """Parse a Newick tree into {(parent, child): branch_length} and the node set.

    Uses dendropy when available (it is in envs/py.yaml), falling back to a
    minimal parser so the script is testable without the full env.
    """
    try:
        import dendropy
        t = dendropy.Tree.get(path=str(path), schema="newick",
                             preserve_underscores=True)
        edges = {}
        for e in t.edges():
            if e.tail_node is None:
                continue
            edges[(e.tail_node.label or "", e.head_node.label or "")] = (
                e.length if e.length is not None else 0.0)
        return edges
    except ImportError:
        return _minimal_newick_edges(path)


def _minimal_newick_edges(path):
    """Very small Newick reader: enough for internal-node-labelled trees."""
    text = Path(path).read_text().strip()

    edges = {}

    pos = [0]

    def parse_node():
        name = ""
        if text[pos[0]] == "(":
            pos[0] += 1
            children = []
            while True:
                children.append(parse_node())
                if text[pos[0]] == ",":
                    pos[0] += 1
                    continue
                if text[pos[0]] == ")":
                    pos[0] += 1
                    break
            # optional internal label
            start = pos[0]
            while pos[0] < len(text) and text[pos[0]] not in ",():":
                pos[0] += 1
            name = text[start:pos[0]]
            length = 0.0
            if pos[0] < len(text) and text[pos[0]] == ":":
                pos[0] += 1
                start = pos[0]
                while pos[0] < len(text) and text[pos[0]] not in ",()":
                    pos[0] += 1
                try:
                    length = float(text[start:pos[0]])
                except ValueError:
                    length = 0.0
            for child_name, _cl in children:
                edges[(name, child_name)] = 0.0
            for cname, _l in children:
                pass
            return name, length
        else:
            start = pos[0]
            while pos[0] < len(text) and text[pos[0]] not in ",():;":
                pos[0] += 1
            name = text[start:pos[0]]
            length = 0.0
            if pos[0] < len(text) and text[pos[0]] == ":":
                pos[0] += 1
                start = pos[0]
                while pos[0] < len(text) and text[pos[0]] not in ",()":
                    pos[0] += 1
                try:
                    length = float(text[start:pos[0]])
                except ValueError:
                    length = 0.0
            return name, length

    parse_node()
    return edges


def read_states_table(path):
    """PastML states tab -> {character: {node: state_token}}.

    Values may be a bare state or 'state:probability'; only the state is kept.
    """
    table = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        nodes = [h.strip() for h in header[1:]]
        for row in reader:
            if not row:
                continue
            gene = row[0].strip()
            if not gene:
                continue
            states = {}
            for node, raw in zip(nodes, row[1:]):
                tok = raw.strip().split(":")[0].strip()
                states[node] = tok
            table[gene] = states
    return table


def classify(token):
    t = str(token).strip().lower()
    if t in PSEUDO_TOKENS:
        return "pseudogene"
    if t in LOSS_TOKENS:
        return "absent"
    if t in PRESENT_TOKENS:
        return "present"
    try:
        v = float(t)
    except ValueError:
        return "unknown"
    if v == 0:
        return "absent"
    if v == 1:
        return "present"
    if v == 2:
        return "pseudogene"
    return "unknown"


SET_SPLIT = re.compile(r"[|,;/]+")


def split_state_set(raw):
    """PastML MPPA can return a SET of states for a node (Brier-optimal set),
    not a single state. Return the list of state tokens.

        '1'        -> ['1']
        '1:0.98'   -> ['1']      (probability suffix dropped)
        '1|0'      -> ['1','0']  (ambiguous)
        '{1,0}'    -> ['1','0']
    """
    tok = str(raw).strip()
    if not tok:
        return []
    tok = tok.strip("{}")
    parts = []
    for chunk in SET_SPLIT.split(tok):
        chunk = chunk.split(":")[0].strip()
        if chunk:
            parts.append(chunk)
    return parts


def state_class(raw):
    """-> (class, is_ambiguous). A non-singleton set means the ASR did not
    resolve that node, and any transition through it must not be counted as
    an independent event."""
    parts = split_state_set(raw)
    if not parts:
        return "unknown", True
    classes = {classify(p) for p in parts}
    if len(parts) > 1:
        return (classes.pop() if len(classes) == 1 else "unknown"), True
    return classify(parts[0]), False


def main():
    if len(sys.argv) < 2:
        raise SystemExit("usage: extract_events.py <pastml_work_dir>")
    work = Path(sys.argv[1])
    if not work.is_dir():
        raise SystemExit(f"FATAL: {work} is not a directory")

    tree_file = find_first(work, ["*.nwk", "*.newick", "*.tre", "*.tree"])
    if tree_file is None:
        raise SystemExit(f"FATAL: no tree file found under {work}")

    states_file = find_first(work, [
        "*combined_ancestral_states.tab", "*ancestral_states.tab",
        "*_states.tab", "*marginal_states.tab", "*.tab",
    ])
    if states_file is None:
        raise SystemExit(f"FATAL: no ancestral states table found under {work}")

    edges = read_tree_edges(tree_file)
    states = read_states_table(states_file)

    if not edges:
        raise SystemExit(f"FATAL: parsed no edges from {tree_file}")
    if not states:
        raise SystemExit(f"FATAL: parsed no characters from {states_file}")

    out = csv.writer(sys.stdout, delimiter="\t")
    out.writerow(["gene", "parent", "child", "branch_length", "event",
                  "confidence"])

    n_events = {"loss": 0, "gain": 0, "pseudogenisation": 0}
    n_ambiguous = 0
    n_transitions = 0
    n_geneless = 0
    for gene, node_states in states.items():
        any_edge = False
        for (parent, child), length in edges.items():
            ps = node_states.get(parent)
            cs = node_states.get(child)
            if ps is None or cs is None:
                continue
            any_edge = True
            n_transitions += 1
            p_class, p_amb = state_class(ps)
            c_class, c_amb = state_class(cs)
            # An unresolved endpoint makes the transition uninterpretable.
            # Counting it would inflate the per-gene loss count, and the null
            # model estimates each gene's rate FROM that count -- numerator
            # and rate move together and the test loses power.
            if p_amb or c_amb:
                n_ambiguous += 1
                out.writerow([gene, parent, child, f"{length:.6g}",
                              "ambiguous", "ambiguous"])
                continue
            event = None
            if p_class == "present" and c_class == "absent":
                event = "loss"
            elif p_class == "present" and c_class == "pseudogene":
                event = "pseudogenisation"
            elif p_class in ("absent", "pseudogene") and c_class == "present":
                event = "gain"
            if event:
                n_events[event if event in n_events else "loss"] += 1
                out.writerow([gene, parent, child, f"{length:.6g}", event,
                              "unambiguous"])
        if not any_edge:
            n_geneless += 1

    total = sum(n_events.values())
    amb_frac = n_ambiguous / max(n_transitions, 1)
    sys.stderr.write(
        f"[extract_events] genes={len(states)} edges={len(edges)} "
        f"transitions={n_transitions} events={total} "
        f"(loss={n_events['loss']} gain={n_events['gain']} "
        f"pseudo={n_events['pseudogenisation']}) "
        f"ambiguous={n_ambiguous} ({amb_frac:.1%})\n"
    )
    if amb_frac > 0.20:
        sys.stderr.write(
            "[extract_events] WARNING: >20% of transitions are ambiguous. "
            "The ancestral reconstruction is too unresolved to support the "
            "parallelism test; revisit the PastML model or reduce taxa "
            "rather than proceeding.\n"
        )
    if total == 0:
        raise SystemExit(
            "FATAL: zero state transitions. Node labels in the states table "
            "and in the tree did not join -- check that both use the same "
            "internal node naming."
        )


if __name__ == "__main__":
    main()
