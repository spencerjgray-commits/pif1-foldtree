#!/usr/bin/env python3
"""
46_ale_anchor_posterior.py  --  read an ALEml_undated .uml_rec file and report where the ScPif1/ScRrm3
duplication lands in the sampled reconciled gene trees.

For every sampled reconciliation, the anchor event is the most recent common ancestor of the two
S. cerevisiae anchors (P07271 = Pif1, P38766 = Rrm3). ALE labels it ').D@<e>...' (duplication on species
branch e) or ').<e>...' (speciation). Species branches are numbered by ALE; they are mapped back to
the node names of the named species tree (workflow/44, v2fix: node_0 = Fungi, node_448 = crown
Saccharomycotina, sac_1 = Saccharomycotina minus Lipomyces) through leaf sets of S_with_ranks.

Outputs (stdout and --out TSV):
  ML rates and log-likelihood, number of samples, posterior frequency of the anchor event per
  species node, and the per-branch expected duplication counts for a few named nodes.

Usage (repo root, pif1 env):
  python workflow/46_ale_anchor_posterior.py --uml-rec X.ale.uml_rec --species-tree NAMED.nwk --out X.anchor.tsv
"""
import argparse
import csv
import gzip
import re
import sys
from collections import Counter

import dendropy

ANCHORS = ("P07271", "P38766")
SACCHAROMYCOTINA = "node_448"


def parse_newick(s):
    """Minimal newick parser: returns nested dicts {name, children}. Labels may contain . @ | > -"""
    i, n = 0, len(s)
    stack, root, cur = [], None, None
    node = {"name": "", "children": []}
    root = node
    while i < n:
        c = s[i]
        if c == "(":
            child = {"name": "", "children": []}
            node["children"].append(child)
            stack.append(node)
            node = child
            i += 1
        elif c == ",":
            parent = stack[-1]
            child = {"name": "", "children": []}
            parent["children"].append(child)
            node = child
            i += 1
        elif c == ")":
            node = stack.pop()
            i += 1
        elif c == ";":
            break
        else:
            j = i
            while j < n and s[j] not in ",();":
                j += 1
            tok = s[i:j]
            node["name"] = tok.split(":")[0]
            i = j
    # ALE writes the root as an extra wrapper "(...)...;": unwrap single-child roots
    while len(root["children"]) == 1 and not root["name"]:
        root = root["children"][0]
    return root


def leaves(nd):
    out, st = [], [nd]
    while st:
        x = st.pop()
        if x["children"]:
            st.extend(x["children"])
        else:
            out.append(x["name"])
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--uml-rec", required=True)
    ap.add_argument("--species-tree", required=True, help="the named binary species tree the ALE tree was made from")
    ap.add_argument("--out", required=True)
    ap.add_argument("--label", default="")
    args = ap.parse_args()

    op = gzip.open if args.uml_rec.endswith(".gz") else open
    lines = op(args.uml_rec, "rt").read().split("\n")
    S = next(l.split("\t", 1)[1] for l in lines if l.startswith("S:\t"))
    ll = float(next(l.split(":")[1] for l in lines if l.startswith(">logl:")))
    ml = next(l for l in lines if l.startswith("ML \t")).split("\t")
    dup, tra, loss = float(ml[1]), float(ml[2]), float(ml[3])
    i0 = next(i for i, l in enumerate(lines) if "reconciled G-s" in l)
    samples = []
    for l in lines[i0 + 1:]:
        if l.startswith("# of"):
            break
        if l.strip().startswith("("):
            samples.append(l.strip())
    branch_rows = {}
    for l in lines:
        p = l.split("\t")
        if p[0] in ("S_internal_branch", "S_terminal_branch") and len(p) > 3:
            branch_rows[p[1].split("(")[0]] = p

    # ALE branch ids -> named node by leaf set
    st = dendropy.Tree.get(path=args.species_tree, schema="newick", preserve_underscores=True, rooting="force-rooted")
    by_set = {}
    for nd in st.postorder_node_iter():
        ls = frozenset(l.taxon.label for l in nd.leaf_iter())
        by_set[ls] = nd.label if not nd.is_leaf() else nd.taxon.label
    sac = frozenset(l.taxon.label for l in st.find_node_with_label(SACCHAROMYCOTINA).leaf_iter())
    aleS = parse_newick(S)
    id2name, id2size = {}, {}
    st_ = [aleS]
    while st_:
        x = st_.pop()
        ls = frozenset(leaves(x))
        nm = x["name"]
        if x["children"]:
            id2name[nm] = by_set.get(ls, f"?{nm}")
            id2size[nm] = len(ls)
            st_.extend(x["children"])
        else:
            id2name[nm] = nm
            id2size[nm] = 1

    events, sizes, insac = Counter(), Counter(), 0
    for s in samples:
        t = parse_newick(s)
        # smallest clade holding both anchors
        best = None
        st_ = [t]
        while st_:
            x = st_.pop()
            if not x["children"]:
                continue
            lv = leaves(x)
            acc = {v.split("_")[1].split(".")[0] for v in lv}
            if all(a in acc for a in ANCHORS):
                if best is None or len(lv) < best[1]:
                    best = (x, len(lv))
                st_.extend(x["children"])
        lab = best[0]["name"]
        m = re.match(r"\.(D@|T@)?([^.|>-]+)", lab)
        ev = {"D@": "D", "T@": "T", None: "S"}[m.group(1)]
        e = m.group(2)
        nm = id2name.get(e, f"?{e}")
        events[(ev, nm, id2size.get(e, 0))] += 1
        sizes[best[1]] += 1
        sp_leaves = {v.split("_")[0] for v in leaves(best[0])}
        insac += sp_leaves <= sac

    n = len(samples)
    print(f"{args.label} {args.uml_rec}\n  ML rates D={dup:.4g} T={tra:.3g} L={loss:.4g}; logL={ll:.3f}; samples={n}")
    print(f"  anchor-clade gene content all Saccharomycotina: {insac}/{n}")
    print(f"  anchor clade sizes: {sizes.most_common(6)}")
    rows = []
    for (ev, nm, k), c in events.most_common():
        print(f"  anchor event {ev} @ {nm} ({k} species leaves): {c}/{n} = {c / n:.3f}")
        rows.append({"label": args.label, "event": ev, "species_node": nm, "species_leaves": k, "count": c, "samples": n,
                     "freq": round(c / n, 4), "D": dup, "L": loss, "T": tra, "logL": ll, "anchor_clade_all_sacch": insac})
    name2id = {v: k for k, v in id2name.items()}
    for nm in ("node_0", "node_1", "node_3", "node_448", "sac_1"):
        b = name2id.get(nm)
        if b in branch_rows:
            p = branch_rows[b]
            print(f"  expected per-branch counts at {nm} (ALE id {b}): D={float(p[2]):.3f} L={float(p[4]):.3f}")
    with open(args.out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()
