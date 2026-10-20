#!/usr/bin/env python3
"""
33_split_supports.py  --  read SH-aLRT / UFBoot off a tree BY BIPARTITION, never by node after rerooting.

Why: IQ-TREE writes branch supports as internal-node labels ("98.2/29"). A label belongs to the
branch above that node *in the tree as written*. Rerooting with ete3/dendropy while the label stays
on the node object moves a support onto a different branch (2026-10-15 review, MET-M6 / RES-R7).
Here every label is converted to an unrooted split (the tip set on one side of the branch) the
moment the tree is parsed, so a "clade" in any rooting is looked up as a split and gets the
support of that exact branch.

Clades are defined by rooting on a single outgroup tip (default hPIF1, Q9H611): a clade is the
side of a split that does not contain the root tip. This is equivalent to outgroup rooting but
never moves a label.

Usage (from repo root):
    python workflow/33_split_supports.py --tree results/seq_tree/pif1_aa3di_supported.treefile
    python workflow/33_split_supports.py --tree results/lba/lba_noEDF.treefile --root A0A074S981
    python workflow/33_split_supports.py --tree T --write-clade out.txt   # tip list of the anchor MRCA
    python workflow/33_split_supports.py --tree T --query-clade tips.txt  # support of a given tip set
"""
import argparse
import csv
import json
import re
import sys
from collections import Counter


# ---------------------------------------------------------------- newick
def parse_newick(s):
    """Minimal newick parser -> nested dicts {name, label, children}. Labels kept verbatim."""
    s = s.strip().rstrip(";")
    pos = 0

    def read_name():
        nonlocal pos
        if pos < len(s) and s[pos] == "'":
            end = s.index("'", pos + 1)
            name = s[pos + 1:end]
            pos = end + 1
            return name
        m = re.compile(r"[^,():;\[]*").match(s, pos)
        pos = m.end()
        return m.group(0)

    def skip_branch():
        nonlocal pos
        if pos < len(s) and s[pos] == "[":           # NHX / comment
            pos = s.index("]", pos) + 1
        if pos < len(s) and s[pos] == ":":
            m = re.compile(r":[^,();\[]*").match(s, pos)
            pos = m.end()
        if pos < len(s) and s[pos] == "[":
            pos = s.index("]", pos) + 1

    def node():
        nonlocal pos
        if s[pos] == "(":
            pos += 1
            kids = [node()]
            while s[pos] == ",":
                pos += 1
                kids.append(node())
            assert s[pos] == ")", f"parse error at {pos}"
            pos += 1
            label = read_name()
            skip_branch()
            return {"name": None, "label": label, "children": kids}
        name = read_name()
        skip_branch()
        return {"name": name, "label": None, "children": []}

    return node()


def parse_label(label):
    """'98.2/29' -> (98.2, 29); '95' -> (None, 95.0) is ambiguous, so return raw too."""
    if not label:
        return None, None
    parts = label.split("/")
    try:
        vals = [float(p) for p in parts]
    except ValueError:
        return None, None
    if len(vals) >= 2:
        return vals[0], vals[1]
    return None, vals[0]


def splits_of(tree, all_tips=None):
    """Return {frozenset(side without ref tip): label} for every internal branch, plus the tip list."""
    out = {}
    tips = []

    def walk(n):
        if not n["children"]:
            tips.append(n["name"])
            return frozenset([n["name"]])
        below = frozenset().union(*(walk(c) for c in n["children"]))
        n["_below"] = below
        return below

    walk(tree)
    universe = frozenset(tips)

    def collect(n, is_root):
        if not n["children"]:
            return
        if not is_root:
            out[n["_below"]] = n["label"]
        for c in n["children"]:
            collect(c, False)

    collect(tree, True)
    return out, universe


class SplitTree:
    def __init__(self, path, root_tip):
        txt = open(path).read().strip()
        self.tree = parse_newick(txt.split(";")[0])
        raw, self.universe = splits_of(self.tree)
        self.root_tip = self._resolve(root_tip)
        # orient every split as the side NOT containing the root tip -> a rooted clade
        self.clades = {}
        for side, lab in raw.items():
            clade = side if self.root_tip not in side else self.universe - side
            if 1 < len(clade) < len(self.universe) - 1:
                # the two root-adjacent branches of an unrooted (basal-trifurcation) tree are one split;
                # keep the first non-empty label
                if clade not in self.clades or not self.clades[clade]:
                    self.clades[clade] = lab

    def _resolve(self, acc_or_tip):
        if acc_or_tip in self.universe:
            return acc_or_tip
        hits = [t for t in self.universe if t.split("_")[0] == acc_or_tip]
        if len(hits) != 1:
            sys.exit(f"cannot resolve tip '{acc_or_tip}' ({len(hits)} hits)")
        return hits[0]

    def support(self, tipset):
        tipset = frozenset(tipset)
        clade = tipset if self.root_tip not in tipset else self.universe - tipset
        if clade in self.clades:
            return True, self.clades[clade]
        return False, None

    def mrca(self, tips):
        tips = {self._resolve(t) for t in tips}
        best = None
        for c in self.clades:
            if tips <= c and (best is None or len(c) < len(best)):
                best = c
        if best is None and self.root_tip not in tips:
            # only the root tip lies outside: the MRCA is the whole ingroup (its "support" is the
            # root tip's pendant edge, which has no label)
            best = self.universe - {self.root_tip}
        return best

    def children(self, clade):
        """Children of a rooted clade = maximal proper sub-clades (+ singleton tips not covered)."""
        subs = [c for c in self.clades if c < clade]
        maximal = [c for c in subs if not any(c < d for d in subs)]
        covered = frozenset().union(*maximal) if maximal else frozenset()
        maximal += [frozenset([t]) for t in clade - covered]
        return sorted(maximal, key=len, reverse=True)


# ---------------------------------------------------------------- taxonomy
def load_tax(tip_map, manifest):
    acc_of = {}
    for r in csv.DictReader(open(tip_map), delimiter="\t"):
        acc_of[r["tip_label"]] = r["accession"]
    tax = {}
    for r in csv.DictReader(open(manifest)):
        tax[r["accession"]] = r
    return acc_of, tax


def compose(tips, acc_of, tax, field):
    return Counter((tax.get(acc_of.get(t, t.split("_")[0]), {}).get(field) or "(none)") for t in tips)


def fmt(label):
    a, b = parse_label(label)
    if a is None and b is None:
        return "NA"
    return f"{a:g}/{b:g}" if a is not None else f"{b:g}"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--root", default="Q9H611", help="outgroup tip (accession or full tip label)")
    ap.add_argument("--anchors", default="P07271,P38766", help="ScPif1,ScRrm3 accessions")
    ap.add_argument("--tip-map", default="data/seqs/tip_map.tsv")
    ap.add_argument("--manifest", default="manifest.csv")
    ap.add_argument("--write-clade", help="write the anchor-MRCA tip list here")
    ap.add_argument("--query-clade", action="append", default=[],
                    help="file with one tip per line; report whether it is a split and its support")
    ap.add_argument("--json", help="write a JSON summary here")
    args = ap.parse_args()

    st = SplitTree(args.tree, args.root)
    acc_of, tax = load_tax(args.tip_map, args.manifest)
    a1, a2 = args.anchors.split(",")
    m = st.mrca([a1, a2])
    res = {"tree": args.tree, "root_tip": st.root_tip, "n_tips": len(st.universe)}
    print(f"tree: {args.tree}\nrooted on: {st.root_tip}   tips: {len(st.universe)}")
    if m is None:
        print("anchor MRCA: none below root (anchors straddle the root)")
    else:
        ok, lab = st.support(m)
        sub = compose(m, acc_of, tax, "subphylum")
        print(f"anchor MRCA: {len(m)} tips | {len(st.universe) - len(m)}  support (SH-aLRT/UFBoot) = {fmt(lab)}")
        print(f"  subphyla: {dict(sub.most_common())}")
        res["mrca"] = {"size": len(m), "support": fmt(lab), "subphyla": dict(sub)}
        kids = st.children(m)
        res["daughters"] = []
        for k in kids:
            _, kl = st.support(k)
            which = "Pif1" if st._resolve(a1) in k else ("Rrm3" if st._resolve(a2) in k else "-")
            fam = compose(k, acc_of, tax, "family")
            print(f"  daughter ({which}): {len(k)} tips  support = {fmt(kl) if len(k) > 1 else 'tip'}"
                  f"  families: {len([f for f in fam if f != '(none)'])} {dict(fam.most_common())}")
            res["daughters"].append({"paralog": which, "size": len(k), "support": fmt(kl),
                                     "families": dict(fam)})
        if args.write_clade:
            with open(args.write_clade, "w") as fh:
                fh.write("\n".join(sorted(m)) + "\n")
    res["queries"] = []
    for q in args.query_clade:
        tips = [l.strip() for l in open(q) if l.strip()]
        missing = [t for t in tips if t not in st.universe]
        tips = [t for t in tips if t in st.universe]
        ok, lab = st.support(tips)
        print(f"query {q}: {len(tips)} tips (missing {len(missing)}) -> "
              f"{'split present, support ' + fmt(lab) if ok else 'NOT a split in this tree'}")
        res["queries"].append({"file": q, "n": len(tips), "missing": len(missing),
                               "present": ok, "support": fmt(lab) if ok else None})
    if args.json:
        json.dump(res, open(args.json, "w"), indent=1)


if __name__ == "__main__":
    main()
