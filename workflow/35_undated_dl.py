#!/usr/bin/env python3
"""
35_undated_dl.py  --  maximum-likelihood UndatedDL reconciliation of one unrooted gene tree with a
rooted binary species tree, re-implementing what `generax --strategy EVAL -r UndatedDL
--unrooted-gene-tree` does with fixed global D/L rates.

Why this exists: the only GeneRax on this Mac (2.1.3, Intel build under Rosetta) segfaults in
inferReconciliation before writing event counts, and 2.1.3 also re-optimises D/L rates, whereas the
archived 2.0.4 runs used fixed D = L = 0.2. This script follows GeneRax 2.0.4's UndatedDLModel:

  per species node e:  PD = d/(1+d+l), PL = l/(1+d+l), PS = 1/(1+d+l)
  extinction uE[e]:    smaller root of PD x^2 - x + (PL + PS uE[f] uE[g])        (postorder)
  gene leaf u, leaf e: P = PS[e] if u maps to e else 0
  otherwise:           P(u,e) = [ PS(P(l,f)P(r,g) + P(l,g)P(r,f))      speciation  (e internal)
                                + PD P(l,e)P(r,e)                       duplication (u internal)
                                + PS(P(u,f)uE[g] + P(u,g)uE[f]) ]       speciation-loss
                                / (1 - 2 PD uE[e])                      duplication-loss
  root:   every gene-tree edge is tried as the root; the ML scenario starts from the
          (root edge, species node) pair with the largest P(virtual root, e) (uniform origination)
  scenario: max backtrace; at each (u,e) the event with the largest term wins (first max, in the
          order S(l->f) S(l->g) D SL(->f) SL(->g), as std::max_element does in GeneRax).

Vectors are stored per gene clade with a log scale so nothing underflows.

Validated against the archived GeneRax 2.0.4 runs before use (see RERUN_2026-10_SUMMARY.md).

Outputs in --outdir:
  speciesEventCounts.txt   <species node> <S+SL+leaf> <D> <SL> 0   (GeneRax 2.0.4 column layout)
  eventCounts.txt          S, SL, D, Leaf totals
  summary.json             totals, root, anchor-duplication mapping, top nodes by D

Usage (repo root):
  python workflow/35_undated_dl.py --gene-tree G.nwk --species-tree S.nwk \
      --map data/species_tree/gene_species.map --outdir OUT [--dup 0.2 --loss 0.2]
"""
import argparse
import json
import math
import os
import sys
from collections import Counter, defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from importlib import import_module

parse_newick = import_module("33_split_supports").parse_newick


# ---------------------------------------------------------------- species tree
class Species:
    def __init__(self, path):
        root = parse_newick(open(path).read().strip().split(";")[0])
        self.names, self.left, self.right, self.leaves = [], [], [], {}
        self.height = []

        def walk(n):  # postorder indexing
            kids = [walk(c) for c in n["children"]]
            if kids and len(kids) != 2:
                sys.exit(f"species tree must be binary (node with {len(kids)} children)")
            i = len(self.names)
            nm = n["name"] if not kids else (n["label"] or f"internal_{i}")
            self.names.append(nm)
            self.left.append(kids[0] if kids else -1)
            self.right.append(kids[1] if kids else -1)
            self.height.append(0 if not kids else 1 + max(self.height[kids[0]], self.height[kids[1]]))
            if not kids:
                self.leaves[nm] = i
            return i

        sys.setrecursionlimit(100000)
        self.root = walk(root)
        self.n = len(self.names)
        self.left = np.array(self.left)
        self.right = np.array(self.right)
        self.height = np.array(self.height)
        self.internal = np.where(self.left >= 0)[0]
        self.levels = [np.where(self.height == h)[0] for h in range(1, self.height.max() + 1)]
        # descendant leaf sets for labelling
        self.desc = [None] * self.n
        for i in range(self.n):
            if self.left[i] < 0:
                self.desc[i] = frozenset([self.names[i]])
            else:
                self.desc[i] = self.desc[self.left[i]] | self.desc[self.right[i]]


# ---------------------------------------------------------------- gene tree (unrooted)
class GeneTree:
    def __init__(self, path, keep=None):
        root = parse_newick(open(path).read().strip().split(";")[0])
        self.adj = defaultdict(list)
        self.leafname = {}
        cnt = [0]

        def walk(n, parent):
            i = cnt[0]
            cnt[0] += 1
            if not n["children"]:
                self.leafname[i] = n["name"]
            if parent is not None:
                self.adj[i].append(parent)
                self.adj[parent].append(i)
            for c in n["children"]:
                walk(c, i)

        sys.setrecursionlimit(100000)
        walk(root, None)
        if keep is not None:
            self._prune(set(keep))
        # suppress degree-2 nodes (e.g. a rooted input)
        changed = True
        while changed:
            changed = False
            for v in list(self.adj):
                if v not in self.leafname and len(self.adj[v]) == 2:
                    a, b = self.adj[v]
                    self.adj[a].remove(v); self.adj[b].remove(v)
                    self.adj[a].append(b); self.adj[b].append(a)
                    del self.adj[v]
                    changed = True
        bad = [v for v in self.adj if v not in self.leafname and len(self.adj[v]) != 3]
        if bad:
            sys.exit(f"gene tree must be binary-unrooted ({len(bad)} nodes of degree != 3)")
        self.edges = sorted({(min(a, b), max(a, b)) for a in self.adj for b in self.adj[a]})

    def _prune(self, keep):
        for v, nm in list(self.leafname.items()):
            if nm not in keep:
                for u in self.adj[v]:
                    self.adj[u].remove(v)
                del self.adj[v]
                del self.leafname[v]
        # remove internal nodes that became leaves
        changed = True
        while changed:
            changed = False
            for v in list(self.adj):
                if v not in self.leafname and len(self.adj[v]) <= 1:
                    for u in self.adj[v]:
                        self.adj[u].remove(v)
                    del self.adj[v]
                    changed = True


# ---------------------------------------------------------------- rate optimisation
def nelder_mead(f, x0, step=0.3, tol=1e-4, max_iter=200):
    """Minimal Nelder-Mead minimiser (no scipy in the pif1 env)."""
    n = len(x0)
    pts = [list(x0)] + [[x0[j] + (step if j == i else 0.0) for j in range(n)] for i in range(n)]
    vals = [f(p) for p in pts]
    for _ in range(max_iter):
        order = sorted(range(n + 1), key=lambda i: vals[i])
        pts, vals = [pts[i] for i in order], [vals[i] for i in order]
        if abs(vals[-1] - vals[0]) < tol:
            break
        cen = [sum(p[j] for p in pts[:-1]) / n for j in range(n)]
        xr = [cen[j] + (cen[j] - pts[-1][j]) for j in range(n)]; fr = f(xr)
        if fr < vals[0]:
            xe = [cen[j] + 2 * (cen[j] - pts[-1][j]) for j in range(n)]; fe = f(xe)
            pts[-1], vals[-1] = (xe, fe) if fe < fr else (xr, fr)
        elif fr < vals[-2]:
            pts[-1], vals[-1] = xr, fr
        else:
            xc = [cen[j] + 0.5 * (pts[-1][j] - cen[j]) for j in range(n)]; fc = f(xc)
            if fc < vals[-1]:
                pts[-1], vals[-1] = xc, fc
            else:
                pts = [pts[0]] + [[pts[0][j] + 0.5 * (p[j] - pts[0][j]) for j in range(n)] for p in pts[1:]]
                vals = [vals[0]] + [f(p) for p in pts[1:]]
    i = min(range(n + 1), key=lambda k: vals[k])
    return pts[i], vals[i]


def optimise_rates(sp, gt, gene2sp, d0=0.2, l0=0.2):
    """ML global D/L rates (as GeneRax optimises them before reconciling): maximise sum_loglik."""
    cache = {}

    def negll(x):
        key = (round(x[0], 6), round(x[1], 6))
        if key not in cache:
            cache[key] = -UndatedDL(sp, gt, gene2sp, math.exp(x[0]), math.exp(x[1])).sum_loglik()
        return cache[key]

    x, v = nelder_mead(negll, [math.log(d0), math.log(l0)])
    return math.exp(x[0]), math.exp(x[1]), -v, len(cache)


# ---------------------------------------------------------------- reconciliation
class UndatedDL:
    def __init__(self, sp, gt, gene2sp, dup, loss):
        self.sp, self.gt = sp, gt
        s = 1.0 + dup + loss
        self.PD, self.PL, self.PS = dup / s, loss / s, 1.0 / s
        uE = np.zeros(sp.n)
        for e in range(sp.n):  # postorder
            c = self.PL
            if sp.left[e] >= 0:
                c += self.PS * uE[sp.left[e]] * uE[sp.right[e]]
            a = self.PD
            uE[e] = (1.0 - math.sqrt(1.0 - 4.0 * a * c)) / (2.0 * a)
        self.uE = uE
        self.denom = 1.0 - 2.0 * self.PD * uE
        self.gsp = {}
        for v, nm in gt.leafname.items():
            if nm not in gene2sp or gene2sp[nm] not in sp.leaves:
                sys.exit(f"gene {nm} has no species-tree leaf")
            self.gsp[v] = sp.leaves[gene2sp[nm]]
        self.clv = {}    # (a,b) directed: subtree at b seen from a -> (vec, logscale, lognorm)

    def children(self, a, b):
        if b in self.gt.leafname:
            return None
        k = [c for c in self.gt.adj[b] if c != a]
        return (b, k[0]), (b, k[1])

    def _combine(self, vl, vr, vec_self_needed=True, leaf_species=None):
        sp = self.sp
        L, R = sp.left, sp.right
        PS, PD = self.PS, self.PD
        out = np.zeros(sp.n)
        if vl is not None:
            out += PD * vl * vr
            I = sp.internal
            out[I] += PS * (vl[L[I]] * vr[R[I]] + vl[R[I]] * vr[L[I]])
            leaf_mask = L < 0
            out[leaf_mask] /= self.denom[leaf_mask]
        else:
            out[leaf_species] = PS     # gene leaf on its species leaf (no DL correction)
        uE = self.uE
        for lev in sp.levels:          # add SL, children first
            out[lev] = (out[lev] + PS * (out[L[lev]] * uE[R[lev]] + out[R[lev]] * uE[L[lev]])) / self.denom[lev]
        return out

    def get(self, a, b):
        key = (a, b)
        if key in self.clv:
            return self.clv[key]
        stack = [key]
        while stack:
            k = stack[-1]
            if k in self.clv:
                stack.pop(); continue
            ch = self.children(*k)
            if ch is None:
                vec = self._combine(None, None, leaf_species=self.gsp[k[1]])
                m = vec.max()
                self.clv[k] = (vec / m, math.log(m), math.log(m))
                stack.pop(); continue
            todo = [c for c in ch if c not in self.clv]
            if todo:
                stack.extend(todo); continue
            (vl, sl, _), (vr, sr, _) = self.clv[ch[0]], self.clv[ch[1]]
            vec = self._combine(vl, vr)
            m = vec.max()
            if m <= 0:
                sys.exit("zero likelihood clade (gene/species mapping inconsistent?)")
            self.clv[k] = (vec / m, sl + sr + math.log(m), math.log(m))
            stack.pop()
        return self.clv[key]

    def root_vec(self, a, b):
        (vl, sl, _), (vr, sr, _) = self.get(a, b), self.get(b, a)
        vec = self._combine(vl, vr)
        m = vec.max()
        return vec / m, sl + sr + math.log(m), math.log(m)

    def sum_loglik(self):
        """log of sum over gene-tree roots and species nodes of P(virtual root, e): the reconciliation
        log-likelihood GeneRax reports (matches GeneRax 2.1.3's printed ll to 0.01 log units)."""
        tot = []
        for (a, b) in self.gt.edges:
            vec, ls, _ = self.root_vec(a, b)
            tot.append(ls + math.log(vec.sum()))
        mx = max(tot)
        return mx + math.log(sum(math.exp(t - mx) for t in tot))

    def infer(self):
        best = None
        for (a, b) in self.gt.edges:
            vec, ls, _ = self.root_vec(a, b)
            e = int(np.argmax(vec))
            ll = ls + math.log(vec[e])
            if best is None or ll > best[0] + 1e-12:
                best = (ll, (a, b), e)
        # total (sum over roots and species) log-likelihood, uniform origination, for reporting
        tot = []
        for (a, b) in self.gt.edges:
            vec, ls, _ = self.root_vec(a, b)
            tot.append(ls + math.log(vec.sum()))
        mx = max(tot)
        self.loglik_sum = mx + math.log(sum(math.exp(t - mx) for t in tot))
        self.best_root = best
        return self.backtrace(best[1], best[2])

    def backtrace(self, root_edge, e_root):
        sp, PS, PD, uE = self.sp, self.PS, self.PD, self.uE
        L, R = sp.left, sp.right
        counts = defaultdict(lambda: [0, 0, 0])   # presence(S+SL+leaf), D, SL
        tot = Counter()
        node_event = {}                           # gene clade -> (event, species index)
        a, b = root_edge
        rv = self.root_vec(a, b)
        stack = [("ROOT", (a, b), e_root)]
        while stack:
            kind, u, e = stack.pop()
            if kind == "ROOT":
                vu, _, nu = rv
                lch, rch = self.get(*u), self.get(u[1], u[0])
                l_id, r_id = u, (u[1], u[0])
                ch = (l_id, r_id)
                clade_id = ("root", u)
            else:
                vu, _, nu = self.get(*u)
                ch = self.children(*u)
                clade_id = u
            is_leaf_gene = ch is None
            is_leaf_sp = L[e] < 0
            if is_leaf_gene and is_leaf_sp:
                assert self.gsp[u[1]] == e
                counts[e][0] += 1; tot["Leaf"] += 1
                node_event[clade_id] = ("Leaf", e)
                continue
            vals = [0.0] * 5
            if not is_leaf_gene:
                vl = self.get(*ch[0])[0]
                vr = self.get(*ch[1])[0]
                if not is_leaf_sp:
                    f, g = L[e], R[e]
                    vals[0] = PS * vl[f] * vr[g]
                    vals[1] = PS * vl[g] * vr[f]
                vals[2] = PD * vl[e] * vr[e]
                k = math.exp(nu)            # vu is in units of (scale_l*scale_r*norm_u)
            else:
                k = None
            if not is_leaf_sp:
                f, g = L[e], R[e]
                if is_leaf_gene:
                    vals[3] = PS * vu[f] * uE[g]
                    vals[4] = PS * vu[g] * uE[f]
                else:
                    vals[3] = PS * vu[f] * uE[g] * k
                    vals[4] = PS * vu[g] * uE[f] * k
            i = int(np.argmax(vals))
            if vals[i] <= 0:
                sys.exit(f"invalid scenario at gene clade {u}, species {sp.names[e]}")
            if i in (0, 1):
                counts[e][0] += 1; tot["S"] += 1
                node_event[clade_id] = ("S", e)
                f, g = L[e], R[e]
                if i == 0:
                    stack += [("N", ch[0], f), ("N", ch[1], g)]
                else:
                    stack += [("N", ch[0], g), ("N", ch[1], f)]
            elif i == 2:
                counts[e][1] += 1; tot["D"] += 1
                node_event[clade_id] = ("D", e)
                stack += [("N", ch[0], e), ("N", ch[1], e)]
            else:
                counts[e][0] += 1; counts[e][2] += 1; tot["SL"] += 1
                stack.append((kind, u, L[e] if i == 3 else R[e]))
        self.counts, self.tot, self.node_event = counts, tot, node_event
        return counts, tot

    def leaves_of(self, a, b):
        out, stack = [], [(a, b)]
        while stack:
            x, y = stack.pop()
            if y in self.gt.leafname:
                out.append(self.gt.leafname[y])
            else:
                stack += [(y, c) for c in self.gt.adj[y] if c != x]
        return out

    def clade_tips(self, cid):
        if cid[0] == "root":
            a, b = cid[1]
            return self.leaves_of(a, b) + self.leaves_of(b, a)
        return self.leaves_of(*cid)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gene-tree", required=True)
    ap.add_argument("--species-tree", required=True)
    ap.add_argument("--map", default="data/species_tree/gene_species.map")
    ap.add_argument("--dup", type=float, default=0.2)
    ap.add_argument("--loss", type=float, default=0.2)
    ap.add_argument("--optimize-rates", action="store_true",
                    help="estimate global D/L rates by ML first (GeneRax default behaviour), then reconcile")
    ap.add_argument("--anchors", default="P07271_Saccharomyces,P38766_Saccharomyces")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    gene2sp = {}
    for line in open(args.map):
        p = line.split()
        if len(p) >= 2:
            gene2sp[p[0]] = p[1]
    sp = Species(args.species_tree)
    gt = GeneTree(args.gene_tree, keep=[g for g in gene2sp if gene2sp[g] in sp.leaves])
    if args.optimize_rates:
        args.dup, args.loss, ll_opt, nev = optimise_rates(sp, gt, gene2sp, args.dup, args.loss)
        print(f"[35] ML rates D={args.dup:.6f} L={args.loss:.6f} ll={ll_opt:.3f} ({nev} evaluations)",
              file=sys.stderr)
    m = UndatedDL(sp, gt, gene2sp, args.dup, args.loss)
    counts, tot = m.infer()

    os.makedirs(args.outdir, exist_ok=True)
    with open(os.path.join(args.outdir, "speciesEventCounts.txt"), "w") as fh:
        for e in range(sp.n):
            c = counts.get(e, [0, 0, 0])
            fh.write(f"{sp.names[e]} {c[0]} {c[1]} {c[2]} 0 \n")
    with open(os.path.join(args.outdir, "eventCounts.txt"), "w") as fh:
        for k in ("S", "SL", "D", "Leaf"):
            fh.write(f"{k}:{tot.get(k, 0)}\n")

    # anchor duplication: smallest reconciled gene clade holding both anchors
    a1, a2 = args.anchors.split(",")
    anchor = None
    for cid, (ev, e) in m.node_event.items():
        if ev == "Leaf":
            continue
        tips = m.clade_tips(cid)
        if a1 in tips and a2 in tips and (anchor is None or len(tips) < anchor[0]):
            anchor = (len(tips), ev, sp.names[e], len(sp.desc[e]))
    topD = sorted(((c[1], sp.names[e], len(sp.desc[e])) for e, c in counts.items()), reverse=True)[:10]
    summ = {"gene_tree": args.gene_tree, "species_tree": args.species_tree, "rates": [args.dup, args.loss], "rates_optimised": args.optimize_rates,
            "n_genes": len(gt.leafname), "n_species_leaves": len(sp.leaves),
            "totals": dict(tot), "ml_scenario_loglik": m.best_root[0], "sum_loglik": m.loglik_sum,
            "root_species": sp.names[m.best_root[2]],
            "anchor_duplication": {"gene_clade_size": anchor[0], "event": anchor[1],
                                   "species_node": anchor[2], "species_node_leaves": anchor[3]} if anchor else None,
            "top_D": [{"node": n, "D": d, "leaves": k} for d, n, k in topD]}
    json.dump(summ, open(os.path.join(args.outdir, "summary.json"), "w"), indent=1)
    if not args.quiet:
        print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
