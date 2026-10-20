#!/usr/bin/env python3
"""
36_compare_event_counts.py  --  compare two per-species event-count files (GeneRax
pif1_speciesEventCounts.txt layout: <node> <S+SL+leaf> <D> <SL> <T>) row by row.

Used to validate workflow/35_undated_dl.py against archived GeneRax 2.0.4 runs.

Usage: python workflow/36_compare_event_counts.py OLD NEW [--nodes node_448,node_578,node_183,node_0]
"""
import argparse


def load(p):
    out = {}
    for line in open(p):
        f = line.split()
        if len(f) >= 4 and not f[0].startswith("#"):
            out[f[0]] = tuple(int(x) for x in f[1:4])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("old")
    ap.add_argument("new")
    ap.add_argument("--nodes", default="node_448,node_578,node_183,node_0")
    a = ap.parse_args()
    o, n = load(a.old), load(a.new)
    keys = sorted(set(o) | set(n))
    diff = [k for k in keys if o.get(k) != n.get(k)]
    tot = lambda d, i: sum(v[i] for v in d.values())
    print(f"rows old/new: {len(o)}/{len(n)}; rows differing: {len(diff)}")
    print(f"totals old (pres,D,SL): {tot(o,0)},{tot(o,1)},{tot(o,2)}   new: {tot(n,0)},{tot(n,1)},{tot(n,2)}")
    for k in a.nodes.split(","):
        print(f"  {k}: old {o.get(k)}  new {n.get(k)}")
    for k in diff[:10]:
        print(f"  differs: {k} old {o.get(k)} new {n.get(k)}")


if __name__ == "__main__":
    main()
