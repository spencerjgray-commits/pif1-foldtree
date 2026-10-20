#!/usr/bin/env python3
"""
34_core_plddt_fixed.py  --  recompute per-tip helicase-core mean pLDDT (the R2 filter) without the
core-only slicing bug, and write the corrected R2 alignment subset.

Bug (2026-10-15 review, R2 row / CRIT-K10): results/reviewer/core_plddt.tsv sliced residues
core_from..core_to (UniProt numbering) out of every ColabFold PDB. Three ColabFold models are
core-only and renumbered 1..N, so for A0ABX8I2X5 / A0A2V1AIT0 (core 166-369, file 1-204) the slice
hit residues 166-204 only: 72.81 instead of the true core mean 88.73. Both then fell below 87 and
were dropped from R2 (748 tips instead of 750).

Here:
  AFDB tips      -> CA B-factors (= pLDDT) from data/structures/afdb/AF-<acc>-F1-model_v6.cif,
                    residues core_from..core_to (label_seq_id is UniProt numbering).
  ColabFold tips -> CA B-factors from data/structures/af3/<acc>.pdb; if the model length equals the
                    full-length protein, slice core_from..core_to; if it is a core-only model
                    (data/structures/af3/_cores_note.txt), average the whole model.

Outputs:
  results/rerun_2026-10/robustness/core_plddt_fixed.tsv   (old value alongside for comparison)
  results/rerun_2026-10/robustness/aln_plddt87_fixed.fasta (alnAA_3di.fasta rows with core pLDDT >= 87)

Run (repo root):  python workflow/34_core_plddt_fixed.py
"""
import argparse
import csv
import os
import sys


def cif_ca_plddt(path):
    cols, out, in_loop = [], {}, False
    with open(path) as fh:
        for line in fh:
            if line.startswith("_atom_site."):
                cols.append(line.strip().split(".", 1)[1])
                in_loop = True
                continue
            if in_loop and line.startswith(("ATOM", "HETATM")):
                f = line.split()
                rec = dict(zip(cols, f))
                if rec.get("label_atom_id") == "CA":
                    out[int(rec["label_seq_id"])] = float(rec["B_iso_or_equiv"])
            elif in_loop and cols and not line.startswith(("ATOM", "HETATM")) and out:
                break
    return out


def pdb_ca_plddt(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                out[int(line[22:26])] = float(line[60:66])
    return out


def mean(v):
    return sum(v) / len(v) if v else None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--aln", default="results/seq_tree/alnAA_3di.fasta")
    ap.add_argument("--tip-map", default="data/seqs/tip_map.tsv")
    ap.add_argument("--manifest", default="manifest.csv")
    ap.add_argument("--old", default="results/reviewer/core_plddt.tsv")
    ap.add_argument("--core-note", default="data/structures/af3/_cores_note.txt")
    ap.add_argument("--threshold", type=float, default=87.0)
    ap.add_argument("--outdir", default="results/rerun_2026-10/robustness")
    args = ap.parse_args()

    tips = {r["tip_label"]: r for r in csv.DictReader(open(args.tip_map), delimiter="\t")}
    src = {r["accession"]: r["structure_source"] for r in csv.DictReader(open(args.manifest))}
    old = {r["tip"]: r["core_plddt"] for r in csv.DictReader(open(args.old), delimiter="\t")}
    core_only = {r["accession"] for r in csv.DictReader(open(args.core_note), delimiter="\t")}

    aln = {}
    name = None
    for line in open(args.aln):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            aln[name] = []
        elif name:
            aln[name].append(line)

    rows = []
    for tip in aln:
        t = tips[tip]
        acc = t["accession"]
        lo, hi, full = int(t["core_from"]), int(t["core_to"]), int(t["full_len"])
        s = src.get(acc, "")
        if s == "AFDB":
            ca = cif_ca_plddt(f"data/structures/afdb/AF-{acc}-F1-model_v6.cif")
            how = "afdb_cif_core_slice"
            val = mean([ca[i] for i in range(lo, hi + 1) if i in ca])
        else:
            ca = pdb_ca_plddt(f"data/structures/af3/{acc}.pdb")
            if acc in core_only:
                how = "colabfold_core_only_whole_model"
                val = mean(list(ca.values()))
            elif len(ca) == full:
                how = "colabfold_core_slice"
                val = mean([ca[i] for i in range(lo, hi + 1) if i in ca])
            else:
                how = f"colabfold_len_mismatch({len(ca)}vs{full})"
                val = mean([ca[i] for i in range(lo, hi + 1) if i in ca])
        rows.append({"tip": tip, "accession": acc, "source": s, "method": how,
                     "core_plddt_fixed": f"{val:.4f}" if val is not None else "NA",
                     "core_plddt_old": old.get(tip, "NA")})

    os.makedirs(args.outdir, exist_ok=True)
    with open(os.path.join(args.outdir, "core_plddt_fixed.tsv"), "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    keep = [r["tip"] for r in rows if r["core_plddt_fixed"] != "NA" and float(r["core_plddt_fixed"]) >= args.threshold]
    with open(os.path.join(args.outdir, "aln_plddt87_fixed.fasta"), "w") as fh:
        for tip in keep:
            fh.write(f">{tip}\n" + "".join(aln[tip]) + "\n")

    changed = [r for r in rows if r["core_plddt_old"] != "NA" and r["core_plddt_fixed"] != "NA"
               and abs(float(r["core_plddt_old"]) - float(r["core_plddt_fixed"])) > 0.01]
    sys.stderr.write(f"[34] {len(rows)} tips; {len(keep)} with core pLDDT >= {args.threshold:g}\n")
    sys.stderr.write(f"[34] {len(changed)} tips differ from {args.old} by > 0.01:\n")
    for r in changed:
        sys.stderr.write(f"     {r['tip']}  old {r['core_plddt_old']}  fixed {r['core_plddt_fixed']}  ({r['method']})\n")
    near = [r for r in rows if r["core_plddt_fixed"] != "NA" and abs(float(r["core_plddt_fixed"]) - args.threshold) < 0.05]
    for r in near:
        sys.stderr.write(f"[34] at threshold: {r['tip']} {r['core_plddt_fixed']}\n")


if __name__ == "__main__":
    main()
