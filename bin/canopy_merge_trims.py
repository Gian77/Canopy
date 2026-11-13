#!/usr/bin/env python3
"""Combine the adapter trims of two writing rounds into one log, and correct the stats to match.

A second round trims a genome that was already trimmed, so its log and stats describe only that
round. Each output row gives one contig's total trim, its final kept length and its final decision:
the decision table's, or too_short for a sequence the length filter left out of the final FASTA.
"""
import argparse
import csv

ap = argparse.ArgumentParser()
ap.add_argument("--prior", required=True, help="trim log of the first round")
ap.add_argument("--new", required=True, help="trim log of the second round; rewritten in place")
ap.add_argument("--decisions", required=True, help="final decision table")
ap.add_argument("--stats", required=True, help="stats of the second round; rewritten in place")
ap.add_argument("--final", required=True, help="final FASTA of the second round")
a = ap.parse_args()

with open(a.final) as fh:
    in_final = {line[1:].split()[0] for line in fh if line.startswith(">")}


def rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


with open(a.decisions) as fh:
    decision = {r["name"].replace("+", "plus"): r["decision"] for r in csv.DictReader(fh, delimiter="\t")}

fields = ["sample", "sequence", "contig", "trimmed_start_bp", "trimmed_end_bp", "kept_bp", "decision", "internal_cut_bp"]
merged = {}
for r in rows(a.prior) + rows(a.new):
    key = (r["sequence"], r["contig"])
    if key in merged:
        merged[key]["trimmed_start_bp"] = int(merged[key]["trimmed_start_bp"]) + int(r["trimmed_start_bp"])
        merged[key]["trimmed_end_bp"] = int(merged[key]["trimmed_end_bp"]) + int(r["trimmed_end_bp"])
        merged[key]["internal_cut_bp"] = int(merged[key]["internal_cut_bp"]) + int(r["internal_cut_bp"])
        merged[key]["kept_bp"] = r["kept_bp"]
    else:
        merged[key] = dict(r)
    d = decision[r["sequence"]]
    merged[key]["decision"] = "too_short" if d != "remove" and r["sequence"] not in in_final else d

with open(a.new, "w", newline="") as fh:
    w = csv.DictWriter(fh, fields, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(merged.values())

kept = [r for r in merged.values() if r["decision"] not in ("remove", "too_short")]
updates = {
    "adapter_contig_ends_trimmed": sum(1 for r in kept if int(r["trimmed_start_bp"]) + int(r["trimmed_end_bp"])),
    "adapter_bases_trimmed": sum(int(r["trimmed_start_bp"]) + int(r["trimmed_end_bp"]) for r in kept),
    "adapter_internal_cuts": sum(1 for r in kept if int(r["internal_cut_bp"])),
    "adapter_internal_bases_cut": sum(int(r["internal_cut_bp"]) for r in kept),
}
with open(a.stats) as fh:
    lines = [line.rstrip("\n").split("=", 1) for line in fh if "=" in line]
with open(a.stats, "w") as fh:
    for key, value in lines:
        fh.write(f"{key}={updates.get(key, value)}\n")
