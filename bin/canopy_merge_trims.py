#!/usr/bin/env python3
"""Combine the adapter trims of two writing rounds into one log, and correct the stats to match.

A second round trims a genome that was already trimmed, so its log and stats describe only that
round. Each output row gives one contig's total trim, its final kept length and its final decision.
"""
import argparse
import csv

ap = argparse.ArgumentParser()
ap.add_argument("--prior", required=True, help="trim log of the first round")
ap.add_argument("--new", required=True, help="trim log of the second round; rewritten in place")
ap.add_argument("--decisions", required=True, help="final decision table")
ap.add_argument("--stats", required=True, help="stats of the second round; rewritten in place")
a = ap.parse_args()


def rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


with open(a.decisions) as fh:
    decision = {r["name"].replace("+", "plus"): r["decision"] for r in csv.DictReader(fh, delimiter="\t")}

fields = ["sample", "sequence", "contig", "trimmed_start_bp", "trimmed_end_bp", "kept_bp", "decision"]
merged = {}
for r in rows(a.prior) + rows(a.new):
    key = (r["sequence"], r["contig"])
    if key in merged:
        merged[key]["trimmed_start_bp"] = int(merged[key]["trimmed_start_bp"]) + int(r["trimmed_start_bp"])
        merged[key]["trimmed_end_bp"] = int(merged[key]["trimmed_end_bp"]) + int(r["trimmed_end_bp"])
        merged[key]["kept_bp"] = r["kept_bp"]
    else:
        merged[key] = dict(r)
    merged[key]["decision"] = decision[r["sequence"]]

with open(a.new, "w", newline="") as fh:
    w = csv.DictWriter(fh, fields, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(merged.values())

kept = [r for r in merged.values() if r["decision"] != "remove"]
updates = {
    "adapter_contig_ends_trimmed": len(kept),
    "adapter_bases_trimmed": sum(int(r["trimmed_start_bp"]) + int(r["trimmed_end_bp"]) for r in kept),
}
with open(a.stats) as fh:
    lines = [line.rstrip("\n").split("=", 1) for line in fh if "=" in line]
with open(a.stats, "w") as fh:
    for key, value in lines:
        fh.write(f"{key}={updates.get(key, value)}\n")
