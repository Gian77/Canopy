#!/usr/bin/env python3
"""Apply Canopy's contaminant decision rule to every sequence in the scaffolded genome.

Only unplaced contigs can be removed. A contig is removed when any of these holds:
  R1_fcs_exclude      FCS-GX action is EXCLUDE
  R2_fcs_review       FCS-GX action is REVIEW and another source agrees it is foreign
  R3_kraken_blast     Kraken/BlobTools rule said remove and BLASTn corroborates it
  R4_fungal_alignment >= min_fungal of the contig aligns to a fungal reference and
                      < max_ref of it aligns to the nuclear reference
  R5_blast            BLASTn corroborates a non-plant hit and < max_ref aligns to the reference

GC and coverage are recorded for interpretation but never trigger removal: in the five sorghum
genomes about 19% of retained unplaced contigs >= 10 kb fell outside the plant GC range, and
contaminant coverage overlapped plant coverage.
"""
import argparse
import collections
import csv

ap = argparse.ArgumentParser()
ap.add_argument("--sample", required=True)
ap.add_argument("--audit", required=True, help="CLASSIFY_CONTAMINANTS audit table")
ap.add_argument("--blob-table", required=True)
ap.add_argument("--fcs-report", required=True)
ap.add_argument("--fungal", required=True)
ap.add_argument("--blast", required=True)
ap.add_argument("--min-fungal", type=float, default=0.5)
ap.add_argument("--max-ref", type=float, default=0.10)
a = ap.parse_args()


def table(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def unsafe_name(name):
    """FCS-GX rejects '+' in IDs, so screening replaces it with 'plus'; undo that here."""
    return name.replace("plus", "+")


cov = {}
with open(a.blob_table) as fh:
    for line in fh:
        if line.startswith("##") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if line.startswith("#"):
            cov_col = next(i for i, c in enumerate(f) if c.startswith("bam") or c == "cov_sum")
            continue
        cov[f[0]] = f[cov_col]

fcs = collections.defaultdict(list)
with open(a.fcs_report) as fh:
    for line in fh:
        if not line.startswith("#"):
            f = line.rstrip("\n").split("\t")
            fcs[unsafe_name(f[0])].append(f)

fungal = {r["contig"]: r for r in table(a.fungal)}
blast = {r["contig"]: r for r in table(a.blast)}
audit = table(a.audit)
known = {r["name"] for r in audit}
missing = sorted(set(fcs) - known)
if missing:
    raise SystemExit(f"ERROR: FCS-GX reported sequences absent from the audit: {missing[:5]}")

extra = ["mean_cov", "rule", "kraken_decision", "fcs_action", "fcs_div", "fcs_top_tax",
         "fungal_species", "fungal_frac", "blast_identity", "blast_qcov",
         "blast_corroborated", "blast_title"]
base = ["name", "length", "GC", "phylum", "placed", "flagged", "ref_aligned_frac", "decision"]

n_removed = bases_removed = n_review = review_bases = placed_fcs = 0
by_rule = collections.Counter()
by_div = collections.Counter()

with open(f"{a.sample}_contig_decisions.tsv", "w", newline="") as out, \
        open(f"{a.sample}_contaminant_ids.txt", "w") as ids, \
        open(f"{a.sample}_review_ids.txt", "w") as rev:
    w = csv.DictWriter(out, base + extra, delimiter="\t", lineterminator="\n")
    w.writeheader()
    for r in audit:
        n = r["name"]
        unplaced = r["placed"] == "unplaced"
        ref = float(r["ref_aligned_frac"]) if unplaced and r["ref_aligned_frac"] != "NA" else 0.0
        F = fcs.get(n, [])
        actions = {x[4] for x in F}
        u = fungal.get(n, {})
        frac = float(u.get("aligned_frac") or 0)
        b = blast.get(n, {})
        corroborated = b.get("corroborated") == "yes"
        kraken_remove = r["decision"] == "remove"
        fungal_hit = frac >= a.min_fungal and ref < a.max_ref
        blast_hit = corroborated and ref < a.max_ref

        rule = ""
        if unplaced:
            if "EXCLUDE" in actions:
                rule = "R1_fcs_exclude"
            elif "REVIEW" in actions and (kraken_remove or fungal_hit or corroborated):
                rule = "R2_fcs_review"
            elif kraken_remove and corroborated:
                rule = "R3_kraken_blast"
            elif fungal_hit:
                rule = "R4_fungal_alignment"
            elif blast_hit:
                rule = "R5_blast"

        unresolved = bool(F) or r["decision"] in ("review", "remove")
        decision = "remove" if rule else "review" if unresolved else "keep"
        length = int(r["length"])
        if F and not unplaced:
            placed_fcs += 1
        if decision == "remove":
            n_removed += 1
            bases_removed += length
            by_rule[rule] += 1
            by_div[F[0][5] if F else "other_evidence"] += 1
            ids.write(n + "\n")
        elif decision == "review":
            n_review += 1
            review_bases += length
            rev.write(n + "\n")

        row = {k: r[k] for k in base}
        row.update(decision=decision, flagged="yes" if decision == "remove" else "no",
                   mean_cov=cov.get(n, "NA"), rule=rule, kraken_decision=r["decision"],
                   fcs_action=";".join(x[4] for x in F), fcs_div=";".join(x[5] for x in F),
                   fcs_top_tax=";".join(x[7] for x in F),
                   fungal_species=u.get("species", ""), fungal_frac=u.get("aligned_frac", ""),
                   blast_identity=b.get("best_identity", ""), blast_qcov=b.get("best_qcov", ""),
                   blast_corroborated=b.get("corroborated", ""), blast_title=b.get("best_title", ""))
        w.writerow(row)

with open(f"{a.sample}_contamination_summary.txt", "w") as s:
    s.write(f"contaminants_flagged={n_removed}\n")
    s.write(f"bases_removed={bases_removed}\n")
    s.write("phyla_detected=" + ",".join(f"{k}:{v}" for k, v in by_div.most_common()) + "\n")
    s.write(f"review_contigs={n_review}\n")
    s.write(f"review_bases={review_bases}\n")
    s.write(f"placed_with_fcs_call={placed_fcs}\n")
    s.write("removed_by_rule=" + ",".join(f"{k}:{v}" for k, v in sorted(by_rule.items())) + "\n")
