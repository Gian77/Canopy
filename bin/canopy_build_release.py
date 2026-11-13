#!/usr/bin/env python3
"""Assemble one self-contained release directory from finished Canopy runs.

  <release>/README.md                how to use and interpret the release (--readme)
  <release>/summary_all_samples.tsv  one row per genome
  <release>/MANIFEST.sha256          checksums of every file
  <release>/<sample>/genome/         final and flagged nuclear genomes, AGPs, organelles
  <release>/<sample>/evidence/       per-contig decisions, screens, BLAST, Sylph, blob plots
  <release>/<sample>/qc/             assembly summary and pipeline_tools (Markdown and HTML), QUAST, BUSCO

Files come from each run's reports/<sample>/<sample>_assembly.zip, which PACKAGE_RESULTS builds from
the task outputs. Published copies under assembly/final/ are not used: a run aborted during
publishing can leave a truncated copy there, and a resumed run does not repair it. Each genome is
checked against its stats before the release is written.
"""
import argparse
import csv
import hashlib
import os
import shutil
import zipfile

from canopy_report_html import convert as report_to_html

ap = argparse.ArgumentParser()
ap.add_argument("--results-pattern", default="results_{sample}_v4")
ap.add_argument("--project", default=".")
ap.add_argument("--release", required=True)
ap.add_argument("--samples", nargs="+",
                default=["F25101", "F07020", "F10702", "BTx2932", "B11077"])
ap.add_argument("--extra-evidence", default=None,
                help="directory with one per-contig evidence table per sample: <sample>.contig_evidence.tsv")
ap.add_argument("--readme", default=None, help="file to install as <release>/README.md")
a = ap.parse_args()


def keyvals(path):
    with open(path) as fh:
        return dict(line.rstrip("\n").split("=", 1) for line in fh if "=" in line)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def fasta_size(path):
    """Length of each sequence."""
    lengths = {}
    name = None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                lengths[name] = 0
            else:
                lengths[name] += len(line.strip())
    return lengths


def check_genome(sample, fasta, agp, expected=None):
    lengths = fasta_size(fasta)
    ends = {}
    with open(agp) as fh:
        for line in fh:
            if line.strip() and not line.startswith("#"):
                f = line.split("\t")
                ends[f[0]] = int(f[2])
    if ends != lengths:
        raise SystemExit(f"ERROR: {sample}: {fasta} does not match its AGP")
    if expected:
        got = (len(lengths), sum(lengths.values()))
        want = (int(expected["final_sequences"]), int(expected["final_bases"]))
        if got != want:
            raise SystemExit(f"ERROR: {sample}: {fasta} has {got}, stats say {want}")


# Where each file in the zip belongs. Anything unlisted stops the build, so nothing is dropped unseen.
def destination(sample, rel):
    name = os.path.basename(rel)
    if rel.startswith("screening/"):
        return os.path.join("evidence", "screening", name)
    if rel.startswith("blobtools/"):
        return os.path.join("evidence", "blobtools", name)
    if name in (f"{sample}_final.fasta", f"{sample}_final.agp", f"{sample}_full_flagged.fasta",
                f"{sample}_full_flagged.agp", f"{sample}_chloroplast.fasta", f"{sample}_mitochondria.fasta"):
        return os.path.join("genome", name)
    if name in (f"{sample}_contig_decisions.tsv", f"{sample}_adapter_trims.tsv", f"{sample}_placed_fcs_cuts.tsv",
                f"{sample}_blast_corroboration.tsv", f"{sample}_sylph_corroboration.tsv"):
        return os.path.join("evidence", name)
    if name.startswith("quast_") or name in (f"{sample}_assembly_summary.md", "busco_nuclear_final_summary.txt",
                                             "pipeline_tools.md", "README.txt"):
        return os.path.join("qc", name)
    raise SystemExit(f"ERROR: {sample}: no place in the release for {rel}")


rows = []
for s in a.samples:
    R = os.path.join(a.project, a.results_pattern.format(sample=s))
    O = os.path.join(a.release, s)
    with zipfile.ZipFile(os.path.join(R, "reports", s, f"{s}_assembly.zip")) as z:
        for info in z.infolist():
            if info.is_dir():
                continue
            rel = info.filename.split("/", 1)[1]
            dst = os.path.join(O, destination(s, rel))
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            with z.open(info) as src, open(dst, "wb") as out:
                shutil.copyfileobj(src, out)
    # In the zip the summary sits next to blobtools/; in the release it is in qc/ and the plots are in
    # evidence/blobtools/, so point its image links there before rendering (the HTML embeds them).
    summary = os.path.join(O, "qc", f"{s}_assembly_summary.md")
    with open(summary, encoding="utf-8") as fh:
        text = fh.read()
    with open(summary, "w", encoding="utf-8") as fh:
        fh.write(text.replace("](blobtools/", "](../evidence/blobtools/"))
    report_to_html(summary)
    report_to_html(os.path.join(O, "qc", "pipeline_tools.md"))

    # The blob plots come with the zip (evidence/blobtools/); they are not copied a second time.
    if a.extra_evidence:
        shutil.copy2(os.path.join(a.extra_evidence, f"{s}.contig_evidence.tsv"), os.path.join(O, "evidence"))

    screen = keyvals(os.path.join(O, "evidence", "screening", f"{s}_ncbi_screen.txt"))
    contam = keyvals(os.path.join(R, "qc", "contamination", s, "decisions", f"{s}_contamination_summary.txt"))
    G = os.path.join(O, "genome")
    check_genome(s, os.path.join(G, f"{s}_final.fasta"), os.path.join(G, f"{s}_final.agp"), screen)
    check_genome(s, os.path.join(G, f"{s}_full_flagged.fasta"), os.path.join(G, f"{s}_full_flagged.agp"))
    with open(os.path.join(O, "evidence", f"{s}_contig_decisions.tsv")) as fh:
        removed = sum(r["decision"] == "remove" for r in csv.DictReader(fh, delimiter="\t"))
    if removed != int(contam["contaminants_flagged"]):
        raise SystemExit(f"ERROR: {s}: decision table and contamination summary disagree")

    rows.append({
        "sample": s,
        "final_bases": screen["final_bases"],
        "final_sequences": screen["final_sequences"],
        "final_contigs": screen["final_contigs"],
        "final_contig_n50": screen["final_contig_n50"],
        "contigs_removed": contam["contaminants_flagged"],
        "bases_removed": contam["bases_removed"],
        "review_contigs": contam["review_contigs"],
        "review_bases": contam["review_bases"],
        "adapter_bases_trimmed": screen["adapter_bases_trimmed"],
        "fcs_gx_exclude_final": screen["gx_exclude"],
        "fcs_gx_other_final": screen["gx_other_actions"],
        "fcs_adaptor_hits_final": screen["adaptor_hits"],
        "screen_status": screen["screen_status"],
        "fcs_gx_db_build": screen["gx_db_build"],
        "organelle_fcs_gx_actions": screen.get("organelle_gx_actions", "not_screened"),
        "organelle_fcs_adaptor_hits": screen.get("organelle_adaptor_hits", "not_screened"),
        "organelle_screen_status": screen.get("organelle_screen_status", "not_screened"),
    })

with open(os.path.join(a.release, "summary_all_samples.tsv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, rows[0].keys(), delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)

if a.readme:
    shutil.copy2(a.readme, os.path.join(a.release, "README.md"))

manifest = os.path.join(a.release, "MANIFEST.sha256")
with open(manifest, "w") as fh:
    for d, _, files in sorted(os.walk(a.release)):
        for f in sorted(files):
            p = os.path.join(d, f)
            if p != manifest:
                fh.write(f"{sha256(p)}  {os.path.relpath(p, a.release)}\n")
