#!/usr/bin/env python3
"""Cut FCS-GX FIX regions out of placed sequences, in the final genome only.

The contaminant rule removes whole unplaced contigs, so a foreign stretch inside a contig on a
chromosome survives it. A FIX region on a placed sequence is cut when under --max-ref of it aligns
to the nuclear reference. A region that does align is left alone and reported, because FCS-GX can be
wrong about sequence that the reference also has.

The region is cut from <sample>_final.fasta and its AGP. A region within --max-end-dist of a contig
end is trimmed to that end; otherwise the contig is split there and the pieces are joined by a 100 bp
gap, named <contig>_fix<N>. <sample>_full_flagged.fasta keeps every base, and the header of the
sequence gains fcs_fix_cut=<start>-<stop>, so coordinates after a cut differ between the two files.

Every FIX region on a placed sequence is logged in <sample>_placed_fcs_cuts.tsv with its outcome, and
the stats file is corrected. Run this after the writer and the trim merge of the same round.
"""
import argparse
import csv

ap = argparse.ArgumentParser()
ap.add_argument("--sample", required=True)
ap.add_argument("--table", required=True, help="FIX regions and their reference alignment")
ap.add_argument("--decisions", required=True, help="decision table with the IDs used in the FASTAs")
ap.add_argument("--screened-agp", required=True, help="AGP of the genome the regions were screened on")
ap.add_argument("--final-fasta", required=True, help="rewritten in place")
ap.add_argument("--final-agp", required=True, help="rewritten in place")
ap.add_argument("--flagged-fasta", required=True, help="headers rewritten in place")
ap.add_argument("--flagged-agp", required=True)
ap.add_argument("--stats", required=True, help="rewritten in place")
ap.add_argument("--max-ref", type=float, default=0.10)
ap.add_argument("--max-end-dist", type=int, default=1000)
ap.add_argument("--min-length", type=int, default=0)
a = ap.parse_args()


def read_agp(path):
    agp = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            agp.setdefault(f[0], []).append(f)
    return agp


def read_fasta(path):
    records, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        records[name] = "".join(buf)
    return records


def write_record(fh, header, seq):
    fh.write(f">{header}\n")
    for i in range(0, len(seq), 80):
        fh.write(seq[i:i + 80] + "\n")


def n50(lengths):
    total, running = sum(lengths), 0
    for x in sorted(lengths, reverse=True):
        running += x
        if running * 2 >= total:
            return x
    return 0


def is_gap(row):
    return row[4] in ("N", "U")


def to_contig(row, start, stop):
    """Object coordinates to coordinates on the contig, which survive adapter trims between rounds."""
    s, beg, end = int(row[1]), int(row[6]), int(row[7])
    if row[8] == "-":
        return end - (stop - s), end - (start - s)
    return beg + (start - s), beg + (stop - s)


def row_for(rows, contig, cbeg, cend):
    """The row of this contig overlapping the region, and the region clipped to it (or Nones)."""
    for row in rows:
        if is_gap(row) or not (row[5] == contig or row[5].startswith(f"{contig}_split")):
            continue
        beg, end = int(row[6]), int(row[7])
        if beg <= cend and cbeg <= end:
            return row, max(cbeg, beg), min(cend, end)
    return None, None, None


def offsets(row, cbeg, cend):
    """The region as [x, y) offsets along the row, in the direction of the sequence."""
    beg, end = int(row[6]), int(row[7])
    if row[8] == "-":
        return end - cend, end - cbeg + 1
    return cbeg - beg, cend - beg + 1


def to_cut(x, y, length):
    """Extend a region near a contig end to that end, as adapter trimming does."""
    to_left, to_right = x, length - y
    if min(to_left, to_right) > a.max_end_dist:
        return x, y
    return (0, y) if to_left <= to_right else (x, length)


def kept_spans(length, cuts):
    spans, at = [], 0
    for x, y in sorted(cuts):
        if x > at:
            spans.append((at, x))
        at = max(at, y)
    if at < length:
        spans.append((at, length))
    return spans


screened = read_agp(a.screened_agp)
flagged_agp = read_agp(a.flagged_agp)
final_agp = read_agp(a.final_agp)
with open(a.decisions) as fh:
    placed = {r["name"]: r["placed"] == "placed" for r in csv.DictReader(fh, delimiter="\t")}

log, flagged_cuts, cuts_on_final = [], {}, {}
with open(a.table) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        seq = r["sequence"]
        if not placed.get(seq):
            continue    # whole unplaced contigs belong to the contaminant rule
        start, stop = int(r["start"]), int(r["stop"])
        entry = {"sample": a.sample, "sequence": seq, "screened_start": start, "screened_stop": stop,
                 "bases": stop - start + 1, "contig": "", "contig_start": "", "contig_stop": "",
                 "flagged_start": "", "flagged_stop": "", "fcs_action": r["action"], "fcs_div": r["div"],
                 "fcs_top_tax": r["top_tax"], "ref_aligned_frac": r["aligned_frac"], "outcome": ""}
        log.append(entry)
        if float(r["aligned_frac"]) >= a.max_ref:
            entry["outcome"] = "kept_aligns_to_reference"
            continue
        comp = next((x for x in screened.get(seq, []) if not is_gap(x)
                     and int(x[1]) <= start and stop <= int(x[2])), None)
        if comp is None:
            entry["outcome"] = "kept_spans_contigs"
            continue
        contig = comp[5]
        cbeg, cend = to_contig(comp, start, stop)
        entry.update(contig=contig, contig_start=cbeg, contig_stop=cend)
        fl_row, fl_beg, fl_end = row_for(flagged_agp.get(seq, []), contig, cbeg, cend)
        fin_row, fin_beg, fin_end = row_for(final_agp.get(seq, []), contig, cbeg, cend)
        if fl_row is None or fin_row is None:
            entry["outcome"] = "kept_not_in_final"
            continue
        x, y = offsets(fl_row, fl_beg, fl_end)
        entry.update(flagged_start=int(fl_row[1]) + x, flagged_stop=int(fl_row[1]) + y - 1)
        flagged_cuts.setdefault(seq, []).append(f"{entry['flagged_start']}-{entry['flagged_stop']}")
        x, y = offsets(fin_row, fin_beg, fin_end)
        length = int(fin_row[2]) - int(fin_row[1]) + 1
        cuts_on_final.setdefault(seq, {}).setdefault(fin_row[3], []).append(to_cut(x, y, length))
        entry["outcome"] = "cut"

too_short = set()
if cuts_on_final:
    final = read_fasta(a.final_fasta)
    for obj, by_part in cuts_on_final.items():
        seq = final[obj]
        parts, rows, pos, part, pending_gap = [], [], 1, 1, None
        for r in final_agp[obj]:
            if is_gap(r):
                pending_gap = r
                continue
            s, e = int(r[1]), int(r[2])
            spans = kept_spans(e - s + 1, by_part.get(r[3], []))
            for i, (x, y) in enumerate(spans):
                gap = ["U", "100", "scaffold", "yes", "unspecified"] if i else (pending_gap[4:] if rows and pending_gap else None)
                if gap:
                    size = int(gap[1])
                    parts.append("N" * size)
                    rows.append([obj, str(pos), str(pos + size - 1), str(part)] + gap)
                    pos += size
                    part += 1
                parts.append(seq[s - 1 + x:s - 1 + y])
                beg, end = int(r[6]), int(r[7])
                beg, end = (end - y + 1, end - x) if r[8] == "-" else (beg + x, beg + y - 1)
                name = r[5] if len(spans) == 1 else f"{r[5]}_fix{i + 1}"
                rows.append([obj, str(pos), str(pos + y - x - 1), str(part), r[4], name, str(beg), str(end), r[8]])
                pos += y - x
                part += 1
            if spans:
                pending_gap = None
        final[obj] = "".join(parts)
        final_agp[obj] = rows
        if len(final[obj]) != pos - 1:
            raise SystemExit(f"ERROR: rebuilt {obj} has inconsistent length")
        if not final[obj] or len(final[obj]) < a.min_length:
            too_short.add(obj)
    for obj in too_short:
        del final[obj], final_agp[obj]

    with open(a.final_fasta, "w") as fh:
        for obj, seq in final.items():
            write_record(fh, obj, seq)
    with open(a.final_agp, "w") as fh:
        fh.write("##agp-version 2.1\n")
        for rows in final_agp.values():
            for r in rows:
                fh.write("\t".join(r) + "\n")

    with open(a.flagged_fasta) as fh:
        lines = fh.readlines()
    with open(a.flagged_fasta, "w") as fh:
        for line in lines:
            if line.startswith(">") and line[1:].split()[0] in flagged_cuts:
                line = f"{line.rstrip()};fcs_fix_cut={','.join(flagged_cuts[line[1:].split()[0]])}\n"
            fh.write(line)

    with open(a.stats) as fh:
        stats = dict(line.rstrip("\n").split("=", 1) for line in fh if "=" in line)
    lengths = [int(r[2]) - int(r[1]) + 1 for rows in final_agp.values() for r in rows if not is_gap(r)]
    stats.update(final_sequences=len(final), final_bases=sum(len(s) for s in final.values()),
                 final_contigs=len(lengths), final_contig_n50=n50(lengths))
else:
    with open(a.stats) as fh:
        stats = dict(line.rstrip("\n").split("=", 1) for line in fh if "=" in line)

done = [e for e in log if e["outcome"] == "cut"]
stats["placed_fix_regions_cut"] = len(done)
stats["placed_fix_bases_cut"] = sum(e["bases"] for e in done)
stats["placed_fix_regions_kept"] = len(log) - len(done)
# One line for the report; the coordinates are those of the flagged genome.
stats["placed_fix_detail"] = "; ".join(
    f"{e['sequence']}:{e['flagged_start']}-{e['flagged_stop']} ({e['bases']} bp, {e['fcs_top_tax']}, "
    f"{float(e['ref_aligned_frac']) * 100:.1f}% aligned to the reference)" for e in done)
with open(a.stats, "w") as fh:
    for key, value in stats.items():
        fh.write(f"{key}={value}\n")

fields = ["sample", "sequence", "flagged_start", "flagged_stop", "bases", "contig", "contig_start", "contig_stop",
          "fcs_action", "fcs_div", "fcs_top_tax", "ref_aligned_frac", "screened_start", "screened_stop", "outcome"]
with open(f"{a.sample}_placed_fcs_cuts.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
    w.writeheader()
    w.writerows(log)
