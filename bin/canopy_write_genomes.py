#!/usr/bin/env python3
"""Write the two published genome versions and their AGPs.

  <sample>_final.fasta         contaminants removed; plain headers
  <sample>_full_flagged.fasta  every sequence; decision and evidence in the headers

Both versions have FCS-adaptor hits trimmed from contig ends, so a contig has the same sequence
in both. Trimming happens at the contig level and the AGP is rebuilt, so scaffold gaps stay
100 bp and contig interiors are untouched. A hit far from any contig end stops the run, because
it could mark a chimeric join that trimming would not fix.

Sequences shorter than --min-length after trimming are left out of the final genome, because NCBI
does not accept sequences under 200 bp. They stay in the flagged genome as
canopy_decision=remove;rule=min_length_<N>bp, and their trims are logged with decision too_short.

'+' in sequence IDs becomes 'plus' because NCBI and FCS do not accept it.
"""
import argparse
import csv

ap = argparse.ArgumentParser()
ap.add_argument("--sample", required=True)
ap.add_argument("--fasta", required=True)
ap.add_argument("--agp", required=True)
ap.add_argument("--decisions", required=True)
ap.add_argument("--adaptor-report", required=True)
ap.add_argument("--max-end-dist", type=int, default=1000)
ap.add_argument("--min-length", type=int, default=0, help="shortest sequence kept in the final genome")
a = ap.parse_args()


def safe(name):
    return name.replace("+", "plus")


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


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


seqs = read_fasta(a.fasta)

agp = {}
with open(a.agp) as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        agp.setdefault(f[0], []).append(f)

hits = {}
with open(a.adaptor_report) as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        hits[f[0]] = [tuple(map(int, x.split(".."))) for x in f[3].split(",")]
unknown = set(hits) - {safe(n) for n in seqs}
if unknown:
    raise SystemExit(f"ERROR: adaptor report names unknown sequences: {sorted(unknown)[:5]}")

with open(a.decisions) as fh:
    decisions = {r["name"]: r for r in csv.DictReader(fh, delimiter="\t")}
if set(decisions) != set(seqs):
    raise SystemExit("ERROR: decision table and FASTA describe different sequences")

new_seq, new_agp, contig_lengths = {}, {}, {}
trim_log = []

for obj, rows in agp.items():
    seq = seqs[obj]
    components = [r for r in rows if r[4] not in ("N", "U")]
    left, right = {}, {}
    for start, stop in hits.get(safe(obj), []):
        comp = next((r for r in components if int(r[1]) <= start and stop <= int(r[2])), None)
        if comp is None:
            raise SystemExit(f"ERROR: adaptor hit {obj}:{start}-{stop} is not inside one contig")
        s, e = int(comp[1]), int(comp[2])
        to_left, to_right = start - s, e - stop
        if min(to_left, to_right) > a.max_end_dist:
            raise SystemExit(f"ERROR: adaptor hit {obj}:{start}-{stop} is "
                             f"{min(to_left, to_right)} bp from a contig end")
        key = comp[5]
        if to_left <= to_right:
            left[key] = max(left.get(key, 0), stop - s + 1)
        else:
            right[key] = max(right.get(key, 0), e - start + 1)

    parts, out_rows, pos, part = [], [], 1, 1
    pending_gap = None
    for r in rows:
        if r[4] in ("N", "U"):
            pending_gap = r
            continue
        s, e, key = int(r[1]), int(r[2]), r[5]
        tl, tr = left.get(key, 0), right.get(key, 0)
        kept = e - s + 1 - tl - tr
        if tl or tr:
            trim_log.append((a.sample, safe(obj), safe(key), tl, tr, max(kept, 0),
                             decisions[obj]["decision"]))
        if kept <= 0:
            continue
        if out_rows and pending_gap is not None:
            size = int(pending_gap[5])
            parts.append("N" * size)
            out_rows.append([safe(obj), str(pos), str(pos + size - 1), str(part)] + pending_gap[4:])
            pos += size
            part += 1
        pending_gap = None
        parts.append(seq[s - 1 + tl:e - tr])
        beg, end = int(r[6]), int(r[7])
        if r[8] == "-":
            beg, end = beg + tr, end - tl
        else:
            beg, end = beg + tl, end - tr
        out_rows.append([safe(obj), str(pos), str(pos + kept - 1), str(part), r[4], safe(key),
                         str(beg), str(end), r[8]])
        pos += kept
        part += 1
        contig_lengths.setdefault(obj, []).append(kept)
    if out_rows:
        new_seq[obj] = "".join(parts)
        new_agp[obj] = out_rows
        if len(new_seq[obj]) != pos - 1:
            raise SystemExit(f"ERROR: rebuilt {obj} has inconsistent length")

too_short = {o for o in new_seq if decisions[o]["decision"] != "remove" and len(new_seq[o]) < a.min_length}
final_objs = [o for o in new_seq if decisions[o]["decision"] != "remove" and o not in too_short]
trim_log = [row[:6] + ("too_short",) if row[1] in {safe(o) for o in too_short} else row for row in trim_log]

agp_header = "##agp-version 2.1\n"
with open(f"{a.sample}_final.fasta", "w") as ff, open(f"{a.sample}_full_flagged.fasta", "w") as fl, \
        open(f"{a.sample}_final.agp", "w") as af, open(f"{a.sample}_full_flagged.agp", "w") as al:
    af.write(agp_header)
    al.write(agp_header)
    for obj in new_seq:
        d = decisions[obj]
        lines = "".join("\t".join(r) + "\n" for r in new_agp[obj])
        al.write(lines)
        taxon = (d["fcs_top_tax"] or d["fungal_species"] or d["phylum"]).split(";")[0].replace(" ", "_")
        decision, rule = d["decision"], d["rule"] or "none"
        if obj in too_short:
            decision, rule = "remove", f"min_length_{a.min_length}bp"
        write_record(fl, f"{safe(obj)} canopy_decision={decision};rule={rule};"
                         f"evidence_taxon={taxon or 'none'}", new_seq[obj])
        if obj in final_objs:
            af.write(lines)
            write_record(ff, safe(obj), new_seq[obj])

with open(f"{a.sample}_adapter_trims.tsv", "w") as t:
    t.write("sample\tsequence\tcontig\ttrimmed_start_bp\ttrimmed_end_bp\tkept_bp\tdecision\n")
    for row in trim_log:
        t.write("\t".join(map(str, row)) + "\n")

final_lengths = [x for o in final_objs for x in contig_lengths[o]]
trimmed_final = [x for x in trim_log if x[6] not in ("remove", "too_short")]
with open(f"{a.sample}_final_stats.txt", "w") as s:
    s.write(f"final_sequences={len(final_objs)}\n")
    s.write(f"final_bases={sum(len(new_seq[o]) for o in final_objs)}\n")
    s.write(f"final_contigs={len(final_lengths)}\n")
    s.write(f"final_contig_n50={n50(final_lengths)}\n")
    s.write(f"adapter_contig_ends_trimmed={len(trimmed_final)}\n")
    s.write(f"adapter_bases_trimmed={sum(x[3] + x[4] for x in trimmed_final)}\n")
    s.write(f"min_seq_length={a.min_length}\n")
    s.write(f"short_sequences_removed={len(too_short)}\n")
    s.write(f"short_bases_removed={sum(len(new_seq[o]) for o in too_short)}\n")
