#!/usr/bin/env python3
"""Keep one organelle assembly out of the candidates from each assembler.

An assembly qualifies when it is one sequence, its assembler reports it as circular, and its length
is within --tolerance of the reference. The qualifying assembly closest to the reference length is
kept. If none qualifies, the Flye assembly is kept whole. Without Flye output, the candidate with
the most bases is kept.

Oatk and hifiasm mark a circular sequence by ending its name in "c". Flye reports it in the circ.
column of assembly_info.txt.
"""
import argparse

DISPLAY = {"oatk": "Oatk", "hifiasm": "hifiasm", "flye": "Flye"}
ORDER = ["oatk", "hifiasm", "flye"]  # also breaks ties in length


def read_fasta(path):
    records = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                records.append([line[1:].split()[0], []])
            elif line and records:
                records[-1][1].append(line)
    return [(name, "".join(parts)) for name, parts in records]


def read_flye_circular(path):
    """Names of the contigs Flye reports as circular."""
    circular = set()
    if not path:
        return circular
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if not line.startswith("#") and len(fields) > 3 and fields[3] == "Y":
                circular.add(fields[0])
    return circular


def describe(arm, records, flye_circular, ref_len, tolerance):
    names = [name for name, _ in records]
    circular = [n for n in names if (n in flye_circular if arm == "flye" else n.endswith("c"))]
    lengths = [len(seq) for _, seq in records]
    total = sum(lengths)
    diff = abs(total - ref_len) / ref_len if ref_len else float("inf")

    if not records:
        why = "no sequence"
    elif len(records) > 1:
        why = f"{len(records)} sequences"
    elif not circular:
        why = "not circular"
    elif diff > tolerance:
        why = f"length {diff:.1%} from reference"
    else:
        why = ""
    return {
        "arm": arm,
        "records": records,
        "n": len(records),
        "total": total,
        "largest": max(lengths, default=0),
        "circular": len(circular),
        "diff": diff,
        "qualifies": not why,
        "why": why,
    }


def choose(rows):
    """The chosen row, or None, and the sentence that explains the choice."""
    rank = {arm: i for i, arm in enumerate(ORDER)}
    qualified = [r for r in rows if r["qualifies"]]
    if qualified:
        best = min(qualified, key=lambda r: (r["diff"], rank.get(r["arm"], 99)))
        return best, (f"{DISPLAY[best['arm']]} qualifies as one circular sequence, "
                      f"{best['diff']:.2%} from the reference length.")
    flye = next((r for r in rows if r["arm"] == "flye" and r["n"]), None)
    if flye:
        return flye, "No assembly qualifies, so the Flye assembly is kept whole."
    assembled = [r for r in rows if r["n"]]
    if assembled:
        best = max(assembled, key=lambda r: r["total"])
        return best, (f"No assembly qualifies and Flye produced none, so {DISPLAY[best['arm']]}, "
                      f"the candidate with the most bases, is kept.")
    return None, "No assembler produced a sequence that passed the reference filter."


def write_fasta(path, records):
    with open(path, "w") as out:
        for name, seq in records:
            out.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i + 80] + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--compartment", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--tolerance", type=float, required=True)
    ap.add_argument("--candidate", action="append", required=True, metavar="ARM=FASTA")
    ap.add_argument("--flye-info")
    ap.add_argument("--output", required=True)
    ap.add_argument("--table", required=True)
    ap.add_argument("--markdown", required=True)
    args = ap.parse_args()

    ref_len = sum(len(seq) for _, seq in read_fasta(args.reference))
    flye_circular = read_flye_circular(args.flye_info)
    candidates = dict(item.split("=", 1) for item in args.candidate)
    rows = [describe(arm, read_fasta(candidates[arm]), flye_circular, ref_len, args.tolerance)
            for arm in ORDER if arm in candidates]
    chosen, reason = choose(rows)

    write_fasta(args.output, chosen["records"] if chosen else [])

    with open(args.table, "w") as out:
        out.write("sample\tcompartment\tassembler\tsequences\ttotal_bp\tlargest_bp\tcircular\t"
                  "fraction_from_reference\tqualifies\tchosen\tnot_qualified_because\n")
        for r in rows:
            out.write(f"{args.sample}\t{args.compartment}\t{r['arm']}\t{r['n']}\t{r['total']}\t"
                      f"{r['largest']}\t{r['circular']}\t{r['diff']:.4f}\t"
                      f"{'yes' if r['qualifies'] else 'no'}\t{'yes' if r is chosen else 'no'}\t"
                      f"{r['why'] or 'NA'}\n")

    with open(args.markdown, "w") as out:
        out.write(f"Reference length: {ref_len:,} bp. Allowed difference: {args.tolerance:.0%}.\n\n")
        out.write("| Assembler | Sequences | Total bp | Largest bp | Circular | From reference "
                  "| Qualifies | Kept |\n|---|---|---|---|---|---|---|---|\n")
        for r in rows:
            circular = f"{r['circular']} of {r['n']}" if r["n"] else "-"
            size = f"{r['diff']:.1%}" if r["n"] else "-"
            verdict = "yes" if r["qualifies"] else f"no ({r['why']})"
            out.write(f"| {DISPLAY[r['arm']]} | {r['n']} | {r['total']:,} | {r['largest']:,} | {circular} "
                      f"| {size} | {verdict} | {'**yes**' if r is chosen else ''} |\n")
        out.write(f"\n**Chosen:** {reason}\n")


if __name__ == "__main__":
    main()
