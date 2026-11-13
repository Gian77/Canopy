#!/usr/bin/env python3
"""Give every sequence in a Canopy release a sample-prefixed name. Only names change, never bases.

Run once on a release built by canopy_build_release.py. For each sample:

  ChrNN_RagTag        -> <sample>_ChrNN
  any other nuclear   -> <sample>_<original name>        (unplaced contigs keep their contig ID)
  chloroplast         -> <sample>_chloroplast            (<sample>_chloroplast_N if in several pieces)
  mitochondrion       -> <sample>_mito_1 .. _mito_N      (longest first)

Renamed: the header names of the four FASTAs (the rest of each header is kept), the object names
(column 1) of both AGPs, and the sequence-name columns of the evidence tables. AGP component names
(column 6) are contig IDs and stay as they are. Reports from external tools (FCS, QUAST, BUSCO,
BlobTools) keep the original names; <sample>/genome/<sample>_sequence_names.tsv maps old to new.

Nuclear records keep their order; organelle records are written in the order of their new numbers.
Each FASTA is checked after rewriting: every record's sequence lines must be byte-identical to the
original record they came from, or the script stops before replacing anything.
"""
import argparse
import csv
import hashlib
import os
import re
import sys

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--release", required=True, help="release directory")
ap.add_argument("--samples", nargs="+", required=True)
a = ap.parse_args()


def records(path):
    """Yield (name, header line, list of sequence lines) with lines kept byte-exact."""
    name = header = None
    seq = []
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if header is not None:
                    yield name, header, seq
                header, seq = line, []
                name = line[1:].split()[0].decode()
            else:
                seq.append(line)
    if header is not None:
        yield name, header, seq


def seq_digests(path):
    return [(len(b"".join(l.strip() for l in seq)), hashlib.sha256(b"".join(seq)).hexdigest())
            for _, _, seq in records(path)]


def replace_atomically(path, writer):
    tmp = path + ".rename.tmp"
    with open(tmp, "wb") as out:
        writer(out)
    return tmp


def rename_fasta(path, mapping, order=None):
    """Write the renamed FASTA next to the original; return the temp path after verifying bases.

    order: optional list of original names giving the record order (organelles are written in the
    order of their new numbers); by default the original order is kept.
    """
    before = {n: d for (n, _, _), d in zip(records(path), seq_digests(path))}

    def write(out):
        recs = {name: (header, seq) for name, header, seq in records(path)}
        for name in order or list(recs):
            header, seq = recs[name]
            rest = header[1 + len(name.encode()):]
            out.write(b">" + mapping[name].encode() + rest)
            out.writelines(seq)

    tmp = replace_atomically(path, write)
    after = {n: d for (n, _, _), d in zip(records(tmp), seq_digests(tmp))}
    if len(after) != len(before) or any(after.get(mapping[n]) != d for n, d in before.items()):
        os.remove(tmp)
        sys.exit(f"ERROR: {path}: sequence bytes changed during renaming; nothing replaced")
    return tmp


def rename_agp(path, mapping):
    def write(out):
        with open(path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    out.write(line.encode())
                    continue
                f = line.rstrip("\n").split("\t")
                if f[0] not in mapping:
                    sys.exit(f"ERROR: {path}: object {f[0]} is not in the FASTA")
                f[0] = mapping[f[0]]
                out.write(("\t".join(f) + "\n").encode())
    return replace_atomically(path, write)


def rename_table(path, columns, mapping):
    """Rename the given columns; '+' in a name is read as 'plus', as the FASTA writes it."""
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    head = rows[0]
    idx = [head.index(c) for c in columns]
    for r in rows[1:]:
        for i in idx:
            key = r[i].replace("+", "plus")
            if key not in mapping:
                sys.exit(f"ERROR: {path}: {r[i]} is not a sequence in the genome")
            r[i] = mapping[key]

    def write(out):
        text = "".join("\t".join(r) + "\n" for r in rows)
        out.write(text.encode())
    return replace_atomically(path, write)


for s in a.samples:
    G = os.path.join(a.release, s, "genome")
    E = os.path.join(a.release, s, "evidence")
    fa = {k: os.path.join(G, f"{s}_{k}.fasta") for k in ("final", "full_flagged", "chloroplast", "mitochondria")}

    full = [(n, len(b"".join(l.strip() for l in seq))) for n, _, seq in records(fa["full_flagged"])]
    final = {n for n, _, _ in records(fa["final"])}
    if any(n.startswith(s + "_") for n, _ in full):
        sys.exit(f"ERROR: {s}: names already carry the sample prefix; the release was renamed before")
    if not final <= {n for n, _ in full}:
        sys.exit(f"ERROR: {s}: the final genome has sequences missing from the flagged genome")

    nuclear = {}
    for n, _ in full:
        m = re.fullmatch(r"(Chr\d+)_RagTag", n)
        nuclear[n] = f"{s}_{m.group(1)}" if m else f"{s}_{n}"
    cp = sorted(((n, len(b"".join(l.strip() for l in q))) for n, _, q in records(fa["chloroplast"])),
                key=lambda x: -x[1])
    mt = sorted(((n, len(b"".join(l.strip() for l in q))) for n, _, q in records(fa["mitochondria"])),
                key=lambda x: -x[1])
    cp_map = {n: f"{s}_chloroplast" if len(cp) == 1 else f"{s}_chloroplast_{i}" for i, (n, _) in enumerate(cp, 1)}
    mt_map = {n: f"{s}_mito_{i}" for i, (n, _) in enumerate(mt, 1)}

    new = list(nuclear.values()) + list(cp_map.values()) + list(mt_map.values())
    if len(new) != len(set(new)):
        sys.exit(f"ERROR: {s}: renaming would create duplicate names")
    if max(map(len, new)) > 50:
        sys.exit(f"ERROR: {s}: a new name exceeds NCBI's 50-character limit")

    # Build every replacement first; replace the originals only when all of them succeeded.
    pending = [
        (fa["final"], rename_fasta(fa["final"], nuclear)),
        (fa["full_flagged"], rename_fasta(fa["full_flagged"], nuclear)),
        (fa["chloroplast"], rename_fasta(fa["chloroplast"], cp_map, [n for n, _ in cp])),
        (fa["mitochondria"], rename_fasta(fa["mitochondria"], mt_map, [n for n, _ in mt])),
        (os.path.join(G, f"{s}_final.agp"), rename_agp(os.path.join(G, f"{s}_final.agp"), nuclear)),
        (os.path.join(G, f"{s}_full_flagged.agp"), rename_agp(os.path.join(G, f"{s}_full_flagged.agp"), nuclear)),
    ]
    tables = {
        f"{s}_contig_decisions.tsv": ["name"],
        f"{s}_adapter_trims.tsv": ["sequence"],
        f"{s}_placed_fcs_cuts.tsv": ["sequence"],
        f"{s}_blast_corroboration.tsv": None,
        f"{s}_sylph_corroboration.tsv": None,
    }
    for t, cols in tables.items():
        p = os.path.join(E, t)
        if os.path.exists(p):
            if cols is None:
                with open(p) as fh:
                    cols = [fh.readline().rstrip("\n").split("\t")[0]]
            pending.append((p, rename_table(p, cols, nuclear)))

    screen = os.path.join(E, "screening", f"{s}_ncbi_screen.txt")
    with open(screen) as fh:
        text = fh.read()
    fixed = re.sub(r"^(placed_fix_detail=)(\S+?):", lambda m: m.group(1) + nuclear[m.group(2)] + ":", text, flags=re.M)
    if fixed != text:
        pending.append((screen, replace_atomically(screen, lambda out: out.write(fixed.encode()))))

    names = os.path.join(G, f"{s}_sequence_names.tsv")
    with open(names + ".rename.tmp", "w") as out:
        out.write("new_name\toriginal_name\tfile\tlength_bp\tin_final\n")
        for n, length in full:
            out.write(f"{nuclear[n]}\t{n}\tfull_flagged\t{length}\t{'yes' if n in final else 'no'}\n")
        for n, length in cp:
            out.write(f"{cp_map[n]}\t{n}\tchloroplast\t{length}\tyes\n")
        for n, length in mt:
            out.write(f"{mt_map[n]}\t{n}\tmitochondria\t{length}\tyes\n")
    pending.append((names, names + ".rename.tmp"))

    for dst, tmp in pending:
        os.replace(tmp, dst)
    print(f"{s}: renamed {len(full)} nuclear, {len(cp)} chloroplast, {len(mt)} mitochondrial sequences; "
          f"{len(pending) - 1} files rewritten, bases verified identical")

# Files changed, so rewrite the checksums the same way canopy_build_release.py does.
manifest = os.path.join(a.release, "MANIFEST.sha256")
if os.path.exists(manifest):
    def sha256(path):
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for block in iter(lambda: fh.read(1 << 20), b""):
                h.update(block)
        return h.hexdigest()
    with open(manifest + ".rename.tmp", "w") as fh:
        for d, _, files in sorted(os.walk(a.release)):
            for f in sorted(files):
                p = os.path.join(d, f)
                if p != manifest and not p.endswith(".rename.tmp"):
                    fh.write(f"{sha256(p)}  {os.path.relpath(p, a.release)}\n")
    os.replace(manifest + ".rename.tmp", manifest)
    print(f"rewrote {manifest}")
