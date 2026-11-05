// Contaminant screening (gated behind --run_kraken2).
// Kraken2 classifies each assembled contig against the shared NCBI reference DB, then its
// per-contig taxids are turned into a BlobTools "hits" file so the coverage-vs-GC blob plot
// is coloured by taxonomy — non-Viridiplantae contigs stand out as candidate contamination.
//
// Uses PlusPFP (Standard + protozoa/fungi/plant), so host (sorghum) sequence gets a real match
// instead of being scored against a DB with no plant clade — see nextflow.config for sizing.
// PlusPFP's ~231.5 GB loaded index needs assemble_heavy's 448 GB / big-RAM-node restriction,
// not qc_heavy's 128 GB — see the withLabel: 'assemble_heavy' block in nextflow.config.

process KRAKEN2_CLASSIFY {
    tag           { "${sample_id}_${stage}" }
    label         'assemble_heavy'
    // Advisory on its own, but --flag_contaminants depends on it: ignoring a failure there would
    // exit 0 with no decontam genome, summary or package.
    errorStrategy { !params.flag_contaminants ? 'ignore' : task.exitStatus in [132, 137, 139, 143, 247] ? 'retry' : 'terminate' }
    publishDir    { "${params.outdir}/qc/kraken2/${sample_id}" }, mode: 'copy'
    // NB: biocontainers kraken2 images use a busybox base whose hardlinked applets
    // (linuxrc, usr/bin/[) fail to unpack under this cluster's rootless apptainer
    // ("unpriv.link: operation not permitted"). staphb/kraken2 is Ubuntu-based with
    // separate coreutils binaries (no problematic hardlinks) and unpacks cleanly.
    container     'quay.io/staphb/kraken2:2.1.3'

    input:
    tuple val(sample_id), val(stage), path(assembly), path(bam), path(bai)
    path kraken2_db

    output:
    tuple val(sample_id), val(stage), path(assembly), path(bam), path(bai),
          path("${sample_id}_${stage}.kraken2.hits"),                               emit: hits
    tuple val(sample_id), val(stage), path("${sample_id}_${stage}.kraken2.report"), emit: report

    script:
    def pfx = "${sample_id}_${stage}"
    """
    kraken2 --db ${kraken2_db} \\
        --threads ${task.cpus} \\
        --confidence ${params.kraken2_confidence} \\
        --output ${pfx}.kraken2.out \\
        --report ${pfx}.kraken2.report \\
        ${assembly}

    # BlobTools hits format: seqID <tab> taxID <tab> score. Kraken2 .out columns are
    # 1=C/U, 2=contigID, 3=taxID; keep classified contigs with a real taxid. Score is
    # unused (one hit per contig).
    awk -F'\\t' '\$1=="C" && \$3!="0" {print \$2"\\t"\$3"\\t1"}' ${pfx}.kraken2.out > ${pfx}.kraken2.hits
    """
}

process BLOBTOOLS_TAXONOMY {
    tag           { "${sample_id}_${stage}" }
    label         'qc'
    errorStrategy { !params.flag_contaminants ? 'ignore' : task.exitStatus in [132, 137, 139, 143, 247] ? 'retry' : 'terminate' }
    publishDir    { "${params.outdir}/qc/blobtools/${sample_id}" }, mode: 'copy'
    container     'quay.io/biocontainers/blobtools:1.1.1--py_1'

    input:
    tuple val(sample_id), val(stage), path(assembly), path(bam), path(bai), path(hits)
    path taxdump

    output:
    tuple val(sample_id), path("blobtax_${sample_id}_${stage}*"), emit: results

    script:
    def pfx = "blobtax_${sample_id}_${stage}"
    """
    # Build the BlobDB with taxonomy. --nodes/--names point at the taxdump so the Kraken2
    # taxids resolve to a lineage; default taxrule (bestsum) picks the per-contig hit.
    # --db gives a WRITABLE cwd path for the built nodesDB: without it blobtools tries to
    # cache nodesDB.txt back into its read-only package dir and dies (Errno 30). The file
    # does not exist yet, so blobtools builds it here from --nodes/--names.
    blobtools create \\
        -i ${assembly} \\
        -b ${bam} \\
        -t ${hits} \\
        --nodes ${taxdump}/nodes.dmp \\
        --names ${taxdump}/names.dmp \\
        --db nodesDB.txt \\
        -o ${pfx}

    # Blob plot + table coloured by phylum; prefix outputs so they match the publish glob.
    blobtools plot -i ${pfx}.blobDB.json --rank phylum --out ${pfx}
    blobtools view -i ${pfx}.blobDB.json --rank phylum -o ${pfx}
    """
}

// Contaminant flagging (gated behind --flag_contaminants, which itself requires --run_kraken2
// and --nuclear_ref). A contig is only flagged when ALL THREE hold:
//   (a) non-target-phylum by BlobTools ("no-hit" never counts — absence of a hit is not
//       evidence of contamination, only an explicit non-plant call is),
//   (b) unplaced by RagTag (derived structurally from the AGP, not from a chromosome-name
//       pattern — a contig is unplaced iff it is the sole component of its own AGP object),
//   (c) less than contam_max_ref_aligned of the contig aligns to the nuclear reference.
// Real chromosomes are always RagTag-placed (never flagged regardless of taxonomy noise), and
// plant contigs misclassified by Kraken2 are saved by aligning to the reference. GC% is recorded
// in the audit but is not a criterion: Burkholderia contigs near 68% GC passed a 0.20-0.70 GC
// rule and left thousands of bacterial gene models in the annotation of the "decontam" genome.
//
// Unplaced contigs with no taxonomy call (no-hit/undef) that also fail to align are never removed
// automatically. They may be divergent plant sequence, so they are listed for review instead and
// checked by the same advisory Sylph/BLASTn steps as the removed contigs.
process CLASSIFY_CONTAMINANTS {
    tag           { sample_id }
    label         'qc'
    // The decontam genome is the input to annotation, so a failure here must stop the run rather
    // than silently skip contaminant removal.
    errorStrategy 'terminate'
    publishDir    { "${params.outdir}/qc/contamination/${sample_id}" }, mode: 'copy'
    container     'quay.io/biocontainers/blobtools:1.1.1--py_1'

    input:
    // Pre-joined by sample_id in main.nf so the blob table, AGP, and reference alignment table
    // are guaranteed to be the same sample even with multiple samples.
    tuple val(sample_id), path(blob_files), path(agp), path(ref_aligned)

    output:
    tuple val(sample_id), path("${sample_id}_contaminant_ids.txt"),       emit: ids
    tuple val(sample_id), path("${sample_id}_review_ids.txt"),            emit: review_ids
    tuple val(sample_id), path("${sample_id}_contamination_audit.tsv"),  emit: audit
    tuple val(sample_id), path("${sample_id}_contamination_summary.txt"), emit: summary

    script:
    def target_phylum   = params.contam_target_phylum
    def max_ref_aligned = params.contam_max_ref_aligned
    def sample          = sample_id
    """
    #!/usr/bin/env bash
    set -uo pipefail

    BLOB_TABLE=\$(find . -maxdepth 1 -name '*.blobDB.table.txt' | head -1)
    if [ -z "\$BLOB_TABLE" ]; then
        echo "ERROR: no *.blobDB.table.txt found in BlobTools output" >&2
        exit 1
    fi

    python3 - "\$BLOB_TABLE" "${agp}" "${ref_aligned}" "${target_phylum}" "${max_ref_aligned}" "${sample}" <<'PYEOF'
import sys, collections

blob_table, agp_path, ref_aligned_path, target_phylum, max_ref_aligned, sample = sys.argv[1:7]
max_ref_aligned = float(max_ref_aligned)
no_taxonomy = {'no-hit', 'undef', 'N/A'}

# ---- 1. Parse BlobTools table: name -> (length, GC, phylum_call) ----
contigs = {}
with open(blob_table) as fh:
    for line in fh:
        if line.startswith('#') or not line.strip():
            continue
        f = line.rstrip('\\n').split('\\t')
        name, length, gc, phylum = f[0], int(f[1]), float(f[2]), f[5]
        contigs[name] = {'length': length, 'gc': gc, 'phylum': phylum}

# ---- 2. Parse AGP: derive placement purely from structure (no naming assumptions) ----
# Group rows by object (col1). A contig is unplaced iff it is the SOLE W-component of its
# own object (single-row object whose component_id == object name); anything joined into a
# multi-row scaffold object is placed, regardless of what that object is named.
by_object = collections.defaultdict(list)
with open(agp_path) as fh:
    for line in fh:
        if line.startswith('#') or not line.strip():
            continue
        f = line.rstrip('\\n').split('\\t')
        obj, comp_type, comp_id = f[0], f[4], f[5]
        by_object[obj].append((comp_type, comp_id))

unplaced = set()
for obj, rows in by_object.items():
    if len(rows) == 1 and rows[0][0] == 'W' and rows[0][1] == obj:
        unplaced.add(obj)

# ---- 3. Parse the unplaced-contig alignment table. Absent contigs had no alignment. ----
ref_aligned = {}
with open(ref_aligned_path) as fh:
    next(fh)
    for line in fh:
        f = line.rstrip('\\n').split('\\t')
        ref_aligned[f[0]] = float(f[3])

# ---- 4. Apply the rule, write audit + ID lists + summary ----
flagged_ids = []
review_ids = []
phyla_seen = collections.Counter()
bases_removed = 0
review_bases = 0

audit_path   = f"{sample}_contamination_audit.tsv"
ids_path     = f"{sample}_contaminant_ids.txt"
review_path  = f"{sample}_review_ids.txt"
summary_path = f"{sample}_contamination_summary.txt"

with open(audit_path, 'w') as audit, open(ids_path, 'w') as ids, open(review_path, 'w') as review:
    audit.write("name\\tlength\\tGC\\tphylum\\tplaced\\tflagged\\tref_aligned_frac\\tdecision\\n")
    for name, info in contigs.items():
        is_unplaced = name in unplaced
        frac        = ref_aligned.get(name, 0.0)
        unaligned   = is_unplaced and frac < max_ref_aligned
        has_call    = info['phylum'] not in no_taxonomy
        flagged     = unaligned and has_call and info['phylum'] != target_phylum
        needs_review = unaligned and not has_call
        decision    = 'remove' if flagged else 'review' if needs_review else 'keep'

        placed_str = 'unplaced' if is_unplaced else 'placed'
        frac_str   = f"{frac:.4f}" if is_unplaced else 'NA'
        audit.write(f"{name}\\t{info['length']}\\t{info['gc']}\\t{info['phylum']}\\t{placed_str}\\t{'yes' if flagged else 'no'}\\t{frac_str}\\t{decision}\\n")
        if flagged:
            ids.write(name + "\\n")
            flagged_ids.append(name)
            phyla_seen[info['phylum']] += 1
            bases_removed += info['length']
        elif needs_review:
            review.write(name + "\\n")
            review_ids.append(name)
            review_bases += info['length']

with open(summary_path, 'w') as s:
    s.write(f"contaminants_flagged={len(flagged_ids)}\\n")
    s.write(f"bases_removed={bases_removed}\\n")
    s.write("phyla_detected=" + ",".join(f"{p}:{c}" for p, c in phyla_seen.items()) + "\\n")
    s.write(f"review_contigs={len(review_ids)}\\n")
    s.write(f"review_bases={review_bases}\\n")
PYEOF
    """
}

// Measures how much of each unplaced contig aligns to the nuclear reference, for the flagging
// rule in CLASSIFY_CONTAMINANTS. asm20 tolerates divergence between sorghum lines; overlapping
// alignments are merged so repeats are not counted twice. Placement uses the same AGP structure
// rule as CLASSIFY_CONTAMINANTS.
process ALIGN_UNPLACED_TO_REFERENCE {
    tag        { sample_id }
    label      'map'
    publishDir { "${params.outdir}/qc/contamination/${sample_id}" }, mode: 'copy'
    container  'quay.io/biocontainers/minimap2:2.30--h577a1d6_0'

    input:
    tuple val(sample_id), path(assembly), path(agp)
    path reference

    output:
    tuple val(sample_id), path("${sample_id}_unplaced_ref_aligned.tsv"), emit: aligned

    script:
    """
    awk -F '\\t' '!/^#/ && NF { n[\$1]++; type[\$1]=\$5; comp[\$1]=\$6 }
         END { for (o in n) if (n[o] == 1 && type[o] == "W" && comp[o] == o) print o }' \\
        ${agp} > unplaced_ids.txt

    awk 'NR == FNR { want[\$1] = 1; next }
         /^>/ { name = substr(\$1, 2); keep = (name in want) }
         keep' unplaced_ids.txt ${assembly} > unplaced.fasta

    if [ -s unplaced.fasta ]; then
        minimap2 -x asm20 --secondary=no -t ${task.cpus} ${reference} unplaced.fasta > unplaced_to_ref.paf
    else
        : > unplaced_to_ref.paf
    fi

    sort -k1,1 -k3,3n unplaced_to_ref.paf | awk -F '\\t' -v OFS='\\t' '
        function flush() { if (q != "") print q, qlen, cov, cov / qlen }
        BEGIN { print "contig", "length", "aligned_bp", "aligned_frac" }
        \$1 != q { flush(); q = \$1; qlen = \$2; cov = 0; end = -1 }
        {
            if (\$3 > end) { cov += \$4 - \$3; end = \$4 }
            else if (\$4 > end) { cov += \$4 - end; end = \$4 }
        }
        END { flush() }' > ${sample_id}_unplaced_ref_aligned.tsv
    """
}

// Second-opinion corroboration for the (typically very few) contigs CLASSIFY_CONTAMINANTS
// flagged: an independent containment-ANI method (Sylph) against real GTDB reference genomes,
// rather than Kraken2's k-mer LCA approach. Only queries the flagged candidates — cheap even
// though the GTDB database itself is large — so this is advisory extra confidence, not a gate
// on removal (REMOVE_CONTAMINANTS' behavior is unchanged either way).
process EXTRACT_CANDIDATE_CONTIGS {
    tag           { sample_id }
    label         'qc'
    container     'quay.io/biocontainers/seqkit:2.13.0--he881be0_0'

    input:
    tuple val(sample_id), path(assembly), path(contaminant_ids), path(review_ids)

    output:
    tuple val(sample_id), path("${sample_id}_candidates.fasta"), emit: fasta

    script:
    """
    cat ${contaminant_ids} ${review_ids} > candidate_ids.txt
    if [ -s candidate_ids.txt ]; then
        seqkit grep -f candidate_ids.txt ${assembly} -o ${sample_id}_candidates.fasta
    else
        : > ${sample_id}_candidates.fasta
    fi
    """
}

process SYLPH_VERIFY_CONTAMINANTS {
    tag           { sample_id }
    label         'qc'
    // Empty candidate sets and no-match queries are handled successfully in
    // the script below. A genuine Sylph/runtime failure must remain visible
    // and prevent publication of a package without corroboration output.
    errorStrategy { task.exitStatus in [132, 137, 139, 143, 247] ? 'retry' : 'terminate' }
    publishDir    { "${params.outdir}/qc/contamination/${sample_id}/sylph" }, mode: 'copy'
    container     'quay.io/biocontainers/sylph:0.9.0--ha6fb395_0'

    input:
    tuple val(sample_id), path(candidates)
    path sylph_db

    output:
    tuple val(sample_id), path("${sample_id}_sylph_corroboration.tsv"), emit: corroboration
    tuple val(sample_id), path("${sample_id}_sylph_summary.txt"),       emit: summary

    script:
    def sample  = sample_id
    def min_ani = params.sylph_min_ani
    """
    #!/usr/bin/env bash
    set -euo pipefail

    if [ ! -s ${candidates} ]; then
        printf "contig\\tbest_ani\\tbest_match_genome\\tcorroborated\\n" > ${sample}_sylph_corroboration.tsv
        echo "sylph_corroborated=0" > ${sample}_sylph_summary.txt
        echo "sylph_total_candidates=0" >> ${sample}_sylph_summary.txt
        exit 0
    fi

    # Sketch each candidate contig as its OWN sample (-r: read-mode sketch, not a genome
    # database — 'query' mode's Contig_name column identifies the matched REFERENCE genome's
    # contig, not the query's, so pooling multiple query contigs into one sketch would make
    # results impossible to attribute back to a specific candidate). One .sylsp per contig
    # keeps them distinguishable via the Sample_file column in the query output.
    mkdir -p per_contig sketches
    awk '/^>/{if (f) close(f); f="per_contig/" substr(\$1,2) ".fasta"} {print > f}' ${candidates}
    for f in per_contig/*.fasta; do
        sylph sketch -r "\$f" -d sketches -t ${task.cpus}
    done

    sylph query ${sylph_db} sketches/*.sylsp \\
        -m ${min_ani} -t ${task.cpus} -o ${sample}_sylph_query.tsv

    # Parse the tab-delimited Sylph output with awk so this process does not
    # require Python in the Sylph container. Sample_file identifies the query
    # contig; Contig_name identifies the matched reference contig and must not
    # be used for attribution. Sylph's -m threshold means every retained row
    # is already a qualifying corroboration.
    awk -F '\\t' -v OFS='\\t' '
        FILENAME == ARGV[1] {
            if (FNR == 1) {
                for (i = 1; i <= NF; i++) {
                    if (\$i == "Sample_file") sample_col = i
                    if (\$i == "Genome_file") genome_col = i
                    if (\$i == "Adjusted_ANI") ani_col = i
                }
                next
            }
            if (!sample_col || !genome_col || !ani_col || \$sample_col == "" || \$ani_col == "") next

            contig = \$sample_col
            sub("^.*/", "", contig)
            sub(".fasta\$", "", contig)
            ani = \$ani_col + 0
            if (!(contig in best_ani) || ani > best_ani[contig]) {
                best_ani[contig] = ani
                best_genome[contig] = \$genome_col
            }
            next
        }

        FILENAME == ARGV[2] {
            if (\$0 ~ /^>/) {
                contig = substr(\$0, 2)
                sub(/[[:space:]].*\$/, "", contig)
                contigs[++n] = contig
            }
            next
        }

        END {
            print "contig", "best_ani", "best_match_genome", "corroborated"
            for (i = 1; i <= n; i++) {
                contig = contigs[i]
                if (contig in best_ani)
                    printf "%s\\t%.2f\\t%s\\tyes\\n", contig, best_ani[contig], best_genome[contig]
                else
                    print contig, "NA", "NA", "no"
            }
        }
    ' "${sample}_sylph_query.tsv" "${candidates}" > "${sample}_sylph_corroboration.tsv"

    corroborated=\$(awk -F '\\t' 'NR > 1 && \$4 == "yes" {n++} END {print n + 0}' "${sample}_sylph_corroboration.tsv")
    total_candidates=\$(awk -F '\\t' 'NR > 1 {n++} END {print n + 0}' "${sample}_sylph_corroboration.tsv")
    printf 'sylph_corroborated=%s\\n' "\${corroborated}" > "${sample}_sylph_summary.txt"
    printf 'sylph_total_candidates=%s\\n' "\${total_candidates}" >> "${sample}_sylph_summary.txt"
    """
}

// Advisory nucleotide-level corroboration for the same candidates checked by Sylph (removed
// contigs plus contigs listed for review). Only ever searches a local database: assembly
// sequence must not be sent to a remote BLAST service. BLAST evidence never changes removal.
process BLAST_VERIFY_CONTAMINANTS {
    tag           { sample_id }
    label         'qc_heavy'
    errorStrategy { task.exitStatus in [132, 137, 139, 143, 247] ? 'retry' : 'terminate' }
    publishDir    { "${params.outdir}/qc/contamination/${sample_id}/blast" }, mode: 'copy'
    container     'quay.io/biocontainers/blast:2.15.0--pl5321h6f7f691_1'

    input:
    tuple val(sample_id), path(candidates)
    path blast_db_dir
    val blast_db_name

    output:
    tuple val(sample_id), path("${sample_id}_blast_corroboration.tsv"), emit: corroboration
    tuple val(sample_id), path("${sample_id}_blast_summary.txt"),       emit: summary

    script:
    def sample = sample_id
    def min_id = params.blast_min_identity
    def min_qc = params.blast_min_qcov
    def max_t  = params.blast_max_targets
    """
    #!/usr/bin/env bash
    set -euo pipefail

    if [ ! -s ${candidates} ]; then
        printf 'contig\\tbest_identity\\tbest_qcov\\tbest_subject\\tbest_title\\tevalue\\tbitscore\\tcorroborated\\n' > ${sample}_blast_corroboration.tsv
        echo 'blast_corroborated=0' > ${sample}_blast_summary.txt
        echo 'blast_total_candidates=0' >> ${sample}_blast_summary.txt
        exit 0
    fi

    export BLASTDB="\$PWD/${blast_db_dir}"
    blastn -query ${candidates} -db "${blast_db_name}" \\
        -num_threads ${task.cpus} -max_target_seqs ${max_t} \\
        -evalue 1e-5 -outfmt '6 qseqid sacc stitle pident length qlen evalue bitscore qcovs' \\
        -out ${sample}_blast_query.tsv

    # Select the highest-bitscore hit per query. A hit corroborates the candidate only
    # when both identity and query coverage meet explicit thresholds. Coverage is qcovs (all
    # HSPs to that subject): a single HSP rarely spans a Mb contig, even from the right species.
    awk -F '\\t' -v OFS='\\t' -v minid=${min_id} -v minqc=${min_qc} \\
        'BEGIN { print "contig", "best_identity", "best_qcov", "best_subject", "best_title", "evalue", "bitscore", "corroborated" }
         FILENAME == ARGV[1] { q=\$1; qcov=\$9; if (!(q in best) || \$8 > best[q]) {
             best[q]=\$8; pid[q]=\$4; cov[q]=qcov; sid[q]=\$2; title[q]=\$3; ev[q]=\$7
         }}
         FILENAME == ARGV[2] && /^>/ { q=substr(\$1,2); sub(/[[:space:]].*\$/, "", q); contigs[++n]=q }
         END { for (i=1; i<=n; i++) { q=contigs[i]; if (q in best)
                    printf "%s\\t%.2f\\t%.2f\\t%s\\t%s\\t%s\\t%s\\t%s\\n", q,pid[q],cov[q],sid[q],title[q],ev[q],best[q],(pid[q]>=minid && cov[q]>=minqc ? "yes" : "no")
                else print q, "NA", "NA", "NA", "NA", "NA", "NA", "no" } }' \\
        ${sample}_blast_query.tsv ${candidates} >> ${sample}_blast_corroboration.tsv

    total=\$(awk 'NR > 1 {n++} END {print n+0}' ${sample}_blast_corroboration.tsv)
    corroborated=\$(awk -F '\\t' 'NR > 1 && \$8 == "yes" {n++} END {print n+0}' ${sample}_blast_corroboration.tsv)
    printf 'blast_corroborated=%s\\n' "\${corroborated}" > ${sample}_blast_summary.txt
    printf 'blast_total_candidates=%s\\n' "\${total}" >> ${sample}_blast_summary.txt
    """
}

// Splits the final genome into the decontam FASTA and the removed contigs, so nothing is lost and
// any removal can be checked or reverted. Deliberately NOT errorStrategy 'ignore': the decontam
// FASTA is a hard dependency for the QUAST_NUCLEAR/BUSCO_NUCLEAR 'decontam' columns downstream.
process REMOVE_CONTAMINANTS {
    tag        { sample_id }
    label      'qc'
    publishDir { "${params.outdir}/assembly/decontam/${sample_id}" }, mode: 'copy'
    container  'quay.io/biocontainers/seqkit:2.13.0--he881be0_0'

    input:
    tuple val(sample_id), path(assembly), path(contaminant_ids)

    output:
    tuple val(sample_id), path("${sample_id}_decontam.fasta"),        emit: assembly
    tuple val(sample_id), path("${sample_id}_removed_contigs.fasta"), emit: removed

    script:
    """
    if [ -s ${contaminant_ids} ]; then
        seqkit grep -v -f ${contaminant_ids} ${assembly} -o ${sample_id}_decontam.fasta
        seqkit grep -f ${contaminant_ids} ${assembly} -o ${sample_id}_removed_contigs.fasta
    else
        cp ${assembly} ${sample_id}_decontam.fasta
        : > ${sample_id}_removed_contigs.fasta
    fi

    # The two outputs must account for every input contig and base, and hold every listed ID.
    count() { seqkit fx2tab -n -l "\$1" | awk -F '\\t' '{ n++; b += \$2 } END { printf "%.0f %.0f\\n", n + 0, b + 0 }'; }
    read total_n total_b < <(count ${assembly})
    read kept_n kept_b < <(count ${sample_id}_decontam.fasta)
    read removed_n removed_b < <(count ${sample_id}_removed_contigs.fasta)
    listed_n=\$(grep -c . ${contaminant_ids} || true)
    if [ "\$total_n" != "\$(awk -v a="\$kept_n" -v b="\$removed_n" 'BEGIN { printf "%.0f", a + b }')" ] \\
        || [ "\$total_b" != "\$(awk -v a="\$kept_b" -v b="\$removed_b" 'BEGIN { printf "%.0f", a + b }')" ] \\
        || [ "\$removed_n" != "\$listed_n" ]; then
        echo "ERROR: split mismatch: input \$total_n/\$total_b, kept \$kept_n/\$kept_b, removed \$removed_n/\$removed_b, listed \$listed_n" >&2
        exit 1
    fi
    """
}

// NCBI FCS-GX: genome cross-species screen, run on the full genome ('full', input to the removal
// rule) and again on the published genome ('final', which must come back clean). The ~470 GiB DB
// is copied to /dev/shm because memory-mapping it over cephfs is not viable. FCS rejects '+' in
// IDs, so it is replaced with 'plus' here and restored by canopy_contig_rule.py.
process FCS_GX_SCREEN {
    tag              { "${sample_id}_${stage}" }
    label            'fcs_gx'
    errorStrategy    { task.exitStatus in [132, 137, 139, 143, 247] ? 'retry' : 'terminate' }
    publishDir       { "${params.outdir}/qc/contamination/${sample_id}/fcs_gx" }, mode: 'copy'
    container        'quay.io/biocontainers/ncbi-fcs-gx:0.4.1--h4ac6f70_0'
    containerOptions '--bind /dev/shm'

    input:
    tuple val(sample_id), val(stage), path(assembly)
    path gx_db

    output:
    tuple val(sample_id), val(stage), path("${sample_id}_${stage}.fcs_gx_report.txt"), emit: report
    path "${sample_id}_${stage}.taxonomy.rpt", emit: taxonomy

    script:
    def base = "${sample_id}_${stage}"
    """
    SHM=/dev/shm/gxdb_${base}_\$\$
    trap 'rm -rf "\$SHM"' EXIT
    trap 'exit 143' INT TERM
    # The container's df is BusyBox, which has no --output.
    avail=\$(df -k /dev/shm | tail -1 | awk '{ print int(\$(NF-2) / 1048576) }')
    if [ "\${avail:-0}" -lt 480 ]; then
        echo "ERROR: /dev/shm has \${avail:-0}G free; FCS-GX needs 480G" >&2
        exit 1
    fi
    mkdir -p "\$SHM"
    cp ${gx_db}/all.* "\$SHM"/
    sed '/^>/ { s/ .*//; s/+/plus/g }' ${assembly} > ${base}.input.fasta
    GX_NUM_CORES=${task.cpus} run_gx.py --fasta ${base}.input.fasta \\
        --tax-id ${params.fcs_gx_taxid} --gx-db "\$SHM/all" --out-dir . --out-basename ${base}
    rm -f ${base}.input.fasta
    """
}

// NCBI FCS-adaptor, eukaryote mode. The image has no bash, which Nextflow's task wrapper needs,
// so the task runs on the node and calls the container itself.
process FCS_ADAPTOR_SCREEN {
    tag           { "${sample_id}_${stage}" }
    label         'qc'
    errorStrategy { task.exitStatus in [132, 137, 139, 143, 247] ? 'retry' : 'terminate' }
    publishDir    { "${params.outdir}/qc/contamination/${sample_id}/fcs_adaptor" }, mode: 'copy'

    input:
    tuple val(sample_id), val(stage), path(assembly)

    output:
    tuple val(sample_id), val(stage), path("${sample_id}_${stage}.fcs_adaptor_report.txt"), emit: report

    script:
    def base = "${sample_id}_${stage}"
    """
    export APPTAINER_CACHEDIR="${System.getenv('HOME')}/.apptainer_cache"
    sed '/^>/ { s/ .*//; s/+/plus/g }' ${assembly} > ${base}.input.fasta
    apptainer exec --bind "\$PWD" docker://${params.fcs_adaptor_image} \\
        /app/fcs/bin/av_screen_x -o "\$PWD/adaptor_out" --euk "\$PWD/${base}.input.fasta"
    cp adaptor_out/fcs_adaptor_report.txt ${base}.fcs_adaptor_report.txt
    rm -rf ${base}.input.fasta adaptor_out
    """
}

// Fungal alignment for the removal rule. Kraken's PlusPFP has no Rhizopus, so whole fungal
// genomes supply evidence it cannot give. Alignment alone is not proof: conserved high-copy plant
// contigs also align, which is why the rule also requires little alignment to the sorghum reference.
process ALIGN_UNPLACED_TO_FUNGI {
    tag        { sample_id }
    label      'map'
    publishDir { "${params.outdir}/qc/contamination/${sample_id}" }, mode: 'copy'
    container  'quay.io/biocontainers/minimap2:2.30--h577a1d6_0'
    containerOptions '--unsquash'

    input:
    tuple val(sample_id), path(assembly), path(agp)
    path fungal_refs

    output:
    tuple val(sample_id), path("${sample_id}_unplaced_fungal_aligned.tsv"), emit: aligned

    script:
    """
    awk -F '\\t' '!/^#/ && NF { n[\$1]++; type[\$1]=\$5; comp[\$1]=\$6 }
         END { for (o in n) if (n[o] == 1 && type[o] == "W" && comp[o] == o) print o }' \\
        ${agp} > unplaced_ids.txt

    awk 'NR == FNR { want[\$1] = 1; next }
         /^>/ { name = substr(\$1, 2); keep = (name in want) }
         keep' unplaced_ids.txt ${assembly} > unplaced.fasta

    if [ -s unplaced.fasta ]; then
        minimap2 -x asm20 --secondary=no -t ${task.cpus} ${fungal_refs} unplaced.fasta > fungal.paf
    else
        : > fungal.paf
    fi

    sort -k1,1 -k3,3n fungal.paf | awk -F '\\t' -v OFS='\\t' '
        function flush() { if (q != "") print q, l, c, c / l, sp }
        BEGIN { print "contig", "length", "aligned_bp", "aligned_frac", "species" }
        \$1 != q { flush(); q = \$1; l = \$2; c = 0; e = -1; split(\$6, a, "|"); sp = a[1] }
        {
            if (\$3 > e) { c += \$4 - \$3; e = \$4 }
            else if (\$4 > e) { c += \$4 - e; e = \$4 }
        }
        END { flush() }' > ${sample_id}_unplaced_fungal_aligned.tsv
    rm -f unplaced.fasta
    """
}

// Combines all evidence into one decision per sequence; see bin/canopy_contig_rule.py for the rule.
process APPLY_CONTAMINANT_RULE {
    tag           { sample_id }
    label         'qc_light'
    errorStrategy 'terminate'
    publishDir    { "${params.outdir}/qc/contamination/${sample_id}/decisions" }, mode: 'copy'
    container     'quay.io/biocontainers/blobtools:1.1.1--py_1'

    input:
    tuple val(sample_id), path(audit), path(blob_files), path(fcs_report), path(fungal), path(blast)
    path rule_script

    output:
    tuple val(sample_id), path("${sample_id}_contig_decisions.tsv"),      emit: decisions
    tuple val(sample_id), path("${sample_id}_contaminant_ids.txt"),       emit: ids
    tuple val(sample_id), path("${sample_id}_review_ids.txt"),            emit: review_ids
    tuple val(sample_id), path("${sample_id}_contamination_summary.txt"), emit: summary

    script:
    """
    BLOB_TABLE=\$(find -L . -maxdepth 1 -name 'blobtax_*.blobDB.table.txt' | head -1)
    python3 ${rule_script} --sample ${sample_id} --audit ${audit} --blob-table "\$BLOB_TABLE" \\
        --fcs-report ${fcs_report} --fungal ${fungal} --blast ${blast} \\
        --min-fungal ${params.contam_min_fungal_frac} --max-ref ${params.contam_max_ref_aligned}
    """
}

// Writes <sample>_final and <sample>_full_flagged; see bin/canopy_write_genomes.py.
process WRITE_GENOME_VERSIONS {
    tag           { sample_id }
    label         'qc_light'
    errorStrategy 'terminate'
    publishDir    { "${params.outdir}/assembly/final/${sample_id}" }, mode: 'copy'
    container     'quay.io/biocontainers/blobtools:1.1.1--py_1'

    input:
    tuple val(sample_id), path(assembly), path(agp), path(decisions), val(stage), path(adaptor_report)
    path writer_script

    output:
    tuple val(sample_id), path("${sample_id}_final.fasta"), emit: assembly
    tuple val(sample_id), path("${sample_id}_final.agp"), path("${sample_id}_full_flagged.fasta"),
          path("${sample_id}_full_flagged.agp"), path("${sample_id}_adapter_trims.tsv"), emit: extras
    tuple val(sample_id), path("${sample_id}_final_stats.txt"), emit: stats

    script:
    """
    python3 ${writer_script} --sample ${sample_id} --fasta ${assembly} --agp ${agp} \\
        --decisions ${decisions} --adaptor-report ${adaptor_report}
    """
}

// Second round for a genome whose final screen had hits: the same rule sees both FCS-GX reports,
// then the first-round genome is trimmed again. Its files are staged in round1/ because the
// outputs have the same names. Sequences stay in both rounds' FASTAs under 'plus' IDs, so the
// decision table is converted to match them for the writer.
process REFINE_GENOME_VERSIONS {
    tag           { sample_id }
    label         'qc_light'
    errorStrategy 'terminate'
    publishDir    { "${params.outdir}/qc/contamination/${sample_id}/decisions" }, mode: 'copy', overwrite: true,
                  pattern: "*_{contig_decisions.tsv,contaminant_ids.txt,review_ids.txt,contamination_summary.txt}"
    publishDir    { "${params.outdir}/assembly/final/${sample_id}" }, mode: 'copy', overwrite: true,
                  pattern: "*_{final.fasta,final.agp,full_flagged.fasta,full_flagged.agp,adapter_trims.tsv,final_stats.txt}"
    container     'quay.io/biocontainers/blobtools:1.1.1--py_1'

    input:
    tuple val(sample_id), path(audit), path(blob_files), path(fcs_full), path(fungal), path(blast),
          path(flagged, stageAs: 'round1/*'), path(flagged_agp, stageAs: 'round1/*'),
          path(prior_trims, stageAs: 'round1/*'), path(gx_final), path(adaptor_final)
    path rule_script
    path writer_script
    path merge_script

    output:
    tuple val(sample_id), path("${sample_id}_final.fasta"), emit: assembly
    tuple val(sample_id), path("${sample_id}_final.agp"), path("${sample_id}_full_flagged.fasta"),
          path("${sample_id}_full_flagged.agp"), path("${sample_id}_adapter_trims.tsv"), emit: extras
    tuple val(sample_id), path("${sample_id}_final_stats.txt"), emit: stats
    tuple val(sample_id), path("${sample_id}_contig_decisions.tsv"), emit: decisions
    tuple val(sample_id), path("${sample_id}_contamination_summary.txt"), emit: summary
    path "${sample_id}_contaminant_ids.txt"
    path "${sample_id}_review_ids.txt"

    script:
    """
    BLOB_TABLE=\$(find -L . -maxdepth 1 -name 'blobtax_*.blobDB.table.txt' | head -1)
    cat ${fcs_full} ${gx_final} > gx_both_rounds.txt
    python3 ${rule_script} --sample ${sample_id} --audit ${audit} --blob-table "\$BLOB_TABLE" \\
        --fcs-report gx_both_rounds.txt --fungal ${fungal} --blast ${blast} \\
        --min-fungal ${params.contam_min_fungal_frac} --max-ref ${params.contam_max_ref_aligned}
    awk -F '\\t' -v OFS='\\t' 'NR > 1 { gsub(/\\+/, "plus", \$1) } { print }' \\
        ${sample_id}_contig_decisions.tsv > decisions_safe_ids.tsv
    python3 ${writer_script} --sample ${sample_id} --fasta ${flagged} --agp ${flagged_agp} \\
        --decisions decisions_safe_ids.tsv --adaptor-report ${adaptor_final}
    python3 ${merge_script} --prior ${prior_trims} --new ${sample_id}_adapter_trims.tsv \\
        --decisions ${sample_id}_contig_decisions.tsv --stats ${sample_id}_final_stats.txt
    rm -f gx_both_rounds.txt decisions_safe_ids.tsv
    """
}

// Screens of the published FASTA, as one key=value file for the report. Hits are reported rather
// than failing the run, so the evidence is published for inspection.
process NCBI_SCREEN_SUMMARY {
    tag        { sample_id }
    label      'qc_light'
    publishDir { "${params.outdir}/qc/contamination/${sample_id}" }, mode: 'copy'
    container  'quay.io/biocontainers/seqkit:2.13.0--he881be0_0'

    input:
    tuple val(sample_id), path(stats), path(gx_report), path(adaptor_report)

    output:
    tuple val(sample_id), path("${sample_id}_ncbi_screen.txt"), emit: summary

    script:
    """
    gx_exclude=\$(awk -F '\\t' '!/^#/ && \$5 == "EXCLUDE"' ${gx_report} | wc -l)
    gx_other=\$(awk -F '\\t' '!/^#/ && NF && \$5 != "EXCLUDE"' ${gx_report} | wc -l)
    adaptor=\$(awk '!/^#/ && NF' ${adaptor_report} | wc -l)
    db=\$(head -1 ${gx_report} | grep -o '"build-date": "[0-9-]*"' | grep -o '[0-9][0-9-]*' | head -1)
    status=PASS
    if [ "\$gx_exclude" != 0 ] || [ "\$gx_other" != 0 ] || [ "\$adaptor" != 0 ]; then status=CHECK; fi
    {
        cat ${stats}
        echo "gx_exclude=\$gx_exclude"
        echo "gx_other_actions=\$gx_other"
        echo "adaptor_hits=\$adaptor"
        echo "gx_db_build=\${db:-unknown}"
        echo "screen_status=\$status"
    } > ${sample_id}_ncbi_screen.txt
    """
}
