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
    errorStrategy 'ignore'     // advisory screening — never fail the pipeline
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
    errorStrategy 'ignore'
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
//   (c) GC% outside the normal eukaryotic range.
// This 3-way AND is a deliberate safety net: real chromosomes are always RagTag-placed (never
// flagged regardless of taxonomy noise), and real plant contigs occasionally misclassified by
// Kraken2 are saved by normal GC%.
process CLASSIFY_CONTAMINANTS {
    tag           { sample_id }
    label         'qc'
    errorStrategy 'ignore'   // advisory classification — never fail the pipeline
    publishDir    { "${params.outdir}/qc/contamination/${sample_id}" }, mode: 'copy'
    container     'quay.io/biocontainers/blobtools:1.1.1--py_1'

    input:
    // Pre-joined by sample_id in main.nf (BLOBTOOLS_TAXONOMY.out.results.join(final_agp_ch))
    // so blob table and AGP are guaranteed to be the same sample even with multiple samples.
    tuple val(sample_id), path(blob_files), path(agp)

    output:
    tuple val(sample_id), path("${sample_id}_contaminant_ids.txt"),       emit: ids
    tuple val(sample_id), path("${sample_id}_contamination_audit.tsv"),  emit: audit
    tuple val(sample_id), path("${sample_id}_contamination_summary.txt"), emit: summary

    script:
    def target_phylum = params.contam_target_phylum
    def gc_min         = params.contam_gc_min
    def gc_max         = params.contam_gc_max
    def sample         = sample_id
    """
    #!/usr/bin/env bash
    set -uo pipefail

    BLOB_TABLE=\$(find . -maxdepth 1 -name '*.blobDB.table.txt' | head -1)
    if [ -z "\$BLOB_TABLE" ]; then
        echo "WARNING: no *.blobDB.table.txt found in BlobTools output — skipping classification" >&2
        : > ${sample}_contaminant_ids.txt
        printf "name\\tlength\\tGC\\tphylum\\tplaced\\tflagged\\n" > ${sample}_contamination_audit.tsv
        echo "contaminants_flagged=0 (blob table missing)" > ${sample}_contamination_summary.txt
        exit 0
    fi

    python3 - "\$BLOB_TABLE" "${agp}" "${target_phylum}" "${gc_min}" "${gc_max}" "${sample}" <<'PYEOF'
import sys, collections

blob_table, agp_path, target_phylum, gc_min, gc_max, sample = sys.argv[1:7]
gc_min, gc_max = float(gc_min), float(gc_max)

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

# ---- 3. Apply the 3-way AND rule, write audit + ID list + summary ----
flagged_ids = []
phyla_seen = collections.Counter()
bases_removed = 0

audit_path   = f"{sample}_contamination_audit.tsv"
ids_path     = f"{sample}_contaminant_ids.txt"
summary_path = f"{sample}_contamination_summary.txt"

with open(audit_path, 'w') as audit, open(ids_path, 'w') as ids:
    audit.write("name\\tlength\\tGC\\tphylum\\tplaced\\tflagged\\n")
    for name, info in contigs.items():
        is_unplaced   = name in unplaced
        is_nonplant   = (info['phylum'] != target_phylum) and (info['phylum'] != 'no-hit')
        is_abnormalgc = not (gc_min <= info['gc'] <= gc_max)
        flagged = is_unplaced and is_nonplant and is_abnormalgc

        placed_str = 'unplaced' if is_unplaced else 'placed'
        audit.write(f"{name}\\t{info['length']}\\t{info['gc']}\\t{info['phylum']}\\t{placed_str}\\t{'yes' if flagged else 'no'}\\n")
        if flagged:
            ids.write(name + "\\n")
            flagged_ids.append(name)
            phyla_seen[info['phylum']] += 1
            bases_removed += info['length']

with open(summary_path, 'w') as s:
    s.write(f"contaminants_flagged={len(flagged_ids)}\\n")
    s.write(f"bases_removed={bases_removed}\\n")
    s.write("phyla_detected=" + ",".join(f"{p}:{c}" for p, c in phyla_seen.items()) + "\\n")
PYEOF
    """
}

// Mechanical FASTA subsetting — deliberately NOT errorStrategy 'ignore': the decontam FASTA is
// a hard dependency for the QUAST_NUCLEAR/BUSCO_NUCLEAR 'decontam' columns downstream, so a
// failure here should surface clearly rather than silently produce a misleading empty or
// zero-genome-fraction column.
// Second-opinion corroboration for the (typically very few) contigs CLASSIFY_CONTAMINANTS
// flagged: an independent containment-ANI method (Sylph) against real GTDB reference genomes,
// rather than Kraken2's k-mer LCA approach. Only queries the flagged candidates — cheap even
// though the GTDB database itself is large — so this is advisory extra confidence, not a gate
// on removal (REMOVE_CONTAMINANTS' behavior is unchanged either way).
process EXTRACT_CANDIDATE_CONTIGS {
    tag           { sample_id }
    label         'qc'
    errorStrategy 'ignore'
    container     'quay.io/biocontainers/seqkit:2.13.0--he881be0_0'

    input:
    tuple val(sample_id), path(assembly), path(contaminant_ids)

    output:
    tuple val(sample_id), path("${sample_id}_candidates.fasta"), emit: fasta

    script:
    """
    if [ -s ${contaminant_ids} ]; then
        seqkit grep -f ${contaminant_ids} ${assembly} -o ${sample_id}_candidates.fasta
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
    errorStrategy 'terminate'
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
    set -uo pipefail

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
    awk '/^>/{f="per_contig/" substr(\$1,2) ".fasta"} {print > f}' ${candidates}
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

// Advisory nucleotide-level corroboration for the same flagged candidates checked by Sylph.
// The database is deliberately supplied by the user (typically a locally formatted NCBI nt
// database) so results are reproducible and the five-genome run does not depend on live NCBI
// rate limits or an undocumented remote database version. BLAST evidence never changes removal.
process BLAST_VERIFY_CONTAMINANTS {
    tag           { sample_id }
    label         'qc'
    errorStrategy 'terminate'
    publishDir    { "${params.outdir}/qc/contamination/${sample_id}/blast" }, mode: 'copy'
    container     'quay.io/biocontainers/blast:2.15.0--pl5321h6f7f691_1'

    input:
    tuple val(sample_id), path(candidates)
    val blast_db_dir

    output:
    tuple val(sample_id), path("${sample_id}_blast_corroboration.tsv"), emit: corroboration
    tuple val(sample_id), path("${sample_id}_blast_summary.txt"),       emit: summary

    script:
    def sample = sample_id
    def min_id = params.blast_min_identity
    def min_qc = params.blast_min_qcov
    def max_t  = params.blast_max_targets
    def blast_target = blast_db_dir ? "-db ${blast_db_dir}/nt" : '-db nt -remote'
    def thread_args = blast_db_dir ? "-num_threads ${task.cpus}" : ''
    """
    #!/usr/bin/env bash
    set -euo pipefail

    if [ ! -s ${candidates} ]; then
        printf 'contig\\tbest_identity\\tbest_qcov\\tbest_subject\\tbest_title\\tevalue\\tbitscore\\tcorroborated\\n' > ${sample}_blast_corroboration.tsv
        echo 'blast_corroborated=0' > ${sample}_blast_summary.txt
        echo 'blast_total_candidates=0' >> ${sample}_blast_summary.txt
        exit 0
    fi

    blastn -query ${candidates} ${blast_target} \\
        ${thread_args} -max_target_seqs ${max_t} \\
        -evalue 1e-5 -outfmt '6 qseqid sacc stitle pident length qlen evalue bitscore' \\
        -out ${sample}_blast_query.tsv

    # Select the highest-bitscore hit per query. A hit corroborates the candidate only
    # when both identity and query coverage meet explicit thresholds.
    awk -F '\\t' -v OFS='\\t' -v minid=${min_id} -v minqc=${min_qc} \\
        'BEGIN { print "contig", "best_identity", "best_qcov", "best_subject", "best_title", "evalue", "bitscore", "corroborated" }
         FILENAME == ARGV[1] { q=\$1; qcov=100*\$5/\$6; if (!(q in best) || \$8 > best[q]) {
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

process REMOVE_CONTAMINANTS {
    tag        { sample_id }
    label      'qc'
    publishDir { "${params.outdir}/assembly/decontam/${sample_id}" }, mode: 'copy'
    container  'quay.io/biocontainers/seqkit:2.13.0--he881be0_0'

    input:
    tuple val(sample_id), path(assembly), path(contaminant_ids)

    output:
    tuple val(sample_id), path("${sample_id}_decontam.fasta"), emit: assembly

    script:
    """
    if [ -s ${contaminant_ids} ]; then
        seqkit grep -v -f ${contaminant_ids} ${assembly} -o ${sample_id}_decontam.fasta
    else
        cp ${assembly} ${sample_id}_decontam.fasta
    fi
    """
}
