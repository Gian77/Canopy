// ============================================================
// modules/qc.nf
// QC, read filtering, BUSCO, MultiQC
// ============================================================

process NANOPLOT {
    tag        { sample_id }
    label      'qc_light'
    publishDir { "${params.outdir}/qc/nanoplot/${sample_id}" }, mode: 'copy'
    container  'quay.io/biocontainers/nanoplot:1.43.0--pyhdfd78af_1'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("nanoplot_${sample_id}/*"), emit: report

    script:
    """
    NanoPlot \\
        --threads ${task.cpus} \\
        --fastq ${reads} \\
        --outdir nanoplot_${sample_id} \\
        --tsv_stats
    """
}

process FILTER_READS {
    tag        { sample_id }
    label      'qc_light'
    publishDir "${params.outdir}/reads_filtered", mode: 'copy'
    container  'quay.io/biocontainers/filtlong:0.2.1--hdcf5f25_3'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.filtered.fastq.gz"), emit: reads

    script:
    """
    cat ${reads} > combined.fastq.gz
    filtlong \\
        --min_length 1000 \\
        --min_mean_q 80 \\
        combined.fastq.gz \\
        | gzip > ${sample_id}.filtered.fastq.gz
    rm combined.fastq.gz
    """
}

process BUSCO_NUCLEAR {
    tag        { "${sample_id}_${stage}" }
    label      'qc_heavy'
    publishDir { "${params.outdir}/qc/busco/${sample_id}_${stage}" }, mode: 'copy'
    container  'quay.io/biocontainers/busco:5.8.0--pyhdfd78af_0'

    input:
    // stage ('medaka'|'purge') disambiguates the two candidate genomes BUSCO runs on so
    // their output dirs do not collide under qc/busco/.
    tuple val(sample_id), val(stage), path(assembly)

    output:
    tuple val(sample_id), val(stage), path("busco_${sample_id}_${stage}/*"),                  emit: full
    tuple val(sample_id), val(stage), path("busco_${sample_id}_${stage}/short_summary*.txt"), emit: summary

    script:
    """
    busco \\
        --in ${assembly} \\
        --out busco_${sample_id}_${stage} \\
        --mode genome \\
        --lineage_dataset ${params.busco_lineage} \\
        --cpu ${task.cpus}
    """
}

process MULTIQC {
    label         'qc'
    errorStrategy 'ignore'   // staging timeout on CephFS should not fail the whole run
    publishDir "${params.outdir}/multiqc", mode: 'symlink'
    container 'quay.io/biocontainers/multiqc:1.25.1--pyhdfd78af_0'

    input:
    path '*'

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc .
    """
}

process QUAST_ORGANELLE {
    tag       { "${sample_id}_${compartment}" }
    label     'assemble_small'
    publishDir { "${params.outdir}/qc/quast/${compartment}/${sample_id}" }, mode: 'copy'
    container 'quay.io/biocontainers/quast:5.3.0--py39pl5321h746d604_1'

    input:
    tuple val(sample_id), val(compartment), path(assembly), path(reference)

    output:
    tuple val(sample_id), val(compartment), path("quast_${sample_id}_${compartment}"),               emit: report
    tuple val(sample_id), val(compartment), path("quast_${sample_id}_${compartment}/report.tsv"),    emit: mqc

    script:
    """
    quast.py \\
        --reference ${reference} \\
        --threads ${task.cpus} \\
        --output-dir quast_${sample_id}_${compartment} \\
        --labels "${sample_id}_${compartment},reference" \\
        --min-contig 100 \\
        ${assembly} ${reference}
    """
}

process QUAST_NUCLEAR {
    tag       { sample_id }
    label     'qc'
    publishDir { "${params.outdir}/qc/quast/nuclear/${sample_id}" }, mode: 'copy'
    container 'quay.io/biocontainers/quast:5.3.0--py39pl5321h746d604_1'

    input:
    // Named assembly paths avoid the list-of-paths staging ambiguity. Unused slots
    // (e.g. asm_medaka_scaffold / asm_scaffold / asm_decontam when absent) receive a
    // passthrough and are excluded from the QUAST command via the has_* flags.
    tuple val(sample_id),
          path(asm_flye,            stageAs: 'asm_flye.fasta'),
          path(asm_medaka,          stageAs: 'asm_medaka.fasta'),
          path(asm_medaka_scaffold, stageAs: 'asm_medaka_scaffold.fasta'),
          path(asm_purge,           stageAs: 'asm_purge.fasta'),
          path(asm_scaffold,        stageAs: 'asm_scaffold.fasta'),
          path(asm_decontam,        stageAs: 'asm_decontam.fasta')
    val  has_scaffold          // true → include scaffold + reference
    val  has_medaka_scaffold   // true → also include the pre-purge (Medaka) RagTag scaffold
    val  has_decontam          // true → also include the contaminant-flagged final genome
    path reference             // pass [] when no nuclear_ref

    output:
    tuple val(sample_id), path("quast_${sample_id}_nuclear"),              emit: report
    tuple val(sample_id), path("quast_${sample_id}_nuclear/report.tsv"),   emit: mqc

    script:
    def ref_arg     = reference  ? "--reference ${reference}" : ""
    def ref_asm     = reference  ? " ${reference}" : ""
    def purge_label = "purge_dups"   // purge_dups always runs now
    // Build the column list in order: flye, medaka, [medaka_scaf], purge_dups, [purge_dups_scaf],
    // [decontam], [reference]. Both scaffolded columns carry their source stage in the name
    // ("medaka_scaf" / "purge_dups_scaf") so the report's qval lookups can target either one
    // unambiguously. decontam (contaminant-flagged final genome) lands last, right before
    // reference, since it's derived from whichever final scaffold was selected.
    def labels = ['flye', 'medaka']
    def asms   = ['asm_flye.fasta', 'asm_medaka.fasta']
    if (has_medaka_scaffold) { labels << 'medaka_scaf'; asms << 'asm_medaka_scaffold.fasta' }
    labels << purge_label;     asms << 'asm_purge.fasta'
    if (has_scaffold)        { labels << 'purge_dups_scaf'; asms << 'asm_scaffold.fasta' }
    if (has_decontam)        { labels << 'decontam';        asms << 'asm_decontam.fasta' }
    def label_str = labels.join(',') + (has_scaffold ? ',reference' : '')
    def asm_str   = asms.join(' ')
    """
    quast.py \\
        ${ref_arg} \\
        --threads ${task.cpus} \\
        --output-dir quast_${sample_id}_nuclear \\
        --labels "${label_str}" \\
        --min-contig 500 \\
        --large \\
        ${asm_str}${ref_asm}
    """
}

process BANDAGE_IMAGE {
    tag           { "${sample_id}_${compartment}" }
    cpus          1
    memory        '4 GB'
    time          '1h'
    errorStrategy 'ignore'   // visualization only — never fail the pipeline
    publishDir    { "${params.outdir}/assembly/${compartment}/${sample_id}" }, mode: 'copy'
    container     'quay.io/biocontainers/bandage:0.9.0--h9948957_0'

    input:
    tuple val(sample_id), val(compartment), path(gfa)

    output:
    tuple val(sample_id), val(compartment), path("${sample_id}_${compartment}_graph.png"),      optional: true, emit: image
    tuple val(sample_id), val(compartment), path("${sample_id}_${compartment}_graph_info.txt"), optional: true, emit: info

    script:
    """
    # Qt offscreen rendering — no X display needed on headless cluster nodes
    export QT_QPA_PLATFORM=offscreen
    Bandage image ${gfa} ${sample_id}_${compartment}_graph.png --scope entire 2>/dev/null || true
    Bandage info  ${gfa} > ${sample_id}_${compartment}_graph_info.txt 2>/dev/null || true
    """
}

process ALIGN_FOR_QC {
    tag       { "${sample_id}_${stage}" }
    label     'map'
    container 'quay.io/biocontainers/minimap2:2.31--h118bc1c_0'

    input:
    tuple val(sample_id), val(stage), path(assembly), path(reads)

    output:
    tuple val(sample_id), val(stage), path(assembly), path("${sample_id}_${stage}.sam"), emit: sam

    script:
    """
    minimap2 -t ${task.cpus} -ax map-ont ${assembly} ${reads} > ${sample_id}_${stage}.sam
    """
}

process SORT_FOR_QC {
    tag       { "${sample_id}_${stage}" }
    label     'map'
    publishDir { "${params.outdir}/qc/alignments/${sample_id}" }, mode: 'symlink'
    container 'quay.io/biocontainers/samtools:1.6--h5fe306e_13'

    input:
    tuple val(sample_id), val(stage), path(assembly), path(sam)

    output:
    tuple val(sample_id), val(stage), path(assembly), path("${sample_id}_${stage}.bam"), path("${sample_id}_${stage}.bam.bai"), emit: bam

    script:
    """
    samtools sort -@ ${task.cpus} -o ${sample_id}_${stage}.bam ${sam}
    samtools index ${sample_id}_${stage}.bam
    """
}

// BlobTools 1.x's own "% reads mapped" stat is internally inconsistent on BAMs with many
// secondary/supplementary alignments (routine for ONT + minimap2 on a repetitive genome):
// its numerator (BtIO.py check_mapped_read) correctly excludes secondary/supplementary per
// contig, but its denominator (BtIO.py checkBam: aln.mapped + aln.unmapped) comes straight
// from the BAM index and counts every alignment RECORD, secondary/supplementary included.
// With ~5x more secondary+supplementary records than primary alignments, that mismatch alone
// can make BlobTools report ~15% "mapped" when the true read-level rate is ~100%. Filtering
// to primary alignments only before BlobTools sees the BAM makes both sides of its fraction
// agree — it does not change per-contig coverage/GC/taxonomy (those already excluded
// secondary/supplementary), only the misleading summary percentage.
process FILTER_PRIMARY_BAM {
    tag       { "${sample_id}_${stage}" }
    label     'map'
    container 'quay.io/biocontainers/samtools:1.6--h5fe306e_13'

    input:
    tuple val(sample_id), val(stage), path(assembly), path(bam), path(bai)

    output:
    tuple val(sample_id), val(stage), path(assembly),
          path("${sample_id}_${stage}.primary.bam"),
          path("${sample_id}_${stage}.primary.bam.bai"), emit: bam

    script:
    """
    samtools view -@ ${task.cpus} -F 0x900 -b ${bam} > ${sample_id}_${stage}.primary.bam
    samtools index ${sample_id}_${stage}.primary.bam
    """
}

process QUALIMAP_BAMQC {
    tag       { "${sample_id}_${stage}" }
    label     'qc'
    publishDir { "${params.outdir}/qc/qualimap/${sample_id}" }, mode: 'copy'
    container 'quay.io/biocontainers/qualimap:2.3--hdfd78af_0'

    input:
    tuple val(sample_id), val(stage), path(assembly), path(bam), path(bai)

    output:
    tuple val(sample_id), val(stage), path("qualimap_${sample_id}_${stage}"), emit: report

    script:
    """
    qualimap bamqc \\
        -bam ${bam} \\
        -outdir qualimap_${sample_id}_${stage} \\
        -nt ${task.cpus} \\
        --java-mem-size=32G
    """
}

process BLOBTOOLS_COVERAGE {
    tag           { "${sample_id}_${stage}" }
    label         'qc'
    errorStrategy 'ignore'
    publishDir    { "${params.outdir}/qc/blobtools/${sample_id}" }, mode: 'copy'
    container     'quay.io/biocontainers/blobtools:1.1.1--py_1'

    input:
    tuple val(sample_id), val(stage), path(assembly), path(bam), path(bai)

    output:
    path "blobdb_${sample_id}_${stage}*", emit: results

    script:
    def pfx = "blobdb_${sample_id}_${stage}"
    """
    # blobtools v1 requires a nodesDB even in coverage-only mode (no hits file).
    # Build a minimal 2-node DB so the tool initialises; all contigs will appear as "no-hit".
    printf "# nodes_count = 2\\n1\\tno rank\\troot\\t1\\n0\\tno rank\\tno-hit\\t1\\n" > nodesDB.txt

    blobtools create -i ${assembly} -b ${bam} -o ${pfx} --db nodesDB.txt
    blobtools plot   -i ${pfx}.blobDB.json --out ${pfx}
    blobtools view   -i ${pfx}.blobDB.json -o ${pfx}
    """
}

process FILTER_ORGANELLE_CONTIGS {
    tag       { "${sample_id}_${compartment}" }
    label     'map'
    publishDir { "${params.outdir}/assembly/${compartment}/${sample_id}" }, mode: 'copy'
    container 'quay.io/biocontainers/minimap2:2.30--h577a1d6_0'

    input:
    tuple val(sample_id), val(compartment), path(raw_assembly), path(reference)

    output:
    tuple val(sample_id), val(compartment), path("${sample_id}_${compartment}_filtered.fasta"), emit: filtered
    tuple val(sample_id), val(compartment), path("${sample_id}_${compartment}_alignment.paf"),  emit: paf
    tuple val(sample_id), val(compartment), path("${sample_id}_${compartment}_kept_ids.txt"),    emit: ids

    script:
    """
    # Align raw assembly contigs to organelle reference.
    # PAF columns: 1=qname 2=qlen 3=qstart 4=qend 5=strand 6=tname 7=tlen 8=tstart 9=tend 10=matches 11=alnlen
    if [ -s ${raw_assembly} ]; then
        minimap2 -cx asm10 -t ${task.cpus} ${reference} ${raw_assembly} \\
            > ${sample_id}_${compartment}_alignment.paf
    else
        : > ${sample_id}_${compartment}_alignment.paf
    fi

    # Keep contigs where:
    #   - matched bases / contig length > 0.5 (covers majority of contig)
    #   - matched bases / aln length > 0.7 (high identity)
    # Summed across multiple alignments per contig (a contig may align in pieces).
    awk 'BEGIN{OFS="\\t"}
         {
             cov[\$1] += \$10
             aln[\$1] += \$11
             qlen[\$1] = \$2
         }
         END {
             for (c in cov) {
                 frac_q   = cov[c] / qlen[c]
                 identity = cov[c] / aln[c]
                 if (frac_q > 0.5 && identity > 0.7) print c
             }
         }' ${sample_id}_${compartment}_alignment.paf | sort -u \\
         > ${sample_id}_${compartment}_kept_ids.txt

    # Extract kept contigs (using awk; no seqkit dep needed)
    if [ -s ${sample_id}_${compartment}_kept_ids.txt ]; then
        awk 'BEGIN { while ((getline line < "'${sample_id}_${compartment}_kept_ids.txt'") > 0) keep[line]=1 }
             /^>/ { name = substr(\$1, 2); flag = (name in keep); if (flag) print; next }
             { if (flag) print }' \\
             ${raw_assembly} > ${sample_id}_${compartment}_filtered.fasta
    else
        echo "No contigs passed filter for ${compartment}" >&2
        : > ${sample_id}_${compartment}_filtered.fasta
    fi
    """
}
