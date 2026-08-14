// ============================================================
// modules/mapping.nf
// Map ONT reads to combined cp+mt reference, extract per-compartment
// FASTQ, dedup, and report stats. Single-tool containers only.
// ============================================================
 
// ─────────────────────────────────────────────────────────────
// Step 1: Align reads to combined cp+mt reference with minimap2
// ─────────────────────────────────────────────────────────────
process ALIGN_TO_ORGANELLES {
    tag       { sample_id }
    label     'map'
    container 'quay.io/biocontainers/minimap2:2.30--h577a1d6_0'

    input:
    tuple val(sample_id), path(reads)
    path  cp_ref
    path  mt_ref

    output:
    tuple val(sample_id), path("${sample_id}.organelles.sam"),
                          path("organelle_ref.fasta"),       emit: sam

    script:
    """
    # Build combined organelle reference (cp + mt) with renamed headers
    awk 'BEGIN{n=0} /^>/{n++; print ">cp_"n; next} {print}' ${cp_ref} >  organelle_ref.fasta
    awk 'BEGIN{n=0} /^>/{n++; print ">mt_"n; next} {print}' ${mt_ref} >> organelle_ref.fasta

    minimap2 \\
        -t ${task.cpus} \\
        -ax map-ont \\
        --secondary=no \\
        organelle_ref.fasta ${reads} \\
        > ${sample_id}.organelles.sam
    """
}

// ─────────────────────────────────────────────────────────────
// Step 2: Sort, index, and flagstat the SAM → BAM
// ─────────────────────────────────────────────────────────────
process SORT_INDEX_BAM {
    tag       { sample_id }
    label     'map'
    publishDir { "${params.outdir}/mapping/${sample_id}" }, mode: 'copy', pattern: "*.flagstat"
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'

    input:
    tuple val(sample_id), path(sam), path(ref)

    output:
    tuple val(sample_id), path("${sample_id}.organelles.bam"),
                          path("${sample_id}.organelles.bam.bai"),
                          path("organelle_ref.fasta"),       emit: bam
    path  "${sample_id}.flagstat",                            emit: stats

    script:
    """
    # Sort, index, and flagstat. Reference file passed through for downstream use.
    cp ${ref} organelle_ref.fasta
    samtools sort -@ ${task.cpus} -o ${sample_id}.organelles.bam ${sam}
    samtools index -@ ${task.cpus} ${sample_id}.organelles.bam
    samtools flagstat ${sample_id}.organelles.bam > ${sample_id}.flagstat
    """
}


// ─────────────────────────────────────────────────────────────
// Step 3: extract per-compartment FASTQ (may contain duplicate
//          read IDs due to supplementary alignments).
// ─────────────────────────────────────────────────────────────
process EXTRACT_RAW_READSETS {
    tag       { sample_id }
    label     'map'
    publishDir { "${params.outdir}/reads_partitioned/${sample_id}" }, mode: 'copy', pattern: "*.nuclear.fastq.gz"
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'

    input:
    tuple val(sample_id), path(bam), path(bai), path(ref)

    output:
    tuple val(sample_id), path("${sample_id}.cp.raw.fastq.gz"),  emit: cp_raw
    tuple val(sample_id), path("${sample_id}.mt.raw.fastq.gz"),  emit: mt_raw
    tuple val(sample_id), path("${sample_id}.nuclear.fastq.gz"), emit: nuclear_reads

    script:
    """
    # Helper: extract primary alignments for contigs matching a prefix as FASTQ.
    # -F 0x904 = exclude unmapped (4) + secondary (256) + supplementary (2048).
    extract_raw() {
        local prefix=\$1
        local out=\$2
        local CONTIGS=\$(samtools view -H ${bam} | awk -v p="SN:\${prefix}" '/^@SQ/ && \$2 ~ p {sub(/SN:/,"",\$2); print \$2}')
        if [ -n "\$CONTIGS" ]; then
            samtools view -@ ${task.cpus} -b -F 0x904 ${bam} \$CONTIGS \\
              | samtools fastq -n - 2>/dev/null \\
              | gzip > \$out
        else
            echo "No \${prefix} contigs found in BAM header" >&2
            : | gzip > \$out
        fi
    }

    extract_raw cp_ ${sample_id}.cp.raw.fastq.gz
    extract_raw mt_ ${sample_id}.mt.raw.fastq.gz

    # Unmapped reads = nuclear candidates; no duplicates from supplementary alignments
    samtools fastq -f 4 ${bam} 2>/dev/null | gzip > ${sample_id}.nuclear.fastq.gz
    """
}

// ─────────────────────────────────────────────────────────────
// Step 4: dedup the cp and mt FASTQs by read name with seqkit.
// ─────────────────────────────────────────────────────────────
process DEDUP_ORGANELLE_READS {
    tag       { sample_id }
    label     'qc'
    publishDir { "${params.outdir}/reads_partitioned/${sample_id}" }, mode: 'copy'
    container 'quay.io/biocontainers/seqkit:2.13.0--he881be0_0'

    input:
    tuple val(sample_id), path(cp_raw), path(mt_raw)

    output:
    tuple val(sample_id), path("${sample_id}.cp.fastq.gz"), emit: cp_reads
    tuple val(sample_id), path("${sample_id}.mt.fastq.gz"), emit: mt_reads

    script:
    """
    # Dedup chloroplast reads by name (handles supplementary-alignment duplicates)
    if [ -s ${cp_raw} ]; then
        seqkit rmdup -n ${cp_raw} -o ${sample_id}.cp.fastq.gz
    else
        cp ${cp_raw} ${sample_id}.cp.fastq.gz
    fi

    # Dedup mitochondrial reads by name
    if [ -s ${mt_raw} ]; then
        seqkit rmdup -n ${mt_raw} -o ${sample_id}.mt.fastq.gz
    else
        cp ${mt_raw} ${sample_id}.mt.fastq.gz
    fi
    """
}

// ─────────────────────────────────────────────────────────────
// Step 5: per-compartment stats and estimated coverage.
// ─────────────────────────────────────────────────────────────
process READSET_STATS {
    tag       { sample_id }
    label     'qc'
    publishDir { "${params.outdir}/qc/readset_stats/${sample_id}" }, mode: 'copy'
    container 'quay.io/biocontainers/seqkit:2.13.0--he881be0_0'

    input:
    tuple val(sample_id),
          path(cp_reads),
          path(mt_reads),
          path(nuclear_reads),
          path(cp_ref),
          path(mt_ref)

    output:
    tuple val(sample_id), path("${sample_id}.readset_stats.tsv"),     emit: stats
    tuple val(sample_id), path("${sample_id}.readset_stats_mqc.tsv"), emit: mqc

    script:
    """
    cp_size=\$(seqkit stats -T ${cp_ref} | awk 'NR==2 {gsub(",","",\$5); print \$5}')
    mt_size=\$(seqkit stats -T ${mt_ref} | awk 'NR==2 {gsub(",","",\$5); print \$5}')

    # Convert nuclear genome-size param (e.g. "720m") into base pairs
    nuc_bp=\$(echo "${params.genome_size}" | awk '{
        n=\$1; u=tolower(substr(n, length(n))); v=substr(n, 1, length(n)-1)
        if (u=="k")      print v*1e3
        else if (u=="m") print v*1e6
        else if (u=="g") print v*1e9
        else             print n
    }')

    get_stat() {
        local f=\$1
        local field=\$2
        if [ -s "\$f" ]; then
            seqkit stats -T \$f | awk -v c=\$field 'NR==2 {gsub(",","",\$c); print \$c}'
        else
            echo 0
        fi
    }

    cp_reads_n=\$( get_stat ${cp_reads}      4)
    cp_bases=\$(   get_stat ${cp_reads}      5)
    mt_reads_n=\$( get_stat ${mt_reads}      4)
    mt_bases=\$(   get_stat ${mt_reads}      5)
    nuc_reads_n=\$(get_stat ${nuclear_reads} 4)
    nuc_bases=\$(  get_stat ${nuclear_reads} 5)

    cp_depth=\$( awk -v b=\$cp_bases  -v s=\$cp_size  'BEGIN { if (s>0) printf "%.0f", b/s; else print "NA" }')
    mt_depth=\$( awk -v b=\$mt_bases  -v s=\$mt_size  'BEGIN { if (s>0) printf "%.0f", b/s; else print "NA" }')
    nuc_depth=\$(awk -v b=\$nuc_bases -v s=\$nuc_bp   'BEGIN { if (s>0) printf "%.0f", b/s; else print "NA" }')

    {
        echo -e "compartment\\treads\\tbases\\tref_size\\tdepth_x"
        echo -e "chloroplast\\t\$cp_reads_n\\t\$cp_bases\\t\$cp_size\\t\$cp_depth"
        echo -e "mitochondria\\t\$mt_reads_n\\t\$mt_bases\\t\$mt_size\\t\$mt_depth"
        echo -e "nuclear\\t\$nuc_reads_n\\t\$nuc_bases\\t\$nuc_bp\\t\$nuc_depth"
    } > ${sample_id}.readset_stats.tsv

    {
        echo "# id: 'readset_stats_${sample_id}'"
        echo "# section_name: 'Readset partition (${sample_id})'"
        echo "# description: 'ONT reads partitioned by compartment with estimated coverage'"
        echo "# format: 'tsv'"
        echo "# plot_type: 'table'"
        cat ${sample_id}.readset_stats.tsv
    } > ${sample_id}.readset_stats_mqc.tsv
    """
}
