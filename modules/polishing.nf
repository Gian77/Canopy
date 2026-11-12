// ============================================================
// modules/polishing.nf
// Medaka polishing, purge_dups, HapDup phasing
// ============================================================

process POLISH_MEDAKA {
    tag       { sample_id }
    label     'polish'
    publishDir { "${params.outdir}/polishing/medaka/${sample_id}" }, mode: 'copy'
    container 'quay.io/biocontainers/medaka:2.2.1--py312hc7af5e1_0'

    input:
    tuple val(sample_id), path(assembly), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_polished.fasta"), emit: assembly
    path "medaka_${sample_id}/*",                               emit: log

    script:
    """
    medaka_consensus \\
        -i ${reads} \\
        -d ${assembly} \\
        -o medaka_${sample_id} \\
        -t ${task.cpus} \\
        -m ${params.medaka_model}

    cp medaka_${sample_id}/consensus.fasta ${sample_id}_polished.fasta
    """
}

process PURGE_DUPS {
    tag       { sample_id }
    label     'polish'
    publishDir { "${params.outdir}/polishing/purge_dups/${sample_id}" }, mode: 'copy'
    container 'quay.io/biocontainers/purge_dups:1.2.6--h577a1d6_3'

    input:
    tuple val(sample_id), path(assembly), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_purged.fasta"),                       emit: assembly
    tuple val(sample_id), path("purge_dups_${sample_id}/cutoffs"),                 emit: cutoffs
    tuple val(sample_id), path("purge_dups_${sample_id}/calcuts.log"),             emit: calcuts_log
    tuple val(sample_id), path("purge_dups_${sample_id}/${sample_id}.PB.stat"),    emit: pbstat
    path  "purge_dups_${sample_id}/*",                                              emit: log

    script:
    def calcuts_extra = params.calcuts_args ?: ""
    """
    mkdir -p purge_dups_${sample_id}
    cd purge_dups_${sample_id}

    # Step 1: map reads to assembly to compute depth
    minimap2 -t ${task.cpus} -xmap-ont ../${assembly} ../${reads} \\
        | gzip -c > reads.paf.gz
    pbcstat reads.paf.gz
    cp PB.stat ${sample_id}.PB.stat

    # Set calcuts thresholds: manual args if provided, else autotune from PB.stat mode
    if [ -n "${calcuts_extra}" ]; then
        calcuts ${calcuts_extra} PB.stat > cutoffs 2> calcuts.log
    else
        PEAK=\$(awk '\$1>=10 {if(\$2>max){max=\$2; peak=\$1}} END{print peak}' PB.stat)
        JUNK=\$(( PEAK / 4 ))
        HAP_LOW=\$(( PEAK / 2 ))
        HAP_HIGH=\$(( PEAK + PEAK / 2 ))
        DIP_LOW=\$(( HAP_HIGH + 1 ))
        DIP_HIGH=\$(( PEAK * 3 ))
        REPEAT=\$(( PEAK * 4 ))
        printf "%d\t%d\t%d\t%d\t%d\t%d\n" \${JUNK} \${HAP_LOW} \${HAP_HIGH} \${DIP_LOW} \${DIP_HIGH} \${REPEAT} > cutoffs
        echo "[autotune] Peak: \${PEAK}x  cutoffs: \${JUNK} \${HAP_LOW} \${HAP_HIGH} \${DIP_LOW} \${DIP_HIGH} \${REPEAT}" > calcuts.log
    fi

    # Step 2: self-align contigs to find duplicates
    split_fa ../${assembly} > assembly.split
    minimap2 -t ${task.cpus} -xasm5 -DP assembly.split assembly.split \\
        | gzip -c > assembly.self.paf.gz

    # Step 3: purge
    purge_dups -2 -T cutoffs -c PB.base.cov assembly.self.paf.gz > dups.bed 2> purge_dups.log
    get_seqs -e dups.bed ../${assembly}

    cd ..
    cp purge_dups_${sample_id}/purged.fa ${sample_id}_purged.fasta
    """
}

process ALIGN_FOR_HAPDUP {
    tag       { sample_id }
    label     'map'
    container 'quay.io/biocontainers/minimap2:2.31--h118bc1c_0'

    input:
    tuple val(sample_id), path(assembly), path(reads)

    output:
    tuple val(sample_id), path(assembly), path("aligned.sam"), emit: sam

    script:
    """
    minimap2 -t ${task.cpus} -ax map-ont ${assembly} ${reads} > aligned.sam
    """
}

process SORT_FOR_HAPDUP {
    tag       { sample_id }
    label     'map'
    container 'quay.io/biocontainers/samtools:1.6--h5fe306e_13'

    input:
    tuple val(sample_id), path(assembly), path(sam)

    output:
    tuple val(sample_id), path(assembly), path("aligned.bam"), path("aligned.bam.bai"), emit: bam

    script:
    """
    samtools sort -@ ${task.cpus} -o aligned.bam ${sam}
    samtools index aligned.bam
    """
}

process HAPDUP {
    tag       { sample_id }
    label     'assemble_heavy'
    publishDir { "${params.outdir}/polishing/hapdup/${sample_id}" }, mode: 'copy'
    container 'docker://mkolmogo/hapdup:0.12'

    input:
    tuple val(sample_id), path(assembly), path(bam), path(bai)

    output:
    tuple val(sample_id), path("hapdup_${sample_id}/hapdup_dual_*.fasta"), emit: haplotypes
    path "hapdup_${sample_id}/*",                                           emit: log

    script:
    """
    hapdup \\
        --assembly ${assembly} \\
        --bam ${bam} \\
        --out-dir hapdup_${sample_id} \\
        --threads ${task.cpus} \\
        --rtype ont
    """
}
