// ============================================================
// modules/assembly.nf
// Flye assemblies (cp/mt/nuclear) and OATK alternative for organelles
// ============================================================

process ASSEMBLE_CP_FLYE {
    tag       { sample_id }
    label     'assemble_small'
    publishDir { "${params.outdir}/assembly/chloroplast/${sample_id}/raw" }, mode: 'copy'
    container 'quay.io/biocontainers/flye:2.9.4--py310h2b6aa90_0'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_cp_raw.fasta"), emit: assembly
    path  "flye_cp_${sample_id}/*",                          emit: log

    script:
    """
    flye \\
        --nano-hq ${reads} \\
        --genome-size 150k \\
        --threads ${task.cpus} \\
        --out-dir flye_cp_${sample_id} \\
        --min-overlap 5000 \\
        --asm-coverage 50 || true

    if [ -f flye_cp_${sample_id}/assembly.fasta ]; then
        cp flye_cp_${sample_id}/assembly.fasta ${sample_id}_cp_raw.fasta
    else
        echo "Flye failed on chloroplast — emitting empty assembly" >&2
        touch ${sample_id}_cp_raw.fasta
    fi
    """
}

process ASSEMBLE_MT_FLYE {
    tag       { sample_id }
    label     'assemble_small'
    publishDir { "${params.outdir}/assembly/mitochondria/${sample_id}/raw" }, mode: 'copy'
    container 'quay.io/biocontainers/flye:2.9.4--py310h2b6aa90_0'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_mt_raw.fasta"), emit: assembly
    path  "flye_mt_${sample_id}/*",                          emit: log

    script:
    """
    flye \\
        --nano-hq ${reads} \\
        --genome-size 500k \\
        --threads ${task.cpus} \\
        --out-dir flye_mt_${sample_id} \\
        --min-overlap 5000 || true

    if [ -f flye_mt_${sample_id}/assembly.fasta ]; then
        cp flye_mt_${sample_id}/assembly.fasta ${sample_id}_mt_raw.fasta
    else
        echo "Flye failed on mitochondria — emitting empty assembly" >&2
        touch ${sample_id}_mt_raw.fasta
    fi
    """
}

process ASSEMBLE_ORGANELLES_OATK {
    tag       { sample_id }
    label     'assemble_small'
    publishDir { "${params.outdir}/assembly/oatk/${sample_id}" }, mode: 'copy'
    container 'docker.io/assteindorff/oatk:1.0'

    input:
    tuple val(sample_id), path(reads)
    tuple path(mito_db), path(mito_idx)
    tuple path(pltd_db), path(pltd_idx)

    output:
    tuple val(sample_id), path("${sample_id}_cp_raw.fasta"),                                     emit: cp_assembly
    tuple val(sample_id), path("${sample_id}_mt_raw.fasta"),                                     emit: mt_assembly
    tuple val(sample_id), path("oatk_${sample_id}/${sample_id}.pltd.gfa"), optional: true,      emit: cp_gfa
    tuple val(sample_id), path("oatk_${sample_id}/${sample_id}.mito.gfa"), optional: true,      emit: mt_gfa
    path "oatk_${sample_id}/*", emit: log, optional: true

    script:
    """
    mkdir -p oatk_${sample_id}
    cd oatk_${sample_id}

    # -k 501: ONT-appropriate overlap size; 1001 (HiFi default) fragments ONT graphs
    oatk \\
        -k 501 \\
        -c 30 \\
        -t ${task.cpus} \\
        -m ../${mito_db} \\
        -p ../${pltd_db} \\
        -o ${sample_id} \\
        ../${reads}

    cd ..

    # Copy outputs to canonical names (handle empty/missing gracefully)
    if [ -f oatk_${sample_id}/${sample_id}.pltd.ctg.fasta ]; then
        cp oatk_${sample_id}/${sample_id}.pltd.ctg.fasta ${sample_id}_cp_raw.fasta
    else
        echo "OATK produced no chloroplast assembly" >&2
        touch ${sample_id}_cp_raw.fasta
    fi

    if [ -f oatk_${sample_id}/${sample_id}.mito.ctg.fasta ]; then
        cp oatk_${sample_id}/${sample_id}.mito.ctg.fasta ${sample_id}_mt_raw.fasta
    else
        echo "OATK produced no mitochondrial assembly" >&2
        touch ${sample_id}_mt_raw.fasta
    fi
    """
}

process ASSEMBLE_NUCLEAR {
    tag       { sample_id }
    label     'assemble_heavy'
    publishDir { "${params.outdir}/assembly/nuclear/${sample_id}" }, mode: 'copy'
    container 'quay.io/biocontainers/flye:2.9.4--py310h2b6aa90_0'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_nuclear.fasta"), emit: assembly
    path  "flye_nuclear_${sample_id}/*",                       emit: log

    script:
    """
    flye \\
        --nano-hq ${reads} \\
        --genome-size ${params.genome_size} \\
        --threads ${task.cpus} \\
        --out-dir flye_nuclear_${sample_id}

    cp flye_nuclear_${sample_id}/assembly.fasta ${sample_id}_nuclear.fasta
    """
}

