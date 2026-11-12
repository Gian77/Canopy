// ============================================================
// modules/assembly.nf
// Organelles: Oatk, hifiasm and Flye (one is kept per compartment, see SELECT_ORGANELLE in qc.nf).
// Nuclear: hifiasm (default) or Flye.
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
    tuple val(sample_id), path("flye_cp_${sample_id}/assembly_info.txt"), optional: true, emit: info

    script:
    """
    # Flye finding nothing to assemble is a biological result, so emit an empty assembly. Any
    # other failure (e.g. SIGILL on a bad node) must fail the task so it can be retried.
    flye_status=0
    flye \\
        --nano-hq ${reads} \\
        --genome-size 150k \\
        --threads ${task.cpus} \\
        --out-dir flye_cp_${sample_id} \\
        --min-overlap 5000 \\
        --asm-coverage 50 || flye_status=\$?

    if [ -f flye_cp_${sample_id}/assembly.fasta ]; then
        cp flye_cp_${sample_id}/assembly.fasta ${sample_id}_cp_raw.fasta
    elif grep -qE 'No disjointigs|No contigs' flye_cp_${sample_id}/flye.log; then
        echo "Flye assembled nothing on chloroplast — emitting empty assembly" >&2
        touch ${sample_id}_cp_raw.fasta
    else
        echo "Flye failed on chloroplast (exit \$flye_status)" >&2
        exit \$(( flye_status == 0 ? 1 : flye_status ))
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
    tuple val(sample_id), path("flye_mt_${sample_id}/assembly_info.txt"), optional: true, emit: info

    script:
    """
    # See ASSEMBLE_CP_FLYE: only "nothing to assemble" produces an empty assembly.
    flye_status=0
    flye \\
        --nano-hq ${reads} \\
        --genome-size 500k \\
        --threads ${task.cpus} \\
        --out-dir flye_mt_${sample_id} \\
        --min-overlap 5000 || flye_status=\$?

    if [ -f flye_mt_${sample_id}/assembly.fasta ]; then
        cp flye_mt_${sample_id}/assembly.fasta ${sample_id}_mt_raw.fasta
    elif grep -qE 'No disjointigs|No contigs' flye_mt_${sample_id}/flye.log; then
        echo "Flye assembled nothing on mitochondria — emitting empty assembly" >&2
        touch ${sample_id}_mt_raw.fasta
    else
        echo "Flye failed on mitochondria (exit \$flye_status)" >&2
        exit \$(( flye_status == 0 ? 1 : flye_status ))
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

    # -k 501: ONT-appropriate overlap size; 1001 (HiFi default) fragments ONT graphs.
    # A biological sample may contain no recoverable organelle component.
    # Preserve stderr so that biological no-result can be distinguished from a
    # genuine Oatk/runtime failure.
    set +e
    oatk \\
        -k 501 \\
        -c 30 \\
        -t ${task.cpus} \\
        -m ../${mito_db} \\
        -p ../${pltd_db} \\
        -o ${sample_id} \\
        ../${reads} 2> >(tee oatk.stderr >&2)
    oatk_status=\$?
    set -e

    if [ "\${oatk_status}" -ne 0 ]; then
        if grep -q 'no organelle component found' oatk.stderr; then
            echo "OATK found no organelle component; emitting empty assemblies" >&2
        else
            exit "\${oatk_status}"
        fi
    fi

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

// Organelle reads are thousands of x deep and hifiasm is built for nuclear depth, so cap them at about 200x.
process DOWNSAMPLE_ORGANELLE_READS {
    tag       { "${sample_id}_${compartment}" }
    label     'qc'
    container 'quay.io/biocontainers/seqkit:2.13.0--he881be0_0'

    input:
    tuple val(sample_id), val(compartment), path(reads), val(target_bases)

    output:
    tuple val(sample_id), val(compartment), path("${sample_id}_${compartment}.sub.fastq.gz")

    script:
    """
    set -euo pipefail
    total=\$(seqkit stats -T ${reads} | awk 'NR == 2 { print \$5 }')
    p=\$(awk -v t=${target_bases} -v n="\$total" 'BEGIN { x = t / n; if (x > 1) x = 1; printf "%.6f", x }')
    echo "total bases \$total  sampling proportion \$p"
    seqkit sample -s 11 -p "\$p" -o ${sample_id}_${compartment}.sub.fastq.gz ${reads}
    """
}

// -f0: no bloom filter, as hifiasm's documentation advises for small genomes. -l0: no purging.
process ASSEMBLE_ORGANELLE_HIFIASM {
    tag        { "${sample_id}_${compartment}" }
    label      'assemble_small'
    publishDir { "${params.outdir}/assembly/${compartment}/${sample_id}/candidates/raw" }, mode: 'copy', pattern: '*.{fasta,log}'
    container  'quay.io/biocontainers/hifiasm:0.25.0--h5ca1c30_0'

    input:
    tuple val(sample_id), val(compartment), path(reads)

    output:
    tuple val(sample_id), val(compartment), path("${sample_id}_${compartment}.hifiasm.fasta"), emit: assembly
    path "${sample_id}_${compartment}.hifiasm.log", emit: log

    script:
    """
    set -euo pipefail
    hifiasm -o asm -t ${task.cpus} --ont -l0 -f0 ${reads} 2> ${sample_id}_${compartment}.hifiasm.log
    awk '/^S/ { print ">" \$2; print \$3 }' asm.bp.p_ctg.gfa > ${sample_id}_${compartment}.hifiasm.fasta
    rm -f asm.*.bin
    """
}

// hifiasm ONT mode (R10 simplex reads). Primary contigs only. --hifiasm_purge_level is 0 by default
// because purge_dups runs later and the sorghum lines are inbred.
process ASSEMBLE_NUCLEAR_HIFIASM {
    tag       { sample_id }
    label     'assemble_heavy'
    publishDir { "${params.outdir}/assembly/nuclear/${sample_id}" }, mode: 'copy'
    container 'quay.io/biocontainers/hifiasm:0.25.0--h5ca1c30_0'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_nuclear.fasta"), emit: assembly
    path  "${sample_id}.hifiasm.log",                          emit: log

    script:
    """
    set -euo pipefail
    hifiasm -o asm -t ${task.cpus} --ont -l${params.hifiasm_purge_level} ${reads} 2> ${sample_id}.hifiasm.log
    awk '/^S/ { print ">" \$2; print \$3 }' asm.bp.p_ctg.gfa > ${sample_id}_nuclear.fasta
    [ -s ${sample_id}_nuclear.fasta ]
    rm -f asm.*.bin
    """
}

process ASSEMBLE_NUCLEAR_FLYE {
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
