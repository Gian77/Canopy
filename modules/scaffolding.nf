// Chromosome scaffolding with RagTag.
// Runs after purge_dups: corrects chimeric contigs against the nuclear reference,
// then orders/orients them into reference-based chromosome groups, bridging
// neighbouring contigs with a 100 bp N-gap (RagTag's default scaffold gap).

process RAGTAG_SCAFFOLD {
    tag       { "${sample_id}_${stage}" }
    label     'assemble_small'
    publishDir { "${params.outdir}/assembly/scaffolds/${sample_id}/${stage}" }, mode: 'copy'
    container 'quay.io/biocontainers/ragtag:2.1.0--pyhdfd78af_2'

    input:
    tuple val(sample_id), val(stage), path(assembly)
    path reference

    output:
    tuple val(sample_id), path("${sample_id}_${stage}_scaffold.fasta"),                        emit: scaffold
    tuple val(sample_id), path("${sample_id}_${stage}_scaffold.agp"),                          emit: agp
    tuple val(sample_id), path("ragtag_${sample_id}_${stage}/scaffold/ragtag.scaffold.stats"), emit: stats
    path "ragtag_${sample_id}_${stage}",                                                        emit: log

    script:
    def rt = "ragtag_${sample_id}_${stage}"
    def do_correct   = params.ragtag_correct
    def scaffold_input = do_correct ? "${rt}/correct/ragtag.correct.fasta"
                                    : "${assembly}"
    """
    mkdir -p ${rt}/correct

    ${do_correct ? """
    # Correction: split chimeric contigs against the reference
    ragtag.py correct ${reference} ${assembly} \\
        -t ${task.cpus} \\
        -o ${rt}/correct
    """ : "# correction skipped (params.ragtag_correct = false)"}

    # Scaffold: order & orient contigs into chromosome-scale groups
    ragtag.py scaffold ${reference} ${scaffold_input} \\
        -t ${task.cpus} \\
        -o ${rt}/scaffold

    cp ${rt}/scaffold/ragtag.scaffold.fasta ${sample_id}_${stage}_scaffold.fasta
    cp ${rt}/scaffold/ragtag.scaffold.agp   ${sample_id}_${stage}_scaffold.agp
    """
}
