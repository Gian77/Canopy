// ============================================================
// modules/trimming.nf
// Optional adaptor trimming of raw ONT reads before Filtlong (--read_trimming).
// ============================================================
//
// Off by default (--read_trimming none): the raw reads go straight to Filtlong, as before.
// Each tool removes adaptors (and barcodes) only; Filtlong still does all length and quality
// filtering. On BTx2932 (2026-11-06) all three removed about 98% of the adaptor that otherwise
// reaches contig ends, but none improved the assembly, and FCS-adaptor trimming at the end of the
// pipeline is still needed for NCBI. See README "Read trimming".
//
// --trim_adapters FASTA (porechop_abi, fastplong): records named <set>_start and <set>_end give
// the sequence expected at read starts and read ends. Use the full adaptor + barcode sequences of
// your kits; without them Porechop_ABI infers adaptors ab initio, which on Rapid-kit reads found
// only ~45 bp and left the adaptor on 86% of reads, and fastplong auto-detects one adaptor per file.

// Dorado demux with the kit of each run. --dorado_kit is one kit name, or a comma list of
// barcode=KIT entries plus one bare KIT for all other reads, for files that mix runs of
// different kits (e.g. "barcode18=SQK-NBD114-24,SQK-RBK114-24"; the barcode is matched in the
// read header's barcode= tag). Unclassified reads are kept, so no reads are dropped.
process TRIM_READS_DORADO {
    tag        { sample_id }
    label      'qc_heavy'
    publishDir "${params.outdir}/reads_trimmed", mode: 'copy', pattern: '*.trim_log.txt'
    container  'docker.io/ontresearch/dorado:sha9809639e07a927bcc0f584dadd5e59674cf59f3f'   // Dorado 2.1.1

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.trimmed.fastq.gz"), emit: reads
    path "${sample_id}.trim_log.txt", emit: log

    script:
    """
    set -euo pipefail
    zcat -f ${reads} > all.fastq
    # Split the reads by kit, following --dorado_kit.
    awk -v spec='${params.dorado_kit}' '
        BEGIN { n = split(spec, e, ","); for (i = 1; i <= n; i++) {
                    if (e[i] ~ /=/) { split(e[i], kv, "="); kit[kv[1]] = kv[2] } else { dflt = e[i] } } }
        NR % 4 == 1 { k = dflt; if (match(\$0, /barcode=[^ \\t]+/)) { b = substr(\$0, RSTART + 8, RLENGTH - 8); if (b in kit) k = kit[b] }
                      if (k == "") { print "ERROR: no kit for read " \$1 > "/dev/stderr"; exit 1 }
                      out = "kit_" k ".fastq" }
        { print > out }' all.fastq
    rm all.fastq
    : > ${sample_id}.trimmed.fastq
    for f in kit_*.fastq; do
        k=\${f#kit_}; k=\${k%.fastq}
        echo "dorado demux --kit-name \$k: \$(awk 'END {print NR / 4}' \$f) reads" >> ${sample_id}.trim_log.txt
        # Drop Dorado's progress counter ("Output records written", rewritten with carriage returns).
        dorado demux --kit-name \$k --emit-fastq -t ${task.cpus} -o demux_\$k \$f 2>&1 \\
            | tr '\\r' '\\n' | awk '!/Output records written/ && NF' >> ${sample_id}.trim_log.txt
        # Dorado writes SAM-style tags after a tab in the read name; keep the original header.
        find demux_\$k -name '*.fastq*' -print0 | xargs -0 zcat -f \\
            | awk 'NR % 4 == 1 { sub(/\\t.*/, "") } { print }' >> ${sample_id}.trimmed.fastq
        rm -rf demux_\$k \$f
    done
    n_in=\$(zcat -f ${reads} | awk 'END {print NR / 4}')
    n_out=\$(awk 'END {print NR / 4}' ${sample_id}.trimmed.fastq)
    echo "reads in: \$n_in  out: \$n_out" >> ${sample_id}.trim_log.txt
    [ "\$n_out" -eq "\$n_in" ] || { echo "ERROR: Dorado output has \$n_out of \$n_in reads" >&2; exit 1; }
    gzip ${sample_id}.trimmed.fastq
    """
}

// Porechop_ABI. With --trim_adapters it uses only those adaptors (-ddb -cap) in every chunk
// (--adapter_threshold 0: by default it picks adaptor sets from the first 10,000 reads, which skips
// a kit whose runs come later in the file). Without it, ab initio inference (-abi) plus Porechop's
// own database. Reads with a middle adaptor are split (Porechop's default).
// 0.5.1 builds die with SIGILL on scarcity-8 (no AVX2); 0.5.0post1 runs on every node.
process TRIM_READS_PORECHOP {
    tag        { sample_id }
    label      'qc_heavy'
    publishDir "${params.outdir}/reads_trimmed", mode: 'copy', pattern: '*.trim_log.txt'
    container  'quay.io/biocontainers/porechop_abi:0.5.0post1--py310h275bdba_0'

    input:
    tuple val(sample_id), path(reads)
    path adapters

    output:
    tuple val(sample_id), path("${sample_id}.trimmed.fastq.gz"), emit: reads
    path "${sample_id}.trim_log.txt", emit: log

    script:
    def custom = adapters.name != 'NO_ADAPTERS'
    """
    set -euo pipefail
    zcat -f ${reads} | gzip -1 > combined.fastq.gz
    ${custom ? """
    # FASTA <set>_start / <set>_end records to Porechop's format: name, start sequence, end sequence.
    awk '/^>/ { name = substr(\$1, 2); set = name; sub(/_(start|end)\$/, "", set)
                end = (name ~ /_end\$/); if (!(set in seen)) { seen[set] = 1; order[++n] = set }; next }
         { if (end) e[set] = e[set] \$0; else s[set] = s[set] \$0 }
         END { for (i = 1; i <= n; i++) { k = order[i]; print k; print s[k]; if (e[k] != "") print e[k]; print "" } }' \\
        ${adapters} > adapters_porechop.txt
    """ : ''}
    mkdir -p tmp
    porechop_abi ${custom ? '-ddb -cap adapters_porechop.txt --adapter_threshold 0' : '-abi'} \\
        -i combined.fastq.gz -o ${sample_id}.trimmed.fastq.gz -t ${task.cpus} -tmp tmp -v 1 \\
        > ${sample_id}.trim_log.txt 2>&1
    rm -rf tmp combined.fastq.gz
    """
}

// fastplong, adaptor trimming only: quality and length filtering are off (-Q -L), Filtlong does those.
process TRIM_READS_FASTPLONG {
    tag        { sample_id }
    label      'qc_heavy'
    publishDir "${params.outdir}/reads_trimmed", mode: 'copy', pattern: '*.{trim_log.txt,fastplong.json,fastplong.html}'
    container  'quay.io/biocontainers/fastplong:0.7.0--h43da1c4_0'

    input:
    tuple val(sample_id), path(reads)
    path adapters

    output:
    tuple val(sample_id), path("${sample_id}.trimmed.fastq.gz"), emit: reads
    path "${sample_id}.trim_log.txt", emit: log
    path "${sample_id}.fastplong.{json,html}"

    script:
    def custom = adapters.name != 'NO_ADAPTERS'
    """
    set -euo pipefail
    zcat -f ${reads} | gzip -1 > combined.fastq.gz
    fastplong -i combined.fastq.gz -o ${sample_id}.trimmed.fastq.gz ${custom ? "-a ${adapters}" : ''} \\
        -Q -L -w ${task.cpus} -z 4 -j ${sample_id}.fastplong.json -h ${sample_id}.fastplong.html \\
        > ${sample_id}.trim_log.txt 2>&1
    rm combined.fastq.gz
    """
}
