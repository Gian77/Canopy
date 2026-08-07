// modules/db.nf
// ============================================================
// Reference database fetchers.
// Processes here download external reference data needed by the pipeline.
// They use `storeDir` so each version is cached persistently per
// outdir — re-running the pipeline reuses cached databases instead of
// re-downloading. Versions are pinned via params for reproducibility.
// ============================================================


process DOWNLOAD_OATKDB {
    label    'fetch'
    storeDir "${params.outdir}/databases/oatkdb/${params.oatkdb_version}/raw"
    // no container — uses host curl

    output:
    path "embryophyta_mito.fam", emit: mito
    path "embryophyta_pltd.fam", emit: pltd

    script:
    def base = "https://raw.githubusercontent.com/c-zhou/OatkDB/${params.oatkdb_commit}/${params.oatkdb_version}"
    """
    curl -fsSL -o embryophyta_mito.fam ${base}/embryophyta_mito.fam
    curl -fsSL -o embryophyta_pltd.fam ${base}/embryophyta_pltd.fam
    """
}

process PRESS_OATKDB {
    label     'fetch'
    storeDir  "${params.outdir}/databases/oatkdb/${params.oatkdb_version}/pressed"
    container 'docker://assteindorff/oatk:1.0'

    input:
    path mito_fam
    path pltd_fam

    output:
    tuple path("embryophyta_mito.fam"), path("embryophyta_mito.fam.h3*"), emit: mito
    tuple path("embryophyta_pltd.fam"), path("embryophyta_pltd.fam.h3*"), emit: pltd

    script:
    """
    hmmpress ${mito_fam}
    hmmpress ${pltd_fam}
    """
}

// Wrapper sub-workflow so main.nf still calls one thing
workflow FETCH_OATKDB {
    main:
    DOWNLOAD_OATKDB()
    PRESS_OATKDB(DOWNLOAD_OATKDB.out.mito, DOWNLOAD_OATKDB.out.pltd)

    emit:
    mito = PRESS_OATKDB.out.mito
    pltd = PRESS_OATKDB.out.pltd
}

// ============================================================
// Build OATK HMM databases from user-supplied FASTA sequences
// (chloroplast and mitochondrion complete genome assemblies).
// Emits the same channel format as FETCH_OATKDB so the rest of
// the pipeline is unchanged.
// ============================================================

process BUILD_OATKDB_FROM_FASTA {
    label     'fetch'
    storeDir  "${params.outdir}/databases/custom_oatkdb"
    container 'docker://assteindorff/oatk:1.0'

    input:
    tuple val(prefix), path(fasta)

    output:
    tuple val(prefix), path("${prefix}.fam"), path("${prefix}.fam.h3*")

    script:
    """
    hmmbuild --dna ${prefix}.fam ${fasta}
    hmmpress ${prefix}.fam
    """
}

// ============================================================
// Kraken2 PlusPFP (Standard + protozoa/fungi/plant) — the contamination-screening DB.
// ~172 GB compressed / ~231.5 GB extracted (2026-06-26 build). storeDir caches it under
// Canopy/ itself (not per-outdir like the fetchers above) since it's a large, run-independent
// shared resource: every sample/run reuses the same copy. Only runs once, ever, unless the
// version param changes. Piping curl straight into tar avoids needing 172 GB of scratch space
// for the compressed tarball on top of the 231.5 GB extracted footprint.
//
// This exists so --run_kraken2 works without depending on the shared compbio database mount
// (no write access there for this account) — point --kraken2_db/--taxdump_dir at that shared
// copy instead, once it's provisioned, to skip this and save the local disk footprint.
// ============================================================

process DOWNLOAD_KRAKEN2_PLUSPFP {
    label    'fetch'
    storeDir "${projectDir}/databases"
    // no container — uses host curl/tar

    output:
    path "kraken2_pluspfp_${params.kraken2_pluspfp_version}", emit: db

    script:
    def dir = "kraken2_pluspfp_${params.kraken2_pluspfp_version}"
    def url = "https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_${params.kraken2_pluspfp_version}.tar.gz"
    """
    mkdir -p ${dir}
    curl -fSL ${url} | tar -xz -C ${dir}
    """
}

workflow FETCH_KRAKEN2_PLUSPFP {
    main:
    DOWNLOAD_KRAKEN2_PLUSPFP()

    emit:
    db = DOWNLOAD_KRAKEN2_PLUSPFP.out.db
}

// ============================================================
// Sylph GTDB sketch database — corroborates Kraken2's contaminant calls with an
// independent containment-ANI method against real reference genomes (GTDB r226,
// ~113k bacterial/archaeal representatives). Only the flagged candidate contigs are
// queried against this (typically a handful per sample), so a much smaller/cheaper
// database than Kraken2's PlusPFP suffices: ~17 GB vs ~172 GB compressed.
// ============================================================

process DOWNLOAD_SYLPH_GTDB {
    label    'fetch'
    storeDir "${projectDir}/databases/sylph"
    // no container — uses host curl

    output:
    path "gtdb-r226-c200-dbv1.syldb", emit: db

    script:
    """
    curl -fSL -o gtdb-r226-c200-dbv1.syldb \\
        https://faust.compbio.cs.cmu.edu/sylph-stuff/gtdb-r226-c200-dbv1.syldb
    """
}

workflow FETCH_SYLPH_GTDB {
    main:
    DOWNLOAD_SYLPH_GTDB()

    emit:
    db = DOWNLOAD_SYLPH_GTDB.out.db
}

workflow BUILD_OATKDB {
    take:
    cp_fasta   // Channel<Path>  — chloroplast reference FASTA
    mt_fasta   // Channel<Path>  — mitochondrion reference FASTA

    main:
    mito_in = mt_fasta.map { f -> tuple("custom_mito", f) }
    pltd_in = cp_fasta.map { f -> tuple("custom_pltd", f) }

    BUILD_OATKDB_FROM_FASTA(mito_in.mix(pltd_in))

    emit:
    mito = BUILD_OATKDB_FROM_FASTA.out
               .filter { label, fam, idx -> label == "custom_mito" }
               .map    { label, fam, idx -> tuple(fam, idx) }
    pltd = BUILD_OATKDB_FROM_FASTA.out
               .filter { label, fam, idx -> label == "custom_pltd" }
               .map    { label, fam, idx -> tuple(fam, idx) }
}