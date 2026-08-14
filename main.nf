#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ============================================================
// Canopy plant genome assembly pipeline — ONT
// Phase 1+2: QC → organelle split → assemble → polish → purge → BUSCO
// ============================================================

// ---- Parameters ----
params.reads          = "${projectDir}/reads"
params.cp_ref         = null
params.mt_ref         = null
params.oatkdb_version = "v20230921"
params.oatkdb_commit  = "75e8db0ac4a7d508a9a518d900876003ceb70737"
params.oatk_mito_db   = "https://raw.githubusercontent.com/c-zhou/OatkDB/${params.oatkdb_commit}/${params.oatkdb_version}/embryophyta_mito.fam"
params.oatk_pltd_db   = "https://raw.githubusercontent.com/c-zhou/OatkDB/${params.oatkdb_commit}/${params.oatkdb_version}/embryophyta_pltd.fam"
params.kraken2_pluspfp_version = "20260626"   // genome-idx/Langmead prebuilt Kraken2 index build date
params.oatkdb_recipe  = "v1-pressed"
params.outdir         = "results"
params.genome_size    = "720m"
params.busco_lineage  = "poales_odb10"
params.medaka_model   = "r1041_e82_400bps_sup_v5.0.0"
params.run_hapdup     = false
params.help           = false

params.organelle_assembler = "flye"    // "flye" (default) or "oatk"

// ---- Help message ---- (top-level function declaration — this is allowed)
def helpMessage() {
    log.info """
    ===================================================================
     Canopy — PLANT GENOME ASSEMBLY PIPELINE (ONT)
    ===================================================================
    Usage:
      nextflow run main.nf -profile condor --reads <dir> --cp_ref <fa> --mt_ref <fa> [options]

    Mandatory arguments:
      --reads [path]                Directory of input FASTQ (one subdir per sample)
      --cp_ref [path]               Chloroplast reference FASTA
      --mt_ref [path]               Mitochondrion reference FASTA

    Organelle assembly options:
      --organelle_assembler [str]   Organelle assembler: 'flye' or 'oatk'           [default: flye]
                                    'oatk' uses the embryophyta OatkDB (v20230921).
      --filter_organelles [bool]    Apply reference-based filter post-assembly      [default: true]
      --organelle_min_qcov [float]  Min fraction of contig covered by ref alignment [default: 0.5]
      --organelle_min_ident [float] Min identity (matches / alignment length)       [default: 0.7]

    Nuclear assembly options:
      --genome_size [str]           Estimated nuclear genome size                   [default: 800m]
      --nuclear_ref [path]          Nuclear reference FASTA (for scaffolding)       [default: none]
      --run_hapdup [bool]           Run HapDup phasing step                         [default: false]

    Polishing & QC options:
      --medaka_model [str]          Medaka model (null = auto-detect from reads)    [default: null]
      --busco_lineage [str]         BUSCO lineage dataset                           [default: poales_odb10]

    General options:
      --outdir [path]               Output directory                                [default: results]
      --help                        Show this message and exit
    ===================================================================
    """.stripIndent()
}

// ---- Module imports ----
// Imports are grouped by source file (one include per module), and within each
// group processes are listed in the order they appear in the workflow below.
// Module order itself follows the pipeline stages: QC → mapping → assembly → polishing.

// QC: read stats, filtering, BUSCO completeness, MultiQC aggregation
include {NANOPLOT; FILTER_READS; FILTER_ORGANELLE_CONTIGS; BUSCO_NUCLEAR; MULTIQC; BANDAGE_IMAGE; QUAST_ORGANELLE; QUAST_NUCLEAR; ALIGN_FOR_QC; SORT_FOR_QC; FILTER_PRIMARY_BAM; QUALIMAP_BAMQC; BLOBTOOLS_COVERAGE} from './modules/qc.nf'
include {FETCH_OATKDB; FETCH_KRAKEN2_PLUSPFP; FETCH_SYLPH_GTDB} from './modules/dbs.nf'
include {ALIGN_TO_ORGANELLES; SORT_INDEX_BAM; EXTRACT_RAW_READSETS; DEDUP_ORGANELLE_READS; READSET_STATS} from './modules/mapping.nf'
include {ASSEMBLE_CP_FLYE; ASSEMBLE_MT_FLYE; ASSEMBLE_ORGANELLES_OATK; ASSEMBLE_NUCLEAR} from './modules/assembly.nf'
include {POLISH_MEDAKA; POLISH_MEDAKA_ORGANELLE; PURGE_DUPS; ALIGN_FOR_HAPDUP; SORT_FOR_HAPDUP; HAPDUP} from './modules/polishing.nf'
include {RAGTAG_SCAFFOLD; RAGTAG_SCAFFOLD as RAGTAG_PREPURGE} from './modules/scaffolding.nf'
include {KRAKEN2_CLASSIFY; BLOBTOOLS_TAXONOMY; CLASSIFY_CONTAMINANTS; REMOVE_CONTAMINANTS; EXTRACT_CANDIDATE_CONTIGS; SYLPH_VERIFY_CONTAMINANTS; BLAST_VERIFY_CONTAMINANTS} from './modules/contamination.nf'
include {FINAL_SUMMARY; TOOLS_REPORT; PACKAGE_RESULTS} from './modules/reports.nf'

// ============================================================
//  Workflow
// ============================================================
workflow {

    if (params.help) {
        helpMessage()
        return            // exits the workflow cleanly — no `exit 0` needed
    }

    // ---- Validate --final_assembly (and retire --skip_purge) ----
    // purge_dups always runs now; --skip_purge has been replaced by --final_assembly, which
    // selects the published genome (params is immutable at runtime, so we don't silently remap).
    if (params.skip_purge != null) {
        def repl = (params.skip_purge.toString().toLowerCase() == 'true') ? 'medaka' : 'purge'
        exit 1, "ERROR: --skip_purge is retired. Use --final_assembly ${repl} " +
                "(purge_dups always runs; --final_assembly selects medaka|purge as the final)."
    }
    if (!(params.final_assembly in ['medaka', 'purge'])) {
        exit 1, "ERROR: --final_assembly must be 'medaka' or 'purge' (got '${params.final_assembly}')"
    }
    if (params.flag_contaminants && !(params.run_kraken2 && params.nuclear_ref)) {
        exit 1, "ERROR: --flag_contaminants requires --run_kraken2 (for BlobTools phylum calls) " +
                "and --nuclear_ref (for RagTag placement/AGP). Enable both, or drop --flag_contaminants."
    }
    if (params.verify_sylph && !params.flag_contaminants) {
        exit 1, "ERROR: --verify_sylph requires --flag_contaminants (it corroborates the " +
                "contigs that step already flagged). Enable both, or drop --verify_sylph."
    }
    if (params.verify_blast && !params.flag_contaminants) {
        exit 1, "ERROR: --verify_blast requires --flag_contaminants (it corroborates the " +
                "contigs that step already flagged). Enable both, or drop --verify_blast."
    }
    if (params.verify_blast && !params.blast_db && !params.blast_remote) {
        exit 1, "ERROR: --verify_blast requires --blast_db or --blast_remote true."
    }

    // ---- Banner ----
    log.info """
        ╔══════════════════════════════════════════════════════╗
        ║                        Canopy                        ║
        ║             PLANT GENOME ASSEMBLY PIPELINE           ║
        ║  Chloroplast · Mitochondria · Nuclear separation     ║
        ╚══════════════════════════════════════════════════════╝
        reads dir           : ${params.reads}
        cp reference        : ${params.cp_ref}
        mt reference        : ${params.mt_ref}
        nuclear ref         : ${params.nuclear_ref ?: '(none — scaffolding off)'}
        organelle assembler : ${params.organelle_assembler}${params.organelle_assembler == 'oatk' ? "  (OatkDB ${params.oatkdb_version})" : ''}
        filter organelles   : ${params.filter_organelles}
        output dir          : ${params.outdir}
        genome size         : ${params.genome_size}
        busco lineage       : ${params.busco_lineage}
        medaka model        : ${params.medaka_model}
        run HapDup          : ${params.run_hapdup}
        final assembly      : ${params.final_assembly}  (purge_dups always runs; both genomes compared)
        """.stripIndent()

    // ---- Reference channels ----
    cp_ref_ch = Channel.value(file(params.cp_ref, checkIfExists: true))
    mt_ref_ch = Channel.value(file(params.mt_ref, checkIfExists: true))
    // Nuclear reference is optional — only needed for RagTag scaffolding.
    nuclear_ref_ch = params.nuclear_ref \
        ? Channel.value(file(params.nuclear_ref, checkIfExists: true))
        : Channel.empty()

    // ---- Sample channel ----
    // Accept either one biological sample directory (reads/SAMPLE/) or a parent
    // directory containing multiple complete biological samples (reads/).
    raw_reads_ch = Channel
        .fromPath("${params.reads}/*.fastq.gz")
        .mix(Channel.fromPath("${params.reads}/*/*.fastq.gz"))
        .map { f -> tuple(f.parent.name, f) }
        .groupTuple()
        .ifEmpty { error "No FASTQ files found under --reads ${params.reads}" }

   // 1. QC + filtering
    NANOPLOT(raw_reads_ch)
    FILTER_READS(raw_reads_ch)

    // 2a. Align filtered reads to combined cp+mt reference (minimap2)
    ALIGN_TO_ORGANELLES(FILTER_READS.out.reads, cp_ref_ch, mt_ref_ch)

    // 2b. Sort and index → BAM (samtools)
    SORT_INDEX_BAM(ALIGN_TO_ORGANELLES.out.sam)

    // 3a. Extract per-compartment FASTQ (raw cp/mt may have dup IDs from supplementary alignments)
    EXTRACT_RAW_READSETS(SORT_INDEX_BAM.out.bam)

    // 3b. Dedup organelle FASTQs by read name
    dedup_in = EXTRACT_RAW_READSETS.out.cp_raw
        .join(EXTRACT_RAW_READSETS.out.mt_raw)
    DEDUP_ORGANELLE_READS(dedup_in)

    // 3c. Per-compartment stats and estimated coverage
    readset_stats_in = DEDUP_ORGANELLE_READS.out.cp_reads
        .join(DEDUP_ORGANELLE_READS.out.mt_reads)
        .join(EXTRACT_RAW_READSETS.out.nuclear_reads)
        .combine(cp_ref_ch)
        .combine(mt_ref_ch)
    READSET_STATS(readset_stats_in)

   // 4. Organelle assemblies — branch on assembler choice
    if (params.organelle_assembler == "oatk") {
        // Always use the pre-built embryophyta OatkDB.
        // Building custom HMMs from whole-genome FASTAs (hmmbuild on 140-468 kb
        // sequences) creates profiles too large for HMMER's DP matrix — integer
        // overflow at scan time. The standard DB covers sorghum organelle genes.
        FETCH_OATKDB()
        oatk_mito_db = FETCH_OATKDB.out.mito
        oatk_pltd_db = FETCH_OATKDB.out.pltd
        ASSEMBLE_ORGANELLES_OATK(
            FILTER_READS.out.reads,
            oatk_mito_db,
            oatk_pltd_db
        )
        cp_raw_ch = ASSEMBLE_ORGANELLES_OATK.out.cp_assembly
        mt_raw_ch = ASSEMBLE_ORGANELLES_OATK.out.mt_assembly

        // Visualize assembly graphs with Bandage (published next to filtered assemblies)
        bandage_in = ASSEMBLE_ORGANELLES_OATK.out.cp_gfa
            .map { id, gfa -> tuple(id, 'chloroplast', gfa) }
            .mix(
                ASSEMBLE_ORGANELLES_OATK.out.mt_gfa
                    .map { id, gfa -> tuple(id, 'mitochondria', gfa) }
            )
        BANDAGE_IMAGE(bandage_in)
    } else {
        ASSEMBLE_CP_FLYE(DEDUP_ORGANELLE_READS.out.cp_reads)
        ASSEMBLE_MT_FLYE(DEDUP_ORGANELLE_READS.out.mt_reads)
        cp_raw_ch = ASSEMBLE_CP_FLYE.out.assembly
        mt_raw_ch = ASSEMBLE_MT_FLYE.out.assembly
    }

    // 4b. Reference-based filtering — remove nuclear contamination (NUPTs/NUMTs)
    if (params.filter_organelles) {
        cp_filter_in = cp_raw_ch.map { id, fa -> tuple(id, 'chloroplast', fa) }
                                .combine(cp_ref_ch)
        mt_filter_in = mt_raw_ch.map { id, fa -> tuple(id, 'mitochondria', fa) }
                                .combine(mt_ref_ch)

        all_organelle_in = cp_filter_in.mix(mt_filter_in)
        FILTER_ORGANELLE_CONTIGS(all_organelle_in)

        // Final outputs (filtered)
        cp_final = FILTER_ORGANELLE_CONTIGS.out.filtered
            .filter { id, comp, fa -> comp == 'chloroplast' }
            .map { id, comp, fa -> tuple(id, fa) }
        mt_final = FILTER_ORGANELLE_CONTIGS.out.filtered
            .filter { id, comp, fa -> comp == 'mitochondria' }
            .map { id, comp, fa -> tuple(id, fa) }
    } else {
        cp_final = cp_raw_ch
        mt_final = mt_raw_ch
    }

    // 4c. Polish organelle assemblies with Medaka — corrects ONT base-call error before QC.
    //     Previously only the nuclear assembly was polished; cp/mt genomes went straight from
    //     assembly (+ optional reference filter) into QUAST, carrying raw ONT mismatches/indels.
    //     Uses the deduped per-compartment reads from step 3b (available regardless of which
    //     organelle assembler ran).
    organelle_polish_in = cp_final
        .join(DEDUP_ORGANELLE_READS.out.cp_reads)
        .map { id, fa, fq -> tuple(id, 'chloroplast', fa, fq) }
        .mix(
            mt_final
                .join(DEDUP_ORGANELLE_READS.out.mt_reads)
                .map { id, fa, fq -> tuple(id, 'mitochondria', fa, fq) }
        )
    POLISH_MEDAKA_ORGANELLE(organelle_polish_in)
    cp_final = POLISH_MEDAKA_ORGANELLE.out.assembly
        .filter { id, comp, fa -> comp == 'chloroplast' }
        .map    { id, comp, fa -> tuple(id, fa) }
    mt_final = POLISH_MEDAKA_ORGANELLE.out.assembly
        .filter { id, comp, fa -> comp == 'mitochondria' }
        .map    { id, comp, fa -> tuple(id, fa) }

    // 4d. QUAST on filtered+polished organelle assemblies (uses cp/mt references)
    quast_organelle_in = cp_final
        .map { id, fa -> tuple(id, 'chloroplast', fa) }
        .combine(cp_ref_ch)
        .mix(
            mt_final
                .map { id, fa -> tuple(id, 'mitochondria', fa) }
                .combine(mt_ref_ch)
        )
    QUAST_ORGANELLE(quast_organelle_in)

    // 5. Nuclear assembly (no dedup needed — unmapped reads only appear once)
    ASSEMBLE_NUCLEAR(EXTRACT_RAW_READSETS.out.nuclear_reads)

    // 6. Polish nuclear assembly with Medaka
    polish_in = ASSEMBLE_NUCLEAR.out.assembly
        .join(EXTRACT_RAW_READSETS.out.nuclear_reads)
    POLISH_MEDAKA(polish_in)

    // 7. Purge haplotigs — ALWAYS runs. Both the purged and the Medaka genomes are scaffolded
    //    and compared (QUAST + BUSCO); --final_assembly selects which becomes the published final.
    purge_in = POLISH_MEDAKA.out.assembly
        .join(EXTRACT_RAW_READSETS.out.nuclear_reads)
    PURGE_DUPS(purge_in)
    purge_final      = PURGE_DUPS.out.assembly
    purge_cutoffs_ch = PURGE_DUPS.out.cutoffs
    purge_calcuts_ch = PURGE_DUPS.out.calcuts_log

    // 8. Optional: HapDup phasing
    if (params.run_hapdup) {
        // Align reads to purge-stage assembly
        align_in = purge_final
            .join(EXTRACT_RAW_READSETS.out.nuclear_reads)
        ALIGN_FOR_HAPDUP(align_in)

        // Sort and index BAM
        SORT_FOR_HAPDUP(ALIGN_FOR_HAPDUP.out.sam)

        // Run HapDup with the prepared BAM
        HAPDUP(SORT_FOR_HAPDUP.out.bam)
    }

    // 8b. Chromosome scaffolding (RagTag correct + scaffold) when a nuclear reference is given.
    //     purge_dups always runs, so BOTH the purged and the Medaka genomes are scaffolded and
    //     compared; --final_assembly selects which scaffold becomes the published final.
    //     NB: QUAST_NUCLEAR itself is invoked further down (step 9d), after the decontam
    //     genome (if any) is known — QUAST needs it as an optional 6th column.
    if (params.nuclear_ref) {
        RAGTAG_SCAFFOLD(purge_final.map { id, fa -> tuple(id, 'purge', fa) }, nuclear_ref_ch)
        RAGTAG_PREPURGE(POLISH_MEDAKA.out.assembly.map { id, fa -> tuple(id, 'medaka', fa) }, nuclear_ref_ch)
        purge_candidate  = RAGTAG_SCAFFOLD.out.scaffold
        medaka_candidate = RAGTAG_PREPURGE.out.scaffold
        // RagTag stats/AGP for the selected final genome (drives the FINAL_SUMMARY scaffold
        // section, the unplaced-only contamination screen, and the results package).
        final_ragtag_ch  = (params.final_assembly == 'purge') ? RAGTAG_SCAFFOLD.out.stats
                                                              : RAGTAG_PREPURGE.out.stats
        final_agp_ch     = (params.final_assembly == 'purge') ? RAGTAG_SCAFFOLD.out.agp
                                                              : RAGTAG_PREPURGE.out.agp
    } else {
        // No reference → no scaffolding; the candidates are the unscaffolded genomes.
        purge_candidate  = purge_final
        medaka_candidate = POLISH_MEDAKA.out.assembly
    }

    // Select the published final genome (default: Medaka).
    nuclear_final = (params.final_assembly == 'purge') ? purge_candidate : medaka_candidate

    // 9b. Final-only QC: ONE read-to-assembly BAM for the chosen final genome, shared by
    //     Qualimap / BlobTools / Kraken2 (the discarded candidate gets no BAM-level QC).
    //     Moved ahead of BUSCO/QUAST so the optional decontam genome (produced here, when
    //     --flag_contaminants is on) is available as a channel before those are invoked.
    decontam_final = Channel.empty()
    def has_decontam = params.flag_contaminants
    if (params.run_qualimap || params.run_blobtools || params.run_kraken2) {
        qc_align_in = nuclear_final
            .map { sid, fa -> tuple(sid, params.final_assembly, fa) }
            .combine(EXTRACT_RAW_READSETS.out.nuclear_reads, by: 0)
        ALIGN_FOR_QC(qc_align_in)
        SORT_FOR_QC(ALIGN_FOR_QC.out.sam)

        if (params.run_qualimap) QUALIMAP_BAMQC(SORT_FOR_QC.out.bam)

        // BlobTools' own "% reads mapped" is unreliable on a BAM with many secondary/
        // supplementary alignments (routine for ONT + minimap2 on a repetitive genome) — see
        // FILTER_PRIMARY_BAM in modules/qc.nf for why. Feed BlobTools/Kraken2 a primary-only
        // BAM instead; Qualimap keeps the original (it benefits from seeing multi-mapping).
        if (params.run_blobtools || params.run_kraken2) {
            FILTER_PRIMARY_BAM(SORT_FOR_QC.out.bam)
            primary_bam_ch = FILTER_PRIMARY_BAM.out.bam
        }
        if (params.run_blobtools) BLOBTOOLS_COVERAGE(primary_bam_ch)

        // Contaminant screening on the final genome: Kraken2 classifies contigs against
        // PlusPFP (includes a plant clade, so host sequence gets a real match), then BlobTools
        // renders the blob plot coloured by taxonomy (non-Viridiplantae = candidate contam).
        if (params.run_kraken2) {
            // Falls back to a self-fetched local copy under Canopy/databases/ (see
            // modules/dbs.nf) when the shared compbio path isn't provisioned yet. PlusPFP
            // ships its own taxdump alongside the hash, so one fetched dir serves both roles.
            if (file(params.kraken2_db).exists() && file(params.taxdump_dir).exists()) {
                kraken2_db_ch = Channel.value(file(params.kraken2_db,  checkIfExists: true))
                taxdump_ch    = Channel.value(file(params.taxdump_dir, checkIfExists: true))
            } else {
                FETCH_KRAKEN2_PLUSPFP()
                kraken2_db_ch = FETCH_KRAKEN2_PLUSPFP.out.db
                taxdump_ch    = FETCH_KRAKEN2_PLUSPFP.out.db
            }
            KRAKEN2_CLASSIFY(primary_bam_ch, kraken2_db_ch)
            BLOBTOOLS_TAXONOMY(KRAKEN2_CLASSIFY.out.hits, taxdump_ch)

            // Flag + remove contaminant contigs from the final genome (3-way AND: non-target
            // phylum, unplaced by RagTag, abnormal GC — see modules/contamination.nf).
            if (params.flag_contaminants) {
                classify_in = BLOBTOOLS_TAXONOMY.out.results.join(final_agp_ch)
                CLASSIFY_CONTAMINANTS(classify_in)
                REMOVE_CONTAMINANTS(nuclear_final.join(CLASSIFY_CONTAMINANTS.out.ids))
                decontam_final = REMOVE_CONTAMINANTS.out.assembly

                // Independent second opinion on the flagged candidates only (typically a
                // handful of contigs) — containment ANI against real GTDB genomes, rather than
                // Kraken2's k-mer LCA approach. Advisory: does not change what got removed above.
                if (params.verify_sylph) {
                    if (file(params.sylph_db).exists()) {
                        sylph_db_ch = Channel.value(file(params.sylph_db, checkIfExists: true))
                    } else {
                        FETCH_SYLPH_GTDB()
                        sylph_db_ch = FETCH_SYLPH_GTDB.out.db
                    }
                    EXTRACT_CANDIDATE_CONTIGS(nuclear_final.join(CLASSIFY_CONTAMINANTS.out.ids))
                    SYLPH_VERIFY_CONTAMINANTS(EXTRACT_CANDIDATE_CONTIGS.out.fasta, sylph_db_ch)
                }
                if (params.verify_blast) {
                    blast_db_ch = Channel.value(params.blast_db ?: '')
                    BLAST_VERIFY_CONTAMINANTS(EXTRACT_CANDIDATE_CONTIGS.out.fasta, blast_db_ch)
                }
            }
        }
    }

    // 9c. BUSCO on BOTH candidates (stage-tagged so qc/busco/ dirs don't collide), plus the
    //     decontam genome when contaminant flagging ran. The comparison feeds the report; the
    //     final genome's BUSCO drives the summary verdict.
    busco_in = medaka_candidate.map { id, fa -> tuple(id, 'medaka', fa) }
        .mix( purge_candidate.map { id, fa -> tuple(id, 'purge', fa) } )
    if (has_decontam) {
        busco_in = busco_in.mix( decontam_final.map { id, fa -> tuple(id, 'decontam', fa) } )
    }
    BUSCO_NUCLEAR(busco_in)

    // 9d. QUAST 6-way (when decontam ran): flye / medaka / medaka_scaf / purge / purge_scaffold /
    //     decontam vs reference. decontam_join_ch pads with purge_final (always present) when
    //     contaminant flagging didn't run, same idiom as the no-ref branch below.
    decontam_join_ch = has_decontam ? decontam_final : purge_final
    if (params.nuclear_ref) {
        quast_nuclear_in = ASSEMBLE_NUCLEAR.out.assembly
            .join(POLISH_MEDAKA.out.assembly)
            .join(medaka_candidate)
            .join(purge_final)
            .join(purge_candidate)
            .join(decontam_join_ch)
            .map { id, flye, medaka, mscaf, purge, scaffold, decontam ->
                tuple(id, flye, medaka, mscaf, purge, scaffold, decontam)
            }
        QUAST_NUCLEAR(quast_nuclear_in, Channel.value(true), Channel.value(true),
                      Channel.value(has_decontam), nuclear_ref_ch)
    } else {
        quast_nuclear_in = ASSEMBLE_NUCLEAR.out.assembly
            .join(POLISH_MEDAKA.out.assembly)
            .join(purge_final)
            .map { id, flye, medaka, purge ->
                tuple(id, flye, medaka, purge, purge, purge, purge)  // 4 unused slots padded
            }
        QUAST_NUCLEAR(quast_nuclear_in, Channel.value(false), Channel.value(false),
                      Channel.value(false), Channel.value([]))
    }

    // 10. MultiQC aggregation
    // QUAST reports are staged as their whole output directory (uniquely named
    // per compartment) rather than the bare report.tsv — three files all named
    // report.tsv collide when staged flat. MultiQC recurses into each dir.
    multiqc_in = Channel.empty()
        .mix( NANOPLOT.out.report.map                        { sample, f -> f } )
        .mix( READSET_STATS.out.mqc.map                      { sample, f -> f } )
        .mix( BUSCO_NUCLEAR.out.summary.map                  { sample, stage, f -> f } )
        .mix( QUAST_ORGANELLE.out.report.map                 { sample, comp, d -> d } )
        .mix( QUAST_NUCLEAR.out.report.map                   { sample, d -> d } )
        .mix( PURGE_DUPS.out.pbstat.map                      { sample, f -> f } )
        .mix( PURGE_DUPS.out.cutoffs.map                     { sample, f -> f } )

    if (params.run_qualimap) {
        multiqc_in = multiqc_in.mix(
            QUALIMAP_BAMQC.out.report.map { sid, stage, d -> d }
        )
    }

    MULTIQC(multiqc_in.collect())

    // 11. Human-readable per-sample summary (requires a reference for RagTag scaffold stats).
    //     BUSCO ran on both candidates; pass the final genome's summary plus both candidates'
    //     summaries so the report can show the Medaka-vs-purge comparison.
    if (params.nuclear_ref) {
        nano_stats_ch = NANOPLOT.out.report
            .map { id, files ->
                def list = files instanceof List ? files : [files]
                tuple(id, list.find { it.name == 'NanoStats.txt' })
            }
        busco_summ      = BUSCO_NUCLEAR.out.summary   // (id, stage, file)
        medaka_busco_ch = busco_summ.filter { id, stage, f -> stage == 'medaka' }.map { id, stage, f -> tuple(id, f) }
        purge_busco_ch  = busco_summ.filter { id, stage, f -> stage == 'purge'  }.map { id, stage, f -> tuple(id, f) }
        // Padded with medaka_busco_ch (an always-present file) when contaminant flagging
        // didn't run — same idiom used elsewhere to fill an unused-but-required slot.
        decontam_busco_ch = has_decontam
            ? busco_summ.filter { id, stage, f -> stage == 'decontam' }.map { id, stage, f -> tuple(id, f) }
            : medaka_busco_ch

        // contamination summary is optional (only when --flag_contaminants); pad with an
        // empty file via remainder:true + null-check, same pattern as blob_ch below.
        contam_summary_ch = params.flag_contaminants ? CLASSIFY_CONTAMINANTS.out.summary : Channel.empty()
        sylph_summary_ch  = params.verify_sylph      ? SYLPH_VERIFY_CONTAMINANTS.out.summary : Channel.empty()
        blast_summary_ch  = params.verify_blast      ? BLAST_VERIFY_CONTAMINANTS.out.summary : Channel.empty()

        // FINAL_SUMMARY picks the final genome's BUSCO from medaka/purge in-script.
        summary_in = nano_stats_ch
            .join(QUAST_NUCLEAR.out.report)
            .join(QUAST_ORGANELLE.out.report
                .filter { id, comp, d -> comp == 'chloroplast' }
                .map { id, comp, d -> tuple(id, d) })
            .join(QUAST_ORGANELLE.out.report
                .filter { id, comp, d -> comp == 'mitochondria' }
                .map { id, comp, d -> tuple(id, d) })
            .join(purge_cutoffs_ch)
            .join(purge_calcuts_ch)
            .join(final_ragtag_ch)
            .join(medaka_busco_ch)
            .join(purge_busco_ch)
            .join(decontam_busco_ch)
            .join(contam_summary_ch, remainder: true)
            .join(sylph_summary_ch, remainder: true)
            .join(blast_summary_ch, remainder: true)
            .map { id, nano, quast, qcp, qmt, cutoffs, calcuts, ragtag, bmed, bpurge, bdecon, contam, sylph, blast ->
                tuple(id, nano, quast, cutoffs, calcuts, ragtag, qcp, qmt, bmed, bpurge, bdecon,
                      contam != null ? contam : [],
                      sylph  != null ? sylph  : [],
                      blast  != null ? blast  : [])
            }
        FINAL_SUMMARY(summary_in)
    }
    TOOLS_REPORT()

    // 12. Package the final assembly + reports as one zip per sample (see modules/reports.nf)
    //     so a run's key results can be grabbed over FileZilla without touching the full,
    //     much larger (100+ GB) output directory.
    if (params.nuclear_ref) {
        final_busco_ch = busco_summ.filter { id, stage, f -> stage == params.final_assembly }
                                    .map    { id, stage, f -> tuple(id, f) }

        quast_cp_ch = QUAST_ORGANELLE.out.report
            .filter { id, comp, d -> comp == 'chloroplast' }
            .map    { id, comp, d -> tuple(id, d) }
        quast_mt_ch = QUAST_ORGANELLE.out.report
            .filter { id, comp, d -> comp == 'mitochondria' }
            .map    { id, comp, d -> tuple(id, d) }

        // Blob plot is optional (only when --run_kraken2); `remainder: true` lets samples
        // through with a null placeholder when the screening step didn't run.
        blob_ch = params.run_kraken2 ? BLOBTOOLS_TAXONOMY.out.results : Channel.empty()

        // Decontam FASTA + audit table are optional (only when --flag_contaminants); same
        // remainder:true + null-check idiom as blob_ch. Sylph corroboration likewise optional.
        decontam_pkg_ch = has_decontam ? decontam_final : Channel.empty()
        contam_audit_ch = has_decontam ? CLASSIFY_CONTAMINANTS.out.audit : Channel.empty()
        sylph_pkg_ch    = params.verify_sylph ? SYLPH_VERIFY_CONTAMINANTS.out.corroboration : Channel.empty()
        blast_pkg_ch    = params.verify_blast ? BLAST_VERIFY_CONTAMINANTS.out.corroboration : Channel.empty()

        package_in = FINAL_SUMMARY.out.summary
            .join(nuclear_final)
            .join(final_agp_ch)
            .join(cp_final)
            .join(mt_final)
            .join(QUAST_NUCLEAR.out.report)
            .join(quast_cp_ch)
            .join(quast_mt_ch)
            .join(final_busco_ch)
            .join(blob_ch, remainder: true)
            .join(decontam_pkg_ch, remainder: true)
            .join(contam_audit_ch, remainder: true)
            .join(sylph_pkg_ch, remainder: true)
            .join(blast_pkg_ch, remainder: true)
            .combine(TOOLS_REPORT.out)
            .map { id, summary, scaffold, agp, cp, mt, qn, qcp, qmt, busco, blob, decontam, audit, sylph, blast, tools ->
                def hasBlob     = (blob != null)
                def hasDecontam = (decontam != null)
                def hasSylph    = (sylph != null)
                def hasBlast    = (blast != null)
                tuple(id, summary, tools, scaffold, agp, cp, mt, qn, qcp, qmt, busco,
                      hasBlob ? blob : [], hasBlob,
                      hasDecontam ? decontam : [], hasDecontam ? audit : [], hasDecontam,
                      hasSylph ? sylph : [], hasSylph,
                      hasBlast ? blast : [], hasBlast)
            }
        PACKAGE_RESULTS(package_in)
    }

    // Capture outdir into a local — `params` resolves to null inside the
    // onComplete closure when it fires, so reference the captured value.
    def outdir = params.outdir ?: 'NA'
    workflow.onComplete {
        log.info "Pipeline completed | Output: ${outdir}"
    }
}
