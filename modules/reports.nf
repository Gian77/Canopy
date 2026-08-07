// Human-readable final summary.
// Pulls the key metrics from NanoPlot (read set), QUAST (nuclear assembly), and
// BUSCO (completeness) into one markdown table per sample. Reporting only —
// errorStrategy 'ignore' so a missing/odd input never fails the pipeline.
// Reuses the already-pulled multiqc image (has bash/grep/awk/python).

process FINAL_SUMMARY {
    tag           { sample_id }
    label         'qc'
    errorStrategy 'ignore'
    publishDir    { "${params.outdir}/reports/${sample_id}" }, mode: 'symlink'
    container     'quay.io/biocontainers/multiqc:1.25.1--pyhdfd78af_0'

    input:
    // BUSCO ran on both candidates (plus decontam when --flag_contaminants); the script picks
    // the final genome's summary from these (staged under distinct names to avoid a collision
    // when final == medaka). busco_decontam/contam_summary are padded/empty when contaminant
    // flagging didn't run — see main.nf's decontam_busco_ch/contam_summary_ch.
    tuple val(sample_id), path(nanostats), path(quast_nuclear),
          path(purge_cutoffs), path(purge_calcuts_log), path(ragtag_stats),
          path(busco_medaka,    stageAs: 'busco_medaka_summary.txt'),
          path(busco_purge,     stageAs: 'busco_purge_summary.txt'),
          path(busco_decontam,  stageAs: 'busco_decontam_summary.txt'),
          path(contam_summary,  stageAs: 'contam_summary.txt'),
          path(sylph_summary,   stageAs: 'sylph_summary.txt')

    output:
    tuple val(sample_id), path("${sample_id}_assembly_summary.md"), emit: summary

    script:
    def busco_lin = params.busco_lineage
    def sample    = sample_id
    def gsize     = params.genome_size
    def final_asm = params.final_assembly
    def screened  = (params.run_kraken2 || params.run_blobtools)
    def taxo_screened = params.run_kraken2
    def decontam_ran  = params.flag_contaminants
    def sylph_ran     = params.verify_sylph
    """
    #!/usr/bin/env bash
    set -uo pipefail

    NANO="${nanostats}"
    QUAST="${quast_nuclear}/report.tsv"
    BUSCO_MEDAKA="busco_medaka_summary.txt"
    BUSCO_PURGE="busco_purge_summary.txt"
    BUSCO_DECONTAM="busco_decontam_summary.txt"
    CONTAM_SUMMARY="contam_summary.txt"
    SYLPH_SUMMARY="sylph_summary.txt"
    # Headline BUSCO = the published final genome's (default medaka).
    if [ "${final_asm}" = "purge" ]; then BUSCO="\$BUSCO_PURGE"; else BUSCO="\$BUSCO_MEDAKA"; fi
    CUTOFFS="${purge_cutoffs}"
    CALCUTS="${purge_calcuts_log}"
    RAGTAG="${ragtag_stats}"

    # ── helper: NanoPlot field lookup ───────────────────────────────────────
    nano() { awk -F'\\t' -v k="\$1" 'tolower(\$1) ~ tolower(k) {print \$NF; exit}' "\$NANO" 2>/dev/null || echo NA; }

    # ── helper: QUAST row as markdown ───────────────────────────────────────
    qrow() {
        awk -F'\\t' -v k="\$1" '
            NR==1 { ncols=NF }
            \$1==k { printf "| **%s** |", \$1; for(i=2;i<=ncols;i++) printf " %s |", \$i; print ""; exit }
        ' "\$QUAST" 2>/dev/null
    }

    qheader() {
        awk -F'\\t' '
            NR==1 { printf "| Metric |"; for(i=2;i<=NF;i++) printf " **%s** |", \$i; print "";
                    printf "|---|";       for(i=2;i<=NF;i++) printf "---|";           print "" }
        ' "\$QUAST" 2>/dev/null
    }

    # ── helper: extract single QUAST cell by metric + column-header pattern ─
    qval() {
        local metric="\$1" colpat="\$2"
        awk -F'\\t' -v m="\$metric" -v p="\$colpat" '
            NR==1 { for(i=2;i<=NF;i++) if(tolower(\$i) ~ tolower(p)) { col=i; break } }
            \$1==m && col { print \$col; exit }
        ' "\$QUAST" 2>/dev/null || echo 0
    }

    # ── extract metrics ──────────────────────────────────────────────────────
    TOTAL_BASES=\$(nano 'number.of.bases')
    READ_N50=\$(nano 'n50')
    MEAN_QUAL=\$(nano 'mean.qual')

    GFRAC_FLYE=\$(qval "Genome fraction (%)" "flye")
    GFRAC_PURGE=\$(qval "Genome fraction (%)" "purge_dups")
    GFRAC_SCAFFOLD=\$(qval "Genome fraction (%)" "purge_dups_scaf")
    # Pre-purge (Medaka) RagTag scaffold, when present.
    GFRAC_MEDAKA_SCAF=\$(qval "Genome fraction (%)" "medaka_scaf")
    GFRAC_DECONTAM=\$(qval "Genome fraction (%)" "decontam")
    N50_DECONTAM=\$(qval "N50" "decontam")
    LEN_DECONTAM=\$(qval "Total length" "decontam")
    N50_FLYE=\$(qval "N50" "flye")
    N50_SCAFFOLD=\$(qval "N50" "purge_dups_scaf")
    LEN_FLYE=\$(qval "Total length" "flye")
    LEN_PURGE=\$(qval "Total length" "purge_dups")
    NCONTIGS_FLYE=\$(qval "# contigs" "flye")
    MISASM_FLYE=\$(qval "# misassemblies" "flye")
    MISASM_SCAFFOLD=\$(qval "# misassemblies" "purge_dups_scaf")

    # The published final genome is selected by --final_assembly (default: medaka). Keep the
    # purge-scaffold genome fraction for the comparison table, then repoint the headline
    # *_SCAFFOLD metrics at the final genome's column so the verdict/conclusion describe it.
    GFRAC_PURGE_SCAF="\$GFRAC_SCAFFOLD"
    if [ "${final_asm}" != "purge" ]; then
        GFRAC_SCAFFOLD="\$GFRAC_MEDAKA_SCAF"
        N50_SCAFFOLD=\$(qval "N50" "medaka_scaf")
        MISASM_SCAFFOLD=\$(qval "# misassemblies" "medaka_scaf")
    fi

    COVERAGE=\$(echo "${gsize}" | awk '{
        s=tolower(\$0)
        if (sub(/g\$/, "", s)) s=s*1e9
        else if (sub(/m\$/, "", s)) s=s*1e6
        else if (sub(/k\$/, "", s)) s=s*1e3
        print s
    }' | awk -v tb="\$TOTAL_BASES" 'BEGIN{tb+=0} {if(\$1>0) printf "%.0fx", tb/\$1; else print "NAx"}')

    PURGE_RETAINED=\$(echo "\$LEN_PURGE \$LEN_FLYE" | awk '{if(\$2>0) printf "%.1f", \$1/\$2*100; else print "NA"}')
    PURGE_REMOVED=\$(echo "\$LEN_PURGE \$LEN_FLYE"  | awk '{if(\$2>0) printf "%.1f", (1-\$1/\$2)*100; else print "NA"}')
    C_JUNK=\$(awk '{print \$1}' "\$CUTOFFS" 2>/dev/null || echo "?")
    C_HAP_LOW=\$(awk '{print \$2}' "\$CUTOFFS" 2>/dev/null || echo "?")
    C_HAP_HIGH=\$(awk '{print \$3}' "\$CUTOFFS" 2>/dev/null || echo "?")
    C_DIP_LOW=\$(awk '{print \$4}' "\$CUTOFFS" 2>/dev/null || echo "?")
    C_DIP_HIGH=\$(awk '{print \$5}' "\$CUTOFFS" 2>/dev/null || echo "?")
    C_REPEAT=\$(awk '{print \$6}' "\$CUTOFFS" 2>/dev/null || echo "?")

    PLACED_MB=\$(awk -F'\\t' 'NR==2 {printf "%.0f", \$2/1e6}' "\$RAGTAG" 2>/dev/null || echo 0)
    UNPLACED_MB=\$(awk -F'\\t' 'NR==2 {printf "%.0f", \$4/1e6}' "\$RAGTAG" 2>/dev/null || echo 0)
    PLACED_PCT=\$(echo "\$PLACED_MB \$UNPLACED_MB" | awk '{t=\$1+\$2; if(t>0) printf "%.0f", \$1/t*100; else print "NA"}')

    BUSCO_C=\$(grep 'C:[0-9]' "\$BUSCO" 2>/dev/null | grep -o 'C:[0-9][0-9.]*' | cut -d: -f2 | head -1 || echo 0)
    BUSCO_S=\$(grep 'S:[0-9]' "\$BUSCO" 2>/dev/null | grep -o 'S:[0-9][0-9.]*' | cut -d: -f2 | head -1 || echo 0)
    BUSCO_D=\$(grep 'D:[0-9]' "\$BUSCO" 2>/dev/null | grep -o 'D:[0-9][0-9.]*' | cut -d: -f2 | head -1 || echo 0)
    BUSCO_F=\$(grep 'F:[0-9]' "\$BUSCO" 2>/dev/null | grep -o 'F:[0-9][0-9.]*' | cut -d: -f2 | head -1 || echo 0)
    BUSCO_M=\$(grep 'M:[0-9]' "\$BUSCO" 2>/dev/null | grep -o 'M:[0-9][0-9.]*' | cut -d: -f2 | head -1 || echo 0)

    # BUSCO completeness/duplication for BOTH candidates (for the comparison table).
    BUSCO_C_MEDAKA=\$(grep 'C:[0-9]' "\$BUSCO_MEDAKA" 2>/dev/null | grep -o 'C:[0-9][0-9.]*' | cut -d: -f2 | head -1 || echo 0)
    BUSCO_D_MEDAKA=\$(grep 'D:[0-9]' "\$BUSCO_MEDAKA" 2>/dev/null | grep -o 'D:[0-9][0-9.]*' | cut -d: -f2 | head -1 || echo 0)
    BUSCO_C_PURGE=\$(grep 'C:[0-9]' "\$BUSCO_PURGE" 2>/dev/null | grep -o 'C:[0-9][0-9.]*' | cut -d: -f2 | head -1 || echo 0)
    BUSCO_D_PURGE=\$(grep 'D:[0-9]' "\$BUSCO_PURGE" 2>/dev/null | grep -o 'D:[0-9][0-9.]*' | cut -d: -f2 | head -1 || echo 0)
    BUSCO_C_DECONTAM=\$(grep 'C:[0-9]' "\$BUSCO_DECONTAM" 2>/dev/null | grep -o 'C:[0-9][0-9.]*' | cut -d: -f2 | head -1 || echo 0)
    BUSCO_D_DECONTAM=\$(grep 'D:[0-9]' "\$BUSCO_DECONTAM" 2>/dev/null | grep -o 'D:[0-9][0-9.]*' | cut -d: -f2 | head -1 || echo 0)

    # Contaminant screening counts (only meaningful when decontam_ran).
    CONTAM_FLAGGED=\$(grep '^contaminants_flagged=' "\$CONTAM_SUMMARY" 2>/dev/null | cut -d= -f2 || echo 0)
    CONTAM_BASES=\$(grep '^bases_removed='        "\$CONTAM_SUMMARY" 2>/dev/null | cut -d= -f2 || echo 0)
    CONTAM_PHYLA=\$(grep '^phyla_detected='       "\$CONTAM_SUMMARY" 2>/dev/null | cut -d= -f2 || echo "none")
    CONTAM_MB=\$(echo "\$CONTAM_BASES" | awk '{printf "%.2f", \$1/1e6}')

    # Sylph corroboration counts (only meaningful when sylph_ran).
    SYLPH_CORROBORATED=\$(grep '^sylph_corroborated='      "\$SYLPH_SUMMARY" 2>/dev/null | cut -d= -f2 || echo 0)
    SYLPH_TOTAL=\$(grep '^sylph_total_candidates='         "\$SYLPH_SUMMARY" 2>/dev/null | cut -d= -f2 || echo 0)

    PEAK=\$(grep "autotune" "\$CALCUTS" 2>/dev/null | grep -o 'Peak: [0-9]*x' | grep -o '[0-9]*' || echo "unknown")

    QUALITY_TIER=\$(awk -v c="\$BUSCO_C" -v gf="\$GFRAC_SCAFFOLD" -v n50="\$N50_SCAFFOLD" 'BEGIN {
        c=c+0; gf=gf+0; n50=n50+0; score=0
        if (c  >= 95) score+=2; else if (c  >= 90) score++
        if (gf >= 80) score+=2; else if (gf >= 60) score++
        if (n50 >= 20000000) score+=2; else if (n50 >= 5000000) score++
        if (score >= 5) print "high"
        else if (score >= 3) print "moderate"
        else print "low"
    }')

    # ── interpretation helpers ───────────────────────────────────────────────
    read_interp() {
        awk -v q="\$MEAN_QUAL" -v n50="\$READ_N50" -v cov="\$COVERAGE" 'BEGIN {
            q=q+0; n50=n50+0; cov=cov+0
            if (q >= 20)      qmsg="excellent (Q>=20)"
            else if (q >= 15) qmsg="good (Q15-Q20)"
            else if (q >= 12) qmsg="acceptable (Q12-Q15)"
            else              qmsg="low (Q<12) -- polishing essential"
            if (n50 >= 30000)      nmsg="excellent (>=30 kb)"
            else if (n50 >= 15000) nmsg="good (15-30 kb)"
            else if (n50 >= 8000)  nmsg="acceptable (8-15 kb)"
            else                   nmsg="short (<8 kb) -- may limit contiguity"
            if (cov >= 50)      cmsg="deep -- well-suited for polishing and purging"
            else if (cov >= 25) cmsg="adequate for ONT assembly"
            else if (cov >= 15) cmsg="moderate -- assembly may be fragmented"
            else                cmsg="low -- additional sequencing recommended"
            printf "Mean quality: %s. Read N50: %s. Coverage ~%dx (%s).", qmsg, nmsg, cov, cmsg
        }'
    }

    assembly_interp() {
        awk -v gf="\$GFRAC_SCAFFOLD" -v n50="\$N50_SCAFFOLD" -v misasm="\$MISASM_SCAFFOLD" -v mf="\$MISASM_FLYE" 'BEGIN {
            gf=gf+0; n50=n50+0; misasm=misasm+0; mf=mf+0
            if (gf >= 85)      gfmsg="high -- assembly covers most of the reference"
            else if (gf >= 70) gfmsg="moderate -- some regions absent or divergent"
            else if (gf >= 50) gfmsg="low -- significant sequence missing or highly divergent"
            else               gfmsg="very low -- highly divergent from reference, or assembly fragmented"
            if (n50 >= 50000000)      n50msg="chromosome-scale (>=50 Mb)"
            else if (n50 >= 10000000) n50msg="near-chromosome-scale (10-50 Mb)"
            else if (n50 >= 1000000)  n50msg="Mb-scale (1-10 Mb)"
            else                      n50msg="sub-Mb -- scaffolding may be limited by contig fragmentation"
            printf "Scaffold N50 is %.1f Mb. Reference genome fraction: %.1f%% (%s). Misassemblies: %d (Flye) -> %d (scaffold).", n50/1e6, gf, gfmsg, mf, misasm
        }'
    }

    purge_interp() {
        awk -v ret="\$PURGE_RETAINED" -v gf_pre="\$GFRAC_FLYE" -v gf_post="\$GFRAC_PURGE" -v peak="\$PEAK" \
            -v junk="\$C_JUNK" -v hl="\$C_HAP_LOW" -v hh="\$C_HAP_HIGH" \
            -v dl="\$C_DIP_LOW" -v dh="\$C_DIP_HIGH" -v rep="\$C_REPEAT" 'BEGIN {
            ret=ret+0; gf_pre=gf_pre+0; gf_post=gf_post+0; drop=gf_pre-gf_post
            removed=sprintf("%.1f", 100-ret)
            if (peak != "unknown") {
                pmsg = "Autotune detected haploid coverage peak at " peak "x. "
                zmsg = "Coverage zones: junk (<" junk "x), haploid-kept (" hl "-" hh "x, centered on " peak "x), haplotig-purged (" dl "-" dh "x), repeat (>" rep "x). "
            } else {
                pmsg = "Coverage peak auto-detection unavailable. "
                zmsg = ""
            }
            mech = "Removed " removed "% of assembled sequence (genome fraction: " sprintf("%.1f",gf_pre) "% -> " sprintf("%.1f",gf_post) "%, drop " sprintf("%.1f",drop) "%). "
            mech = mech "Purging is driven primarily by the self-alignment overlap step: contigs that map redundantly against the primary assembly are classified as haplotigs and discarded; the coverage cutoffs assign each contig to a depth class before overlap testing. "
            if (drop > 15)
                verdict = "Genome fraction fell " sprintf("%.1f",drop) " points -- purge_dups removed UNIQUE reference sequence, not just haplotigs (over-purging). A true haplotig purge leaves genome fraction roughly flat. The Medaka (unpurged) genome is the default final; keep --final_assembly medaka (or set manual --calcuts_args), especially at low coverage. See the comparison below."
            else if (ret < 55)
                verdict = "Large removal -- possible over-purging for a homozygous sample. Consider --calcuts_args to set manual cutoffs, or skip purge_dups."
            else if (ret < 75)
                verdict = "Moderate removal, typical for a heterozygous plant genome where alternate-haplotype contigs are assembled separately at similar depth to the primary contigs."
            else
                verdict = "Conservative removal -- sample is likely largely homozygous with minimal haplotig duplication in the assembly."
            print pmsg zmsg mech verdict
        }'
    }

    scaffold_interp() {
        awk -v pct="\$PLACED_PCT" -v n50c="\$N50_FLYE" -v n50s="\$N50_SCAFFOLD" 'BEGIN {
            pct=pct+0; n50c=n50c+0; n50s=n50s+0
            if (pct >= 85)      anch="excellent anchoring"
            else if (pct >= 70) anch="good anchoring"
            else if (pct >= 50) anch="moderate anchoring"
            else                anch="low anchoring -- many contigs unplaced, possibly from highly divergent or novel sequence"
            fold=int(n50s/n50c)
            printf "%s (%d%% of bases placed on chromosomes). N50 improved %d-fold: %.1f kb (contigs) -> %.1f Mb (scaffolds).", anch, pct, fold, n50c/1e3, n50s/1e6
        }'
    }

    busco_interp() {
        awk -v c="\$BUSCO_C" -v s="\$BUSCO_S" -v d="\$BUSCO_D" -v f="\$BUSCO_F" -v m="\$BUSCO_M" 'BEGIN {
            c=c+0; s=s+0; d=d+0; f=f+0; m=m+0
            if (c >= 95 && d <= 3)      qual="Excellent gene-space recovery"
            else if (c >= 90 && d <= 5) qual="Good gene-space recovery"
            else if (c >= 80)           qual="Moderate gene-space recovery"
            else                        qual="Poor gene-space recovery -- significant gene loss detected"
            notes=""
            if (d > 5)  notes=notes " Duplication (" d "%) elevated -- residual haplotigs may remain."
            if (m > 10) notes=notes " Missing BUSCOs (" m "%) -- check assembly completeness."
            if (f > 5)  notes=notes " Fragmented BUSCOs (" f "%) indicate assembly fragmentation."
            printf "%s (complete: %s%%, single: %s%%, duplicated: %s%%, fragmented: %s%%, missing: %s%%).", qual, c, s, d, f, m
            if (notes != "") printf " Notes:%s", notes
            printf "\\n"
        }'
    }

    {
      echo "# Assembly summary — ${sample}"
      echo
      echo "_Generated by Canopy on \$(date -u '+%Y-%m-%d %H:%M UTC')_"
      echo

      # ── 1. Read set ───────────────────────────────────────────────────────
      echo "## 1. Read set (NanoPlot)"
      echo
      echo "| Metric | Value |"
      echo "|---|---|"
      echo "| Number of reads     | \$(nano 'number.of.reads') |"
      echo "| Total bases         | \$TOTAL_BASES |"
      echo "| Estimated coverage  | \$COVERAGE |"
      echo "| Read N50            | \$READ_N50 bp |"
      echo "| Median read length  | \$(nano 'median.read.length') bp |"
      echo "| Mean read quality   | Q\$MEAN_QUAL |"
      echo
      echo "> **Interpretation:** \$(read_interp)"
      echo

      # ── 2. Nuclear assembly progression ──────────────────────────────────
      echo "## 2. Nuclear assembly progression (QUAST)"
      echo
      qheader
      qrow "# contigs"
      qrow "Total length"
      qrow "Largest contig"
      qrow "N50"
      qrow "GC (%)"
      qrow "Genome fraction (%)"
      qrow "# misassemblies"
      qrow "# N's per 100 kbp"
      echo
      echo "> **Interpretation:** \$(assembly_interp)"
      echo

      # ── 3. Purge_dups decision ────────────────────────────────────────────
      echo "## 3. Haplotig removal decision (purge_dups vs Medaka)"
      echo
      echo "Purge_dups always runs; both the Medaka and purged genomes are scaffolded and compared."
      echo "**Final assembly adopted: ${final_asm}.**"
      echo
      echo "Coverage cutoffs (junk / hap-low / hap-high / dip-low / dip-high / repeat):"
      echo
      echo '```'
      cat "\$CUTOFFS" 2>/dev/null || echo "NA"
      echo '```'
      echo
      WARN=\$(grep -i "warn" "\$CALCUTS" 2>/dev/null | head -3 || true)
      if [ -n "\$WARN" ]; then
        echo "> **Warning:** \$WARN"
        echo
      fi
      echo "> **Interpretation:** \$(purge_interp)"
      echo
      echo "| Candidate (scaffolded) | Genome fraction | BUSCO complete | BUSCO duplicated |"
      echo "|---|---|---|---|"
      echo "| Medaka (unpurged) | \${GFRAC_MEDAKA_SCAF}% | \${BUSCO_C_MEDAKA}% | \${BUSCO_D_MEDAKA}% |"
      echo "| purge_dups        | \${GFRAC_PURGE_SCAF}% | \${BUSCO_C_PURGE}% | \${BUSCO_D_PURGE}% |"
      echo
      echo "> **How to read this:** a large genome-fraction gap with little reduction in BUSCO"
      echo "> duplication means purge_dups discarded unique sequence (over-purging) rather than"
      echo "> haplotigs. Medaka is the default final; pass --final_assembly purge to adopt the purged genome."
      echo

      # ── 4. Scaffolding ────────────────────────────────────────────────────
      echo "## 4. Chromosomal scaffolding (RagTag)"
      echo
      echo "| | Sequences | Bases |"
      echo "|---|---|---|"
      awk -F'\\t' 'NR==2 {
          placed_mb   = sprintf("%.1f Mb", \$2/1e6);
          unplaced_mb = sprintf("%.1f Mb", \$4/1e6);
          gap_kb      = sprintf("%.1f kb", \$5/1e3);
          printf "| Placed on chromosomes | %s | %s |\\n", \$1, placed_mb;
          printf "| Unplaced              | %s | %s |\\n", \$3, unplaced_mb;
          printf "| N-gaps introduced     | %s | %s |\\n", \$6, gap_kb
      }' "\$RAGTAG" 2>/dev/null || echo "| NA | NA | NA |"
      echo
      echo "> **Interpretation:** \$(scaffold_interp)"
      echo

      # ── 5. BUSCO ──────────────────────────────────────────────────────────
      echo "## 5. Gene space completeness (BUSCO — ${busco_lin})"
      echo
      echo '```'
      grep -E "C:|Complete|Fragmented|Missing|Total BUSCO" "\$BUSCO" 2>/dev/null | sed 's/^[[:space:]]*//' || echo "NA"
      echo '```'
      echo
      echo "> **Interpretation:** \$(busco_interp)"
      echo

      # ── 6. Contaminant screening (only when --flag_contaminants ran) ───────
      if [ "${decontam_ran}" = "true" ]; then
          echo "## 6. Contaminant screening (BlobTools + RagTag + GC)"
          echo
          echo "| Metric | Value |"
          echo "|---|---|"
          echo "| Contigs flagged as contaminant | \$CONTAM_FLAGGED |"
          echo "| Total bases removed | \${CONTAM_MB} Mb |"
          echo "| Non-target phyla detected | \$CONTAM_PHYLA |"
          echo "| Decontam genome fraction | \${GFRAC_DECONTAM}% |"
          echo "| Decontam N50 | \$(echo \$N50_DECONTAM | awk '{printf "%.1f Mb", \$1/1e6}') |"
          echo "| Decontam BUSCO complete / duplicated | \${BUSCO_C_DECONTAM}% / \${BUSCO_D_DECONTAM}% |"
          if [ "${sylph_ran}" = "true" ]; then
              echo "| Sylph-corroborated (independent ANI vs GTDB) | \$SYLPH_CORROBORATED / \$SYLPH_TOTAL |"
          fi
          echo
          echo "> **Interpretation:** \$CONTAM_FLAGGED contigs (\${CONTAM_MB} Mb) were flagged as" \
               "contaminant -- simultaneously non-${params.contam_target_phylum} by BlobTools," \
               "unplaced by RagTag, and outside the ${params.contam_gc_min}-${params.contam_gc_max}" \
               "GC fraction typical of eukaryotic nuclear sequence. See qc/contamination/ for the" \
               "full per-contig audit table before trusting this genome version, given known" \
               "Kraken2 long-contig misclassification risk on assembled scaffolds."
          if [ "${sylph_ran}" = "true" ]; then
              echo
              echo "> Sylph corroboration: \$SYLPH_CORROBORATED of \$SYLPH_TOTAL flagged contigs" \
                   "also had a confident (>=${params.sylph_min_ani}% adjusted ANI) containment" \
                   "match to a real GTDB reference genome -- an independent method agreeing with" \
                   "Kraken2's taxonomy call. See qc/contamination/*_sylph_corroboration.tsv for" \
                   "per-contig detail. This is advisory only and did not change what was removed." \
                   "Note: short contigs (a few kb) may show 'no' simply because there isn't" \
                   "enough sequence for a confident ANI estimate at Sylph's default sketch" \
                   "density -- that is not the same as Sylph actively disagreeing with Kraken2."
          fi
          echo
      fi

      # ── Conclusion ────────────────────────────────────────────────────────
      echo "---"
      echo
      echo "## Conclusion"
      echo
      echo "The **${sample}** nuclear genome assembly is of **\${QUALITY_TIER} overall quality** based on BUSCO completeness (\${BUSCO_C}%), scaffold N50 (\${N50_SCAFFOLD} bp), and reference genome fraction (\${GFRAC_SCAFFOLD}%)."
      echo
      PURGE_MB=\$(echo \$LEN_PURGE | awk '{printf "%.0f Mb", \$1/1e6}')
      if [ "${final_asm}" = "purge" ]; then
        HAPLO_SENTENCE="**Haplotig removal:** Purge_dups removed \${PURGE_REMOVED}% of the assembly (autotune peak: \${PEAK}x), retaining \${PURGE_MB}; the purged genome was adopted as final."
      else
        HAPLO_SENTENCE="**Haplotig removal:** Purge_dups ran (would remove \${PURGE_REMOVED}%, autotune peak \${PEAK}x) but the Medaka (unpurged) genome was adopted as final -- see the section 3 comparison."
      fi
      echo "**Sequencing:** \$COVERAGE ONT reads at mean quality Q\${MEAN_QUAL} (read N50 \${READ_N50} bp). **Assembly:** Flye produced \$NCONTIGS_FLYE contigs totalling \$(echo \$LEN_FLYE | awk '{printf "%.0f Mb", \$1/1e6}') from a \$(echo \$TOTAL_BASES | awk '{printf "%.1f Gb", \$1/1e9}') read set. **Polishing:** Medaka corrected base-level errors without structural changes. \${HAPLO_SENTENCE} **Scaffolding:** RagTag anchored \${PLACED_PCT}% of assembled bases to chromosomal positions using the reference, raising scaffold N50 to \$(echo \$N50_SCAFFOLD | awk '{printf "%.1f Mb", \$1/1e6}'). **Gene space:** BUSCO completeness \${BUSCO_C}% (${busco_lin}), with \${BUSCO_D}% duplication and \${BUSCO_M}% missing genes."
      echo
      awk -v c="\$BUSCO_C" -v gf="\$GFRAC_SCAFFOLD" -v ret="\$PURGE_RETAINED" -v dup="\$BUSCO_D" -v m="\$BUSCO_M" 'BEGIN {
          c=c+0; gf=gf+0; ret=ret+0; dup=dup+0; m=m+0
          rec=""
          if (ret < 55)
              rec=rec "- Purge_dups removed " sprintf("%.1f",100-ret) "% of the assembly. If the sample is largely homozygous, try --calcuts_args to set manual cutoffs, or consider skipping purge_dups.\\n"
          if (dup > 5)
              rec=rec "- BUSCO duplication (" dup "%) is elevated -- residual haplotigs may remain. Consider stricter purging parameters.\\n"
          if (gf < 60)
              rec=rec "- Reference genome fraction (" gf "%) is low. If BUSCO is high, this likely reflects genuine divergence from reference (not poor assembly quality).\\n"
          if (c < 90)
              rec=rec "- BUSCO completeness (" c "%) is below 90%. Consider deeper sequencing, additional polishing, or verifying the busco lineage param.\\n"
          if (m > 10)
              rec=rec "- " m "% of BUSCO genes are missing. Check for contamination, coverage gaps, or lineage mismatch.\\n"
          if (rec != "")
              printf "**Suggestions:**\\n\\n%s\\n", rec
      }'
      ${taxo_screened
        ? 'echo "_Contamination screening: Kraken2 + BlobTools were run on the whole final assembly against PlusPFP (includes a plant reference, so host sequence should classify correctly); see qc/kraken2/ and qc/blobtools/ and judge non-Viridiplantae phyla on the blob plot as candidate contamination._"'
        : (screened
            ? 'echo "_Coverage-only BlobTools plot was run on the whole final assembly (qc/blobtools/) -- no taxonomy DB, just GC-vs-depth clustering. Enable --run_kraken2 for a taxonomy-based contamination screen._"'
            : 'echo "_No contamination screening was performed. Enable --run_kraken2 / --run_blobtools to screen the final assembly._"')}
    } > ${sample}_assembly_summary.md
    """
}

// Bundles the small set of files worth downloading to sanity-check a run — the narrative
// summary, final assembly + organelle FASTAs, the RagTag AGP (maps each sequence to the
// reference chromosome it was placed on, or shows it as unplaced), QUAST/BUSCO reports, and
// (when contamination screening ran) the unplaced-only blob plot — as one zip under reports/,
// instead of someone having to hunt through the full (100+ GB) output directory over FileZilla.
process PACKAGE_RESULTS {
    tag           { sample_id }
    label         'qc'
    errorStrategy 'ignore'
    publishDir    { "${params.outdir}/reports/${sample_id}" }, mode: 'copy'
    container     'quay.io/biocontainers/multiqc:1.25.1--pyhdfd78af_0'

    input:
    tuple val(sample_id),
          path(assembly_summary),
          path(tools_report),
          path(final_scaffold, stageAs: 'final_scaffold.fasta'),
          path(final_agp,      stageAs: 'final_scaffold.agp'),
          path(cp_fasta,       stageAs: 'chloroplast.fasta'),
          path(mt_fasta,       stageAs: 'mitochondria.fasta'),
          path(quast_nuclear_dir),
          path(quast_cp_dir),
          path(quast_mt_dir),
          path(busco_summary,  stageAs: 'busco_final_summary.txt'),
          path(blob_pngs),
          val(has_blobplots),
          path(decontam_fasta, stageAs: 'decontam.fasta'),
          path(contam_audit,   stageAs: 'contamination_audit.tsv'),
          val(has_decontam),
          path(sylph_corroboration, stageAs: 'sylph_corroboration.tsv'),
          val(has_sylph)

    output:
    tuple val(sample_id), path("${sample_id}_final_package.zip"), emit: zip

    script:
    """
    PKG="${sample_id}_final_package"
    mkdir -p "\$PKG"

    cp ${assembly_summary} "\$PKG/${sample_id}_assembly_summary.md"
    cp ${tools_report}     "\$PKG/pipeline_tools.md"
    cp final_scaffold.fasta "\$PKG/${sample_id}_nuclear_scaffold.fasta"
    cp final_scaffold.agp   "\$PKG/${sample_id}_nuclear_scaffold.agp"
    cp chloroplast.fasta    "\$PKG/${sample_id}_chloroplast.fasta"
    cp mitochondria.fasta   "\$PKG/${sample_id}_mitochondria.fasta"
    cp ${quast_nuclear_dir}/report.html "\$PKG/quast_nuclear_report.html"
    cp ${quast_cp_dir}/report.html      "\$PKG/quast_chloroplast_report.html"
    cp ${quast_mt_dir}/report.html      "\$PKG/quast_mitochondria_report.html"
    cp busco_final_summary.txt          "\$PKG/busco_final_summary.txt"

    if [ "${has_blobplots}" = "true" ]; then
        mkdir -p "\$PKG/blobtools"
        cp ${blob_pngs} "\$PKG/blobtools/" 2>/dev/null || true
    fi

    if [ "${has_decontam}" = "true" ]; then
        cp decontam.fasta            "\$PKG/${sample_id}_nuclear_decontam.fasta"
        cp contamination_audit.tsv   "\$PKG/${sample_id}_contamination_audit.tsv"
    fi

    if [ "${has_sylph}" = "true" ]; then
        cp sylph_corroboration.tsv   "\$PKG/${sample_id}_sylph_corroboration.tsv"
    fi

    cat > "\$PKG/README.txt" << 'EOF'
${sample_id} -- final results package
======================================

${sample_id}_assembly_summary.md   -- narrative QC summary (start here)
pipeline_tools.md                  -- tool/version rationale for the whole pipeline
${sample_id}_nuclear_scaffold.fasta -- final chromosome-scaffolded nuclear genome
${sample_id}_nuclear_scaffold.agp   -- maps every sequence to the reference chromosome
                                       RagTag placed it on (column 1 = "ChrNN_RagTag"),
                                       or shows it under its own contig name if unplaced.
${sample_id}_chloroplast.fasta, ${sample_id}_mitochondria.fasta -- organelle genomes
quast_*_report.html                 -- open in a browser; self-contained, no other files needed
busco_final_summary.txt             -- gene-space completeness for the published final genome
blobtools/ (if present)             -- contamination-screen blob plot(s); "_unplaced" in the
                                       filename means only contigs RagTag could NOT place on a
                                       reference chromosome were screened (see assembly summary).
EOF

    python3 -m zipfile -c "\${PKG}.zip" "\$PKG"
    """
}

process TOOLS_REPORT {
    label         'qc'
    errorStrategy 'ignore'
    publishDir    "${params.outdir}/reports", mode: 'copy'
    container     'quay.io/biocontainers/multiqc:1.25.1--pyhdfd78af_0'

    output:
    path "pipeline_tools.md"

    script:
    """
    cat > pipeline_tools.md << 'TOOLSEOF'
# Pipeline tools -- Canopy

Describes the role and rationale of each tool in the Canopy ONT plant genome assembly pipeline.

---

## Read QC and filtering

### NanoPlot
Generates statistics and plots for the raw ONT read set: read length histogram, N50, total
base count, and per-read quality scores (Phred scale). Run before filtering so the raw
distribution is visible in MultiQC alongside filtered results.

### Filtlong
Filters reads by minimum length (>=1 kb) and minimum mean quality. Removing short, low-quality
reads reduces assembly fragmentation and lowers the chance of incorporating sequencing errors
into the draft assembly. The filtered reads are used for all downstream steps.

---

## Organelle read separation

### minimap2 (organelle alignment)
Aligns filtered reads to the combined chloroplast + mitochondria reference using the map-ont
preset (tuned for noisy long reads against a short reference). Reads that align are extracted
as cp/mt read sets; unmapped reads form the nuclear read set.

### samtools
Sorts and indexes the organelle alignment BAM so reads can be extracted by genomic region.

---

## Organelle assembly

### Flye (organelle mode)
De novo assembler for long reads. Assembled separately on deduplicated cp and mt read sets.
Produces FASTA + assembly graph (GFA format). Used when --organelle_assembler flye (default).

### OATK (--organelle_assembler oatk)
HMM-based organelle assembler. Identifies organelle reads directly from the full filtered
read set using pHMM profiles from the OatkDB (embryophyta plant gene database), then assembles
them. More sensitive than minimap2-based extraction for highly divergent species. Produces GFA
files for both chloroplast and mitochondria simultaneously.

### Bandage
Renders the GFA assembly graph as a PNG image. A complete circular chloroplast typically
appears as a single circular node or a two-bubble structure (reflecting the inverted repeat).
Useful for a quick sanity check on organelle topology.

### FILTER_ORGANELLE_CONTIGS
After assembly, aligns organelle contigs back to the reference and retains only those with
>=50% query coverage and >=70% identity. Removes nuclear-derived contigs (NUPTs/NUMTs) that
Flye or OATK occasionally incorporates into the organelle assembly. Controlled by
--filter_organelles (default: true), --organelle_min_qcov, and --organelle_min_ident.

### Medaka (organelle)
Same polisher as the nuclear stage (see below), run separately on the filtered chloroplast and
mitochondrial contigs using their own deduplicated read sets. Corrects raw ONT base-level errors
before QUAST scoring — without this step, organelle mismatch/indel rates directly reflect
uncorrected basecalling error rather than assembly quality.

---

## Nuclear assembly

### Flye (nuclear mode)
Same assembler as for organelles, but run on the nuclear read set with the full estimated
genome size (--genome_size). Uses a repeat graph approach to handle the high repeat content
of plant genomes.

### Medaka
Neural-network-based polisher trained on ONT signal data. Re-aligns the original reads to
the Flye draft and calls a consensus sequence to correct base-level errors. Improves raw
accuracy from approximately Q20 (Flye output) toward Q30+. The model must match the flowcell
chemistry and basecaller used; set via --medaka_model.

---

## Haplotig removal (always runs; --final_assembly selects the final genome)

### purge_dups
Removes redundant haplotig contigs from diploid assemblies. In a heterozygous organism, Flye
may assemble both haplotypes of a locus as separate contigs. purge_dups identifies redundant
contigs by two criteria:

1. Coverage depth: reads map to a haplotig at roughly half the expected haploid depth because
   reads from both haplotypes align to the primary contig, leaving the haplotig undercovered.
   Coverage cutoffs define depth zones: junk / haploid-kept / haplotig-purged / repeat.

2. Self-alignment overlap: contigs that align redundantly against each other at >50% overlap
   are flagged as haplotigs regardless of depth class. This step drives most of the removal.

The pipeline uses autotune mode by default: the haploid coverage peak is detected from the
depth distribution (PB.stat) and thresholds are set automatically. Use --calcuts_args to
override the cutoffs.

purge_dups ALWAYS runs, but it is a data-dependent judgment call (autotune can over-purge at
low/uneven coverage, discarding unique sequence rather than haplotigs). The pipeline therefore
scaffolds and BUSCOs BOTH the Medaka (unpurged) and purged genomes for side-by-side comparison
in the QUAST table and the assembly summary, and --final_assembly selects which one becomes the
published final assembly (gets the contaminant screen + FINAL_SUMMARY):

- --final_assembly medaka (DEFAULT): the unpurged Medaka genome is the final.
- --final_assembly purge: the purge_dups genome is the final.

(The legacy --skip_purge flag is retired; pipeline errors with the --final_assembly equivalent.)

### HapDup (--run_hapdup)
Phases the purged primary assembly into haplotype-resolved contigs using long-read signal.
Aligns reads to the purged assembly then separates them by haplotype using heterozygous variant
sites. Produces two haplotype FASTA files. More compute-intensive than purge_dups.

---

## Scaffolding (requires --nuclear_ref)

### RagTag
Orders and orients purged contigs into chromosome-scale scaffolds by alignment to a reference
genome. Two steps are run in sequence:

- ragtag.py correct (enabled by default via --ragtag_correct): breaks contigs at positions
  where the assembly disagrees with the reference, reducing misassemblies before joining.
- ragtag.py scaffold: orders corrected contigs into scaffolds, inserting N-gaps at joins.

The reference does not need to be from the same species -- a related genome provides useful
synteny information even at moderate sequence identity.

---

## Assembly QC

### QUAST
Computes standard assembly statistics (N50, total length, contig count, GC content) and,
when a reference is provided, alignment-based metrics: genome fraction covered,
misassemblies, and structural variants. Run simultaneously on all nuclear assembly stages
(Flye, Medaka, purge_dups/medaka_nopurge, scaffold) to track how each step affects quality.

### BUSCO
Searches the final scaffolded assembly for a lineage-specific set of near-universal
single-copy orthologs. Completeness categories:
- Complete single-copy (S): present exactly once -- ideal
- Complete duplicated (D): present in multiple copies -- may indicate residual haplotigs
- Fragmented (F): partially recovered -- may indicate assembly fragmentation
- Missing (M): absent -- may indicate genuinely missing sequence or lineage mismatch

Use --busco_lineage to select the appropriate database for your organism.

### MultiQC
Aggregates QC reports from NanoPlot, BUSCO, QUAST, and (optionally) Qualimap into a single
interactive HTML report. Enables comparison across samples when multiple samples are processed
in the same batch.

---

## Optional BAM-level QC (final genome only)

These run only on the published final assembly (selected by --final_assembly), using a single
read-to-assembly BAM.

### Qualimap bamqc (--run_qualimap)
Generates per-base coverage statistics, GC bias plots, and insert-size distributions from the
read-to-assembly BAM. Results are integrated into the MultiQC report automatically via native
Qualimap support.

### BlobTools (--run_blobtools / --run_kraken2)
Generates GC-vs-coverage blob plots from the read-to-assembly BAM and assembly FASTA. In
coverage-only mode (--run_blobtools, no taxonomy DB) it visualises whether all contigs cluster
at the expected depth and GC content. With --run_kraken2 it additionally renders a
taxonomy-coloured blob plot from the Kraken2 per-contig hits + the NCBI taxdump, so
non-Viridiplantae contigs (candidate contamination) stand out by phylum. Results are published
as PNG plots and TSV tables in the output directory.

### Kraken2 (--run_kraken2)
Classifies each assembled contig of the final genome against the Kraken2 reference database
(see Reference databases), assigning an NCBI taxon per contig. The per-contig taxids are written
as a BlobTools "hits" file that colours the blob plot above. Judge real contaminants on the blob
plot, where non-Viridiplantae phyla separate from the host cloud by GC and coverage.

---

## Reference databases

### Kraken2 PlusPFP (--kraken2_db)
Prebuilt Kraken2 database covering the Standard content (bacteria, archaea, viral, human,
plasmid, UniVec) plus RefSeq protozoa, fungi, and plant (~231.5 GB loaded index, 171.8 GB
compressed as of the 2026-06-26 build). Unlike PlusPF, PlusPFP includes a plant clade, so host
(sorghum) sequence gets a real match instead of being recognised only by absence. Loaded fully
into RAM (no --memory-mapping) under the assemble_heavy label (448 GB, restricted to specific
big-RAM nodes) rather than qc_heavy (128 GB) -- see nextflow.config. PlusPFP ships its own
nodes.dmp/names.dmp, reused as the BlobTools taxdump (--taxdump_dir) so classifier and taxonomy
share one vintage. Source: genome-idx (Langmead prebuilt indexes). The larger NCBI "reference" DB
(~456 GB) is deliberately not used: it exceeds even assemble_heavy's RAM, forcing
--memory-mapping whose per-k-mer random reads over the shared filesystem never finished.

### OatkDB (--organelle_assembler oatk)
HMM profile database of organelle genes for OATK. A sorghum-specific DB is built from the
provided cp/mt reference FASTAs when available; otherwise the generic embryophyta (land-plant)
OatkDB is downloaded.

### Reference genomes (--cp_ref / --mt_ref / --nuclear_ref)
Sorghum bicolor references: chloroplast NC_008602, mitochondrion NC_008360, and nuclear
Sbicolor_730_v5.0. Used for organelle read/contig separation (minimap2 + FILTER_ORGANELLE_CONTIGS),
RagTag scaffolding, and QUAST genome-fraction. A related genome suffices -- exact species match
is not required.

### BUSCO lineage (--busco_lineage)
Lineage-specific near-universal single-copy ortholog set (default poales_odb10) used to score
the gene-space completeness of each candidate assembly.

---

## References and links

Citations for each tool, in order of appearance, with project homepages.

- **NanoPlot / NanoPack** -- De Coster W, et al. NanoPack: visualizing and processing long-read sequencing data. Bioinformatics. 2018;34(15):2666-2669. <https://github.com/wdecoster/NanoPlot>
- **Filtlong** -- Wick RR. Filtlong (no associated publication). <https://github.com/rrwick/Filtlong>
- **minimap2** -- Li H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics. 2018;34(18):3094-3100. <https://github.com/lh3/minimap2>
- **SAMtools** -- Danecek P, et al. Twelve years of SAMtools and BCFtools. GigaScience. 2021;10(2):giab008. <https://www.htslib.org>
- **Flye** -- Kolmogorov M, et al. Assembly of long, error-prone reads using repeat graphs. Nat Biotechnol. 2019;37:540-546. <https://github.com/fenderglass/Flye>
- **OATK** -- Zhou C. Oatk: organelle assembly toolkit (no associated publication). <https://github.com/c-zhou/oatk>
- **Bandage** -- Wick RR, et al. Bandage: interactive visualization of de novo genome assemblies. Bioinformatics. 2015;31(20):3350-3352. <https://rrwick.github.io/Bandage>
- **Medaka** -- Oxford Nanopore Technologies. Medaka (no associated publication). <https://github.com/nanoporetech/medaka>
- **purge_dups** -- Guan D, et al. Identifying and removing haplotypic duplication in primary genome assemblies. Bioinformatics. 2020;36(9):2896-2898. <https://github.com/dfguan/purge_dups>
- **HapDup** -- Kolmogorov M, et al. Scalable nanopore sequencing of human genomes provides a comprehensive view of haplotype-resolved variation and methylation. Nat Methods. 2023;20:1483-1492. <https://github.com/KolmogorovLab/hapdup>
- **RagTag** -- Alonge M, et al. Automated assembly scaffolding using RagTag elevates a new tomato system for high-throughput genome editing. Genome Biol. 2022;23:258. <https://github.com/malonge/RagTag>
- **QUAST** -- Gurevich A, et al. QUAST: quality assessment tool for genome assemblies. Bioinformatics. 2013;29(8):1072-1075. <https://github.com/ablab/quast>
- **BUSCO** -- Manni M, et al. BUSCO update: novel and streamlined workflows along with broader and deeper phylogenetic coverage. Mol Biol Evol. 2021;38(10):4647-4654. <https://busco.ezlab.org>
- **MultiQC** -- Ewels P, et al. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016;32(19):3047-3048. <https://multiqc.info>
- **Qualimap** -- Okonechnikov K, et al. Qualimap 2: advanced multi-sample quality control for high-throughput sequencing data. Bioinformatics. 2016;32(2):292-294. <http://qualimap.conf.es>
- **BlobTools** -- Laetsch DR, Blaxter ML. BlobTools: Interrogation of genome assemblies. F1000Research. 2017;6:1287. <https://github.com/DRL/blobtools>
- **Kraken2** -- Wood DE, Lu J, Langmead B. Improved metagenomic analysis with Kraken 2. Genome Biol. 2019;20:257. <https://github.com/DerrickWood/kraken2> (PlusPF prebuilt DB: <https://benlangmead.github.io/aws-indexes/k2>)
TOOLSEOF
    """
}
