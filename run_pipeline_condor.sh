#!/bin/bash
# Condor-job executable for the Canopy Nextflow head process.
# Submitted via pipeline.condor (universe = local) so it runs ON the submit
# node (scarcity-ap-1) where Nextflow can reach the HTCondor schedd.
#
# Unlike run_pipeline.sh, this does NOT use nohup/backgrounding: Condor owns
# the process lifecycle, so Nextflow must run in the FOREGROUND. If it were
# backgrounded, Condor would see the wrapper exit and tear the job down,
# killing the Nextflow head. Stdout/stderr are captured by Condor (see the
# output/error lines in pipeline.condor).
set -euo pipefail

export HOME=/mnt/cephfs/linuxhome/benucci
CONDA_ENV="${HOME}/.conda/envs/nextflow"
PROJECT_DIR=/mnt/cephfs/linuxhome/benucci/Canopy
export JAVA_HOME="${CONDA_ENV}"
export PATH="${PROJECT_DIR}/bin:${CONDA_ENV}/bin:${PATH}"

cd "${PROJECT_DIR}"

# Each sample has its own work and output directory, so it can be resumed
# independently if the Condor head job is interrupted.
SAMPLES=(B11077 BTx2932 F07020 F10702 F25101)

for sample in "${SAMPLES[@]}"; do
    echo "=== Starting ${sample} ==="
    "${CONDA_ENV}/bin/nextflow" run "${PROJECT_DIR}/main.nf" \
        -profile condor \
        --reads "reads/${sample}/" \
        --cp_ref refs/sorghum/sorghum_cp_NC008602.fasta \
        --mt_ref refs/sorghum/sorghum_mt_NC008360.fasta \
        --nuclear_ref refs/sorghum/Sbicolor_730_v5.0.fa \
        --outdir "results_${sample}" \
        -w "nf-work-${sample}" \
        --organelle_assembler oatk \
        --run_qualimap true \
        --run_blobtools true \
        --run_kraken2 true \
        --flag_contaminants true \
        --verify_sylph true \
        --verify_blast true \
        --blast_remote true \
        --final_assembly "${FINAL_ASSEMBLY:-medaka}" \
        -resume
    echo "=== Finished ${sample} ==="
done
