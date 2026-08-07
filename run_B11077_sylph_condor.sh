#!/bin/bash
# Condor-job executable for the B11077 Sylph/decontamination regression run.
# This is intentionally separate from run_pipeline_condor.sh, which is the
# production F10702 launcher.
set -euo pipefail

export HOME=/mnt/cephfs/linuxhome/benucci
CONDA_ENV="${HOME}/.conda/envs/nextflow"
PROJECT_DIR=/mnt/cephfs/linuxhome/benucci/Canopy
export JAVA_HOME="${CONDA_ENV}"
export PATH="${PROJECT_DIR}/bin:${CONDA_ENV}/bin:${PATH}"

cd "${PROJECT_DIR}"

# Bare -resume uses the latest session for this project/work directory. Set
# RESUME_SESSION to pin a specific session after a preview or other stray run.
RESUME="-resume"
[ -n "${RESUME_SESSION:-}" ] && RESUME="-resume ${RESUME_SESSION}"

exec "${CONDA_ENV}/bin/nextflow" run "${PROJECT_DIR}/main.nf" \
    -profile condor \
    --reads tests/B11077_test/ \
    --cp_ref refs/sorghum/sorghum_cp_NC008602.fasta \
    --mt_ref refs/sorghum/sorghum_mt_NC008360.fasta \
    --nuclear_ref refs/sorghum/Sbicolor_730_v5.0.fa \
    --outdir results_B11077 \
    -w nf-work-B11077 \
    --organelle_assembler oatk \
    --run_qualimap true \
    --run_blobtools true \
    --run_kraken2 true \
    --flag_contaminants true \
    --verify_sylph true \
    --final_assembly "${FINAL_ASSEMBLY:-medaka}" \
    ${RESUME}
