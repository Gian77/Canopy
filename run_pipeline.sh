#!/bin/bash
# Run the Canopy pipeline on the HTCondor submit node (scarcity-ap-1).
# Nextflow must run here — compute nodes cannot reach the HTCondor schedd.
# Usage: bash run_pipeline.sh
set -euo pipefail

export HOME=/mnt/cephfs/linuxhome/benucci
CONDA_ENV="${HOME}/.conda/envs/nextflow"
export JAVA_HOME="${CONDA_ENV}"
export PATH="${CONDA_ENV}/bin:${PATH}"

PROJECT_DIR=/mnt/cephfs/linuxhome/benucci/Canopy
LOG_DIR="${PROJECT_DIR}/condor_logs"
mkdir -p "${LOG_DIR}"

nohup "${CONDA_ENV}/bin/nextflow" run "${PROJECT_DIR}/main.nf" \
    -profile condor \
    --reads F10702_test/ \
    --cp_ref refs/sorghum/sorghum_cp_NC008602.fasta \
    --mt_ref refs/sorghum/sorghum_mt_NC008360.fasta \
    --nuclear_ref refs/sorghum/Sbicolor_730_v5.0.fa \
    --outdir results_F10702 \
    -w nf-work-F10702 \
    --organelle_assembler oatk \
    --run_qualimap \
    --run_blobtools \
    --run_kraken2 \
    --final_assembly ${FINAL_ASSEMBLY:-medaka} \
    -resume \
    > "${LOG_DIR}/pipeline.stdout.txt" \
    2> "${LOG_DIR}/pipeline.stderr.txt" &

echo $! > "${LOG_DIR}/pipeline.pid"
echo "Pipeline started (PID $(cat ${LOG_DIR}/pipeline.pid))"
echo "stdout : ${LOG_DIR}/pipeline.stdout.txt"
echo "stderr : ${LOG_DIR}/pipeline.stderr.txt"
echo "watch  : tail -f ${LOG_DIR}/pipeline.stdout.txt"
