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
export PROJECT_DIR CONDA_ENV LOG_DIR

SAMPLES=(B11077 BTx2932 F07020 F10702 F25101)

run_sample() {
    local sample="$1"
    local run_dir="${PROJECT_DIR}/.canopy-run-${sample}"
    mkdir -p "${run_dir}"
    cd "${run_dir}"
    echo "=== Starting ${sample} ==="
    "${CONDA_ENV}/bin/nextflow" run "${PROJECT_DIR}/main.nf" \
        -profile condor \
        --reads "${PROJECT_DIR}/reads/${sample}/" \
        --cp_ref "${PROJECT_DIR}/refs/sorghum/sorghum_cp_NC008602.fasta" \
        --mt_ref "${PROJECT_DIR}/refs/sorghum/sorghum_mt_NC008360.fasta" \
        --nuclear_ref "${PROJECT_DIR}/refs/sorghum/Sbicolor_730_v5.0.fa" \
        --outdir "${PROJECT_DIR}/results_${sample}" \
        -w "${PROJECT_DIR}/nf-work-${sample}" \
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
}

run_all_samples() {
    local pids=()
    for sample in "${SAMPLES[@]}"; do
        run_sample "${sample}" \
            > "${LOG_DIR}/${sample}.stdout.txt" \
            2> "${LOG_DIR}/${sample}.stderr.txt" &
        echo $! > "${LOG_DIR}/${sample}.pid"
        pids+=("$!")
        echo "Started ${sample} (PID $!)"
    done

    local failed=0
    for pid in "${pids[@]}"; do
        wait "${pid}" || failed=1
    done
    return "${failed}"
}

nohup bash -euo pipefail -c "$(declare -p SAMPLES); $(declare -f run_sample run_all_samples); run_all_samples" \
    > "${LOG_DIR}/pipeline.stdout.txt" \
    2> "${LOG_DIR}/pipeline.stderr.txt" &

echo $! > "${LOG_DIR}/pipeline.pid"
echo "Pipeline started (PID $(cat ${LOG_DIR}/pipeline.pid))"
echo "stdout : ${LOG_DIR}/pipeline.stdout.txt"
echo "stderr : ${LOG_DIR}/pipeline.stderr.txt"
echo "watch  : tail -f ${LOG_DIR}/pipeline.stdout.txt"
