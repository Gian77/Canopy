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

# Resume from the most recent session by default. Override with RESUME_SESSION to
# pin a specific session UUID: a plain `nextflow run -preview` (or any extra run)
# against this work dir starts a NEW session, and bare `-resume` then targets that
# newer (cache-empty) session, causing a needless full re-run. Pin the good
# session's UUID (from `.nextflow/history`) to reuse its cache.
RESUME="-resume"
[ -n "${RESUME_SESSION:-}" ] && RESUME="-resume ${RESUME_SESSION}"

exec "${CONDA_ENV}/bin/nextflow" run "${PROJECT_DIR}/main.nf" \
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
    --flag_contaminants \
    --final_assembly ${FINAL_ASSEMBLY:-medaka} \
    ${RESUME}
