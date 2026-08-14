#!/bin/bash
# Helpers for finding and verifying quay.io biocontainer images.
#
# Functions:
#   quay_tags <tool> [limit]            — list available tags via API (fast)
#   image_exists <image>                — check tag exists via API (fast, no download)
#   verify_image <image> [cmd]          — actually pull image; optionally run cmd inside
#
# Examples:
# Check pullability AND that medaka_consensus works:
#   check_image quay.io/biocontainers/medaka:2.2.1--py310h237e959_0 "medaka --version"
# For Flye:
#   check_image quay.io/biocontainers/flye:2.9.6--py313h7fbb527_1 "flye --version"
#
# Tip: set VERIFY_VERBOSE=1 to see errors from verify_image.

quay_tags() {
    local tool="$1"
    local limit="${2:-20}"
    curl -fsSL \
        "https://quay.io/api/v1/repository/biocontainers/${tool}/tag/?limit=${limit}&onlyActiveTags=true" \
    | python3 -c "
import sys, json
tool = '${tool}'
data = json.load(sys.stdin)
for tag in data.get('tags', []):
    print(f'quay.io/biocontainers/{tool}:{tag[\"name\"]}')
"
}

# Fast API-only existence check (no download)
image_exists() {
    local image="$1"
    local repo="${image#quay.io/}"
    local name="${repo%:*}"
    local tag="${repo##*:}"
    if curl -fsSL \
        "https://quay.io/api/v1/repository/${name}/tag/?specificTag=${tag}&onlyActiveTags=true" \
        | python3 -c "import sys, json; sys.exit(0 if json.load(sys.stdin).get('tags') else 1)"
    then
        echo "[OK] $image"
        return 0
    else
        echo "[FAIL] $image"
        return 1
    fi
}

# Full pull + optional binary check (slow; downloads the image)
# Usage: VERIFY_VERBOSE=1 verify_image quay.io/biocontainers/medaka:2.2.1--py310h237e959_0 "medaka_consensus --help"
verify_image() {
    local img=$1
    local cmd=${2:-}
    local verbose=${VERIFY_VERBOSE:-0}
    local tmpdir="${APPTAINER_TMPDIR:-${HOME}/.apptainer_cache}"
    local sif
    sif=$(mktemp "${tmpdir}/check.XXXXXX.sif")

    if [[ "$verbose" == "1" ]]; then
        apptainer pull --force "$sif" "docker://$img"
    else
        apptainer pull --force "$sif" "docker://$img" >/dev/null 2>&1
    fi
    local rc=$?
    if [[ $rc -ne 0 ]]; then
        echo "[FAIL pull] $img  (re-run with VERIFY_VERBOSE=1 to see errors)"
        rm -f "$sif"
        return 1
    fi

    if [[ -n "$cmd" ]]; then
        if apptainer exec "$sif" bash -c "$cmd" >/dev/null 2>&1; then
            echo "[OK] $img  ($cmd works)"
        else
            echo "[FAIL exec] $img  ($cmd failed; re-run with VERIFY_VERBOSE=1)"
            rm -f "$sif"
            return 2
        fi
    else
        echo "[OK pull] $img  (no command check)"
    fi
    rm -f "$sif"
    return 0
}
