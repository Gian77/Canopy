#!/usr/bin/env bash
# Watch .nextflow.log in real-time and surface task completions, errors, and warnings.
# Usage: ./watch_pipeline.sh [path/to/.nextflow.log]

LOG="${1:-.nextflow.log}"

if [[ ! -f "$LOG" ]]; then
    echo "Log not found: $LOG"
    echo "Start the pipeline first, then re-run this script."
    exit 1
fi

RED='\033[0;31m'; YELLOW='\033[0;33m'; GREEN='\033[0;32m'; CYAN='\033[0;36m'; RESET='\033[0m'

echo -e "${CYAN}=== Watching $LOG (Ctrl-C to stop) ===${RESET}"

tail -F "$LOG" 2>/dev/null | while IFS= read -r line; do
    # Task completed — show exit status coloured
    if [[ "$line" =~ "Task completed" ]]; then
        if [[ "$line" =~ "exit: 0" ]]; then
            name=$(echo "$line" | grep -oP "name: \K[^;]+")
            echo -e "${GREEN}[OK]${RESET}  $name"
        else
            exit_code=$(echo "$line" | grep -oP "exit: \K[^;]+")
            name=$(echo "$line" | grep -oP "name: \K[^;]+")
            echo -e "${RED}[FAIL exit=$exit_code]${RESET}  $name"
        fi

    # Error lines
    elif [[ "$line" =~ "ERROR" ]]; then
        echo -e "${RED}[ERROR]${RESET} $(echo "$line" | sed 's/.*ERROR //')"

    # Submitted — show what just launched
    elif [[ "$line" =~ "Submitted process" ]]; then
        proc=$(echo "$line" | grep -oP "Submitted process > \K.*")
        echo -e "${CYAN}[SUB]${RESET}  $proc"

    # Retry notice
    elif [[ "$line" =~ "Execution will be retried" ]]; then
        echo -e "${YELLOW}[RETRY]${RESET} $(echo "$line" | sed 's/.*INFO //')"

    # Pipeline done
    elif [[ "$line" =~ "Pipeline completed" ]]; then
        if [[ "$line" =~ "SUCCESS" ]]; then
            echo -e "${GREEN}[DONE]${RESET} $(echo "$line" | grep -oP 'Pipeline completed.*')"
        else
            echo -e "${RED}[DONE]${RESET} $(echo "$line" | grep -oP 'Pipeline completed.*')"
        fi

    # Queue warnings
    elif [[ "$line" =~ "WARN" ]] && ! [[ "$line" =~ "DEBUG" ]]; then
        echo -e "${YELLOW}[WARN]${RESET} $(echo "$line" | sed 's/.*WARN[[:space:]]*//')"
    fi
done
