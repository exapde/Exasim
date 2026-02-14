#!/usr/bin/env bash

set -u
set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_LOG_DIR="${SCRIPT_DIR}/logs"
LOG_DIR="${LOG_DIR:-${DEFAULT_LOG_DIR}}"
TIMESTAMP="$(date +"%Y%m%d-%H%M%S")"

mkdir -p "${LOG_DIR}"

pass_count=0
fail_count=0
failed_cases=()
failed_logs=()

run_step() {
    local logfile="$1"
    shift
    if "$@" 2>&1 | tee -a "${logfile}"; then
        return 0
    fi
    return 1
}

run_case() {
    local case_id="$1"
    local pdeapp_path="$2"
    local solver_mode="$3"
    local np="$4"

    local safe_case="${case_id//\//_}"
    local logfile="${LOG_DIR}/${safe_case}-${TIMESTAMP}.log"

    echo "==================================================" | tee -a "${logfile}"
    echo "Running ${case_id}" | tee -a "${logfile}"
    echo "PDE App: ${pdeapp_path}" | tee -a "${logfile}"
    echo "Log file: ${logfile}" | tee -a "${logfile}"
    echo "==================================================" | tee -a "${logfile}"

    local case_ok=1
    if ! run_step "${logfile}" ../build/text2code "${pdeapp_path}"; then
        case_ok=0
    fi

    if [ "${case_ok}" -eq 1 ]; then
        if [ "${solver_mode}" = "cpu-single" ]; then
            if ! run_step "${logfile}" ../build/cput2cEXASIM "${pdeapp_path}"; then
                case_ok=0
            fi
        else
            if ! run_step "${logfile}" mpirun -np "${np}" ../build/cpumpit2cEXASIM "${pdeapp_path}"; then
                case_ok=0
            fi
        fi
    fi

    if [ "${case_ok}" -eq 1 ]; then
        pass_count=$((pass_count + 1))
        echo "STATUS: PASS (${case_id})" | tee -a "${logfile}"
    else
        fail_count=$((fail_count + 1))
        failed_cases+=("${case_id}")
        failed_logs+=("${logfile}")
        echo "STATUS: FAIL (${case_id})" | tee -a "${logfile}"
    fi

    echo | tee -a "${logfile}"
}

# case_id | pdeapp_path | solver_mode | np
run_case "poisson/poisson1d" "poisson/poisson1d/pdeapp.txt" "cpu-single" "1"
run_case "poisson/poisson2d" "poisson/poisson2d/pdeapp.txt" "mpi-cpu" "4"
run_case "poisson/poisson3d" "poisson/poisson3d/pdeapp.txt" "mpi-cpu" "4"
run_case "poisson/periodic" "poisson/periodic/pdeapp.txt" "mpi-cpu" "2"
run_case "poisson/orion" "poisson/orion/pdeapp.txt" "mpi-cpu" "4"
run_case "poisson/lshape" "poisson/lshape/pdeapp.txt" "mpi-cpu" "4"
run_case "poisson/cone" "poisson/cone/pdeapp.txt" "mpi-cpu" "8"
run_case "navierstokes/isoq" "navierstokes/isoq/pdeapp.txt" "mpi-cpu" "4"
run_case "navierstokes/orion" "navierstokes/orion/pdeapp.txt" "mpi-cpu" "4"
run_case "navierstokes/nsmach8" "navierstokes/nsmach8/pdeapp.txt" "mpi-cpu" "4"
run_case "navierstokes/naca0012steady" "navierstokes/naca0012steady/pdeapp.txt" "mpi-cpu" "4"
run_case "navierstokes/naca0012unsteady" "navierstokes/naca0012unsteady/pdeapp.txt" "mpi-cpu" "8"

echo "==================== Final Summary ===================="
echo "Logs directory: ${LOG_DIR}"
echo "Total cases: $((pass_count + fail_count))"
echo "Passed: ${pass_count}"
echo "Failed: ${fail_count}"

if [ "${fail_count}" -gt 0 ]; then
    echo "Failed cases:"
    for c in "${failed_cases[@]}"; do
        echo "  - ${c}"
    done
    echo "Failed case logs:"
    for l in "${failed_logs[@]}"; do
        echo "  - ${l}"
    done
    exit 1
fi

echo "All app runs completed successfully."
