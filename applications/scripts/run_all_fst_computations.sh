#!/bin/bash
# =============================================================================
# FST Computation Pipeline -- Sequential Runner
# Runs all 4 FST computation scripts on ellmos-services after MCMC completes.
#
# Usage:
#   1. Upload to server:  scp run_all_fst_computations.sh root@46.62.243.71:/opt/fst_calculations/
#   2. Run:               ssh root@46.62.243.71 "nohup bash /opt/fst_calculations/run_all_fst_computations.sh > /opt/fst_calculations/pipeline.log 2>&1 &"
#   3. Monitor:           ssh root@46.62.243.71 "tail -f /opt/fst_calculations/pipeline.log"
#
# Estimated total runtime: ~7-8 hours
# =============================================================================

set -euo pipefail
export PYTHONIOENCODING=utf-8

# Paths
BASE="/opt/fst_calculations"
SCRIPTS="$BASE/scripts"
DATA="$BASE/data/pdb"
RESULTS="$BASE/results"
LOG="$BASE/pipeline.log"

# Timestamps
ts() { date '+%Y-%m-%d %H:%M:%S'; }

echo "=============================================="
echo "FST Computation Pipeline"
echo "Started: $(ts)"
echo "=============================================="

TOTAL_START=$SECONDS
FAILED=0
COMPLETED=0

# --- Helper ---
run_step() {
    local STEP_NUM="$1"
    local STEP_NAME="$2"
    local SCRIPT_PATH="$3"
    shift 3
    local ARGS=("$@")

    echo ""
    echo "----------------------------------------------"
    echo "[$STEP_NUM/4] $STEP_NAME"
    echo "Script: $SCRIPT_PATH"
    echo "Start:  $(ts)"
    echo "----------------------------------------------"

    local STEP_START=$SECONDS

    if [ ! -f "$SCRIPT_PATH" ]; then
        echo "ERROR: Script not found: $SCRIPT_PATH"
        FAILED=$((FAILED + 1))
        return 1
    fi

    if python3 "$SCRIPT_PATH" "${ARGS[@]}"; then
        local ELAPSED=$(( SECONDS - STEP_START ))
        local MINS=$(( ELAPSED / 60 ))
        local SECS=$(( ELAPSED % 60 ))
        echo "OK: $STEP_NAME completed in ${MINS}m ${SECS}s"
        COMPLETED=$((COMPLETED + 1))
    else
        local EXIT_CODE=$?
        local ELAPSED=$(( SECONDS - STEP_START ))
        echo "FAILED: $STEP_NAME (exit code $EXIT_CODE) after ${ELAPSED}s"
        FAILED=$((FAILED + 1))
        # Continue with next script instead of aborting
        return 0
    fi
}

# --- Pre-flight checks ---
echo ""
echo "Pre-flight checks..."
echo "  Python:    $(python3 --version 2>&1)"
echo "  NumPy:     $(python3 -c 'import numpy; print(numpy.__version__)' 2>&1)"
echo "  SciPy:     $(python3 -c 'import scipy; print(scipy.__version__)' 2>&1)"
echo "  BioPython: $(python3 -c 'import Bio; print(Bio.__version__)' 2>&1 || echo 'not installed')"
echo "  PDB files: $(ls $DATA/*.pdb 2>/dev/null | wc -l) found"
echo ""

# Create output directories
mkdir -p "$RESULTS/fst_iii/eta_scan"
mkdir -p "$RESULTS/fst_iii/benchmark"
mkdir -p "$RESULTS/fst_iii/extended"
mkdir -p "$RESULTS/fst_i/stellar"

# =============================================================================
# Step 1: Extended Protein Analysis (~30 min)
# =============================================================================
run_step "1" "Extended Nash-Frustration Analysis (3 proteins)" \
    "$SCRIPTS/fst_iii/run_extended_analysis.py" \
    --data-dir "$DATA" \
    --output-dir "$RESULTS/fst_iii/extended"

# =============================================================================
# Step 2: eta-Scan (~3-4 hours)
# =============================================================================
run_step "2" "Systematic eta-Scan (7 values x 3 proteins)" \
    "$SCRIPTS/fst_iii/eta_scan.py" \
    --data-dir "$DATA" \
    --output-dir "$RESULTS/fst_iii/eta_scan"

# =============================================================================
# Step 3: Frustration Benchmark (~1 hour)
# =============================================================================
run_step "3" "ROC/AUPRC Benchmark (Nash vs. functional residues)" \
    "$SCRIPTS/fst_iii/frustration_benchmark.py" \
    --data-dir "$DATA" \
    --output-dir "$RESULTS/fst_iii/benchmark"

# =============================================================================
# Step 4: Stellar Testbed (~2 hours)
# =============================================================================
run_step "4" "Stellar Testbed (alpha variation)" \
    "$SCRIPTS/fst_i/stellar_testbed.py" \
    --output-dir "$RESULTS/fst_i/stellar"

# =============================================================================
# Summary
# =============================================================================
TOTAL_ELAPSED=$(( SECONDS - TOTAL_START ))
TOTAL_HOURS=$(( TOTAL_ELAPSED / 3600 ))
TOTAL_MINS=$(( (TOTAL_ELAPSED % 3600) / 60 ))

echo ""
echo "=============================================="
echo "Pipeline Complete: $(ts)"
echo "=============================================="
echo "  Completed: $COMPLETED / 4"
echo "  Failed:    $FAILED / 4"
echo "  Runtime:   ${TOTAL_HOURS}h ${TOTAL_MINS}m"
echo ""

if [ $FAILED -eq 0 ]; then
    echo "ALL STEPS SUCCESSFUL"
    echo ""
    echo "Results:"
    echo "  $RESULTS/fst_iii/extended/"
    echo "  $RESULTS/fst_iii/eta_scan/"
    echo "  $RESULTS/fst_iii/benchmark/"
    echo "  $RESULTS/fst_i/stellar/"
    echo ""
    echo "Next: Download results to local machine:"
    echo "  scp -r -i ~/.ssh/id_ed25519_mcmc root@46.62.243.71:$RESULTS/ ./results/"
else
    echo "WARNING: $FAILED step(s) failed. Check log above."
fi
