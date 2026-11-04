#!/bin/bash
# Simulate, run BEAST inference and make validation graphs for both setups.
# Run from the skyline-SRFBD-validation folder:
#   bash pipeline.sh
#
# This script finishes once the jobs are submitted; the SLURM jobs keep running
# afterwards. Check progress with:  squeue -u $USER
set -e

mkdir -p runs/slurm_logs

# The run simulate.R just wrote to is always the highest-numbered run_NUM
# (overwrite = FALSE creates highest + 1; overwrite = TRUE reuses the highest).
latest_run() {
    ls -1 runs | grep -E '^run_[0-9]+$' | sed 's/^run_//' | sort -n | tail -n 1
}

# For one run, submit two chained jobs:
#   1. BEAST inference (array job)
#   2. make graphs, once every array task has finished (afterany: also runs if some failed;
#      the graphing script reads logs directly from runs/run_NUM/inf/ and skips missing ones)
submit_run() {
    local run=$1
    local logs="runs/run_${run}/slurm_logs"
    mkdir -p "$logs"

    local inf_job plot_job
    inf_job=$(sbatch --parsable run_BEAST_inferences.sh "$run")
    inf_job=${inf_job%%;*}   # strip cluster name if present

    plot_job=$(sbatch --parsable \
           --dependency=afterany:"$inf_job" \
           --job-name="plots_${run}" \
           --time=04:00:00 --mem=8000M \
           --output="$logs/plots-%j.out" \
           --error="$logs/plots-%j.err" \
           --wrap="Rscript graphing/validation_plots.R ${run}")
    plot_job=${plot_job%%;*}

    echo "run_${run}: inference ${inf_job} -> graphs ${plot_job}"
}

# four intervals (the defaults)
Rscript simulate.R
RUN_FOUR=$(latest_run)
echo "Four-interval simulation: run_${RUN_FOUR}"
submit_run "$RUN_FOUR"

# one interval
Rscript simulate.R --times=0 --inference_template=templates/ssRanges_inference_template_one_int.xml
RUN_ONE=$(latest_run)
echo "One-interval simulation: run_${RUN_ONE}"
submit_run "$RUN_ONE"

# record which run is which
echo "$(date '+%Y-%m-%d %H:%M')  four_intervals=run_${RUN_FOUR}  one_interval=run_${RUN_ONE}" >> runs/pipeline_runs.txt
echo "Run numbers saved to runs/pipeline_runs.txt"
echo "Graphs will appear in runs/run_${RUN_FOUR}/figures/ and runs/run_${RUN_ONE}/figures/"
echo "All jobs submitted - they keep running after this script exits. Check with: squeue -u \$USER"