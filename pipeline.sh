#!/bin/bash
# Simulate, run BEAST inference and make validation graphs.
# Run from the skyline-SRFBD-validation folder:
#   bash pipeline.sh                  # all three setups:
#                                     #   1. four intervals (simulate.R defaults)
#                                     #   2. one interval, ONE_TEMPLATE below
#                                     #   3. one interval, CHECK_TEMPLATE below
#   bash pipeline.sh one              # one interval only, using CHECK_TEMPLATE
#   bash pipeline.sh two              # four interval, and one interval using ONE_TEMPLATE
#   bash pipeline.sh one path/to.xml  # one interval only, with a different template
#   bash pipeline.sh submit 3 4 5     # no simulation: submit inference + graphs
#                                     # for existing runs/run_3, run_4, run_5
#
# This script finishes once the jobs are submitted; the SLURM jobs keep running
# afterwards. Check progress with:  squeue -u $USER

# Created using Claude

set -e
shopt -s inherit_errexit 2>/dev/null || true   # make set -e apply inside $(...) too (bash >= 4.4)

MODE=${1:-all}
ONE_TEMPLATE=templates/ssRanges_inference_template_one_int.xml
CHECK_TEMPLATE=${2:-templates/check_single_int.xml}

if [[ "$MODE" != "all" && "$MODE" != "one" && "$MODE" != "two" && "$MODE" != "submit" ]]; then
    echo "Usage: bash pipeline.sh [all|one|two] [inference_template]" >&2
    echo "       bash pipeline.sh submit RUN [RUN...]" >&2
    exit 1
fi
if [[ "$MODE" != "submit" ]]; then
    case "$MODE" in
        all) templates=("$CHECK_TEMPLATE" "$ONE_TEMPLATE") ;;
        one) templates=("$CHECK_TEMPLATE") ;;
        two) templates=("$ONE_TEMPLATE") ;;
    esac
    for tpl in "${templates[@]}"; do
        if [[ ! -f "$tpl" ]]; then
            echo "Inference template not found: $tpl" >&2
            exit 1
        fi
    done
fi

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
# If either sbatch call fails, the script stops with an error message.
submit_run() {
    local run=$1
    local logs="runs/run_${run}/slurm_logs"
    mkdir -p "$logs"

    local inf_job plot_job
    if ! inf_job=$(sbatch --parsable run_BEAST_inferences.sh "$run"); then
        echo "run_${run}: inference submission FAILED" >&2
        exit 1
    fi
    inf_job=${inf_job%%;*}   # strip cluster name if present

    if ! plot_job=$(sbatch --parsable \
           --dependency=afterany:"$inf_job" \
           --job-name="plots_${run}" \
           --time=04:00:00 --mem=8000M \
           --output="$logs/plots-%j.out" \
           --error="$logs/plots-%j.err" \
           --wrap="Rscript validation_plots.R ${run}"); then
        echo "run_${run}: graphing submission FAILED" >&2
        exit 1
    fi
    plot_job=${plot_job%%;*}

    echo "run_${run}: inference ${inf_job} -> graphs ${plot_job}"
}

# Simulate one setup, submit its jobs, and print the run number.
# Usage: simulate_and_submit LABEL [simulate.R arguments...]
# Progress messages go to stderr so that only the run number is captured.
simulate_and_submit() {
    local label=$1; shift
    Rscript simulate.R "$@" >&2
    local run
    run=$(latest_run)
    echo "${label}: run_${run}" >&2
    submit_run "$run" >&2
    echo "$run"
}

if [[ "$MODE" == "submit" ]]; then
    shift
    if [[ $# -eq 0 ]]; then
        echo "Usage: bash pipeline.sh submit RUN [RUN...]" >&2
        exit 1
    fi
    for run in "$@"; do
        run=${run#run_}   # accept "3" or "run_3"
        if [[ ! -d "runs/run_${run}" ]]; then
            echo "No such run: runs/run_${run}" >&2
            exit 1
        fi
        submit_run "$run"
    done
    echo "$(date '+%Y-%m-%d %H:%M')  resubmitted=$(printf 'run_%s ' "${@#run_}")" >> runs/pipeline_runs.txt
    echo "All jobs submitted - they keep running after this script exits. Check with: squeue -u \$USER"
    exit 0
fi

if [[ "$MODE" == "one" ]]; then
    RUN_CHECK=$(simulate_and_submit "One interval (template: $CHECK_TEMPLATE)" \
                --times=0 --inference_template="$CHECK_TEMPLATE")

    echo "$(date '+%Y-%m-%d %H:%M')  one_interval_check=run_${RUN_CHECK}  template=${CHECK_TEMPLATE}" >> runs/pipeline_runs.txt
    echo "Run number saved to runs/pipeline_runs.txt"
    echo "Graphs will appear in runs/run_${RUN_CHECK}/figures/"
else
    # 1. four intervals (the defaults)
    RUN_FOUR=$(simulate_and_submit "Four intervals")

    # 2. one interval, ssRanges_inference_template_one_int template
    RUN_ONE=$(simulate_and_submit "One interval (template: $ONE_TEMPLATE)" \
              --times=0 --inference_template="$ONE_TEMPLATE")

    if [[ "$MODE" == "all" ]]; then
        # 3. one interval, check_single_int template
        RUN_CHECK=$(simulate_and_submit "One interval (template: $CHECK_TEMPLATE)" \
                    --times=0 --inference_template="$CHECK_TEMPLATE")

        echo "$(date '+%Y-%m-%d %H:%M')  four_intervals=run_${RUN_FOUR}  one_interval=run_${RUN_ONE}  one_interval_check=run_${RUN_CHECK}  one_template=${ONE_TEMPLATE}  check_template=${CHECK_TEMPLATE}" >> runs/pipeline_runs.txt
        echo "Run numbers saved to runs/pipeline_runs.txt"
        echo "Graphs will appear in runs/run_${RUN_FOUR}/figures/, runs/run_${RUN_ONE}/figures/ and runs/run_${RUN_CHECK}/figures/"
    else
        echo "$(date '+%Y-%m-%d %H:%M')  four_intervals=run_${RUN_FOUR}  one_interval=run_${RUN_ONE}  one_template=${ONE_TEMPLATE}" >> runs/pipeline_runs.txt
        echo "Run numbers saved to runs/pipeline_runs.txt"
        echo "Graphs will appear in runs/run_${RUN_FOUR}/figures/ and runs/run_${RUN_ONE}/figures/"
    fi
fi

echo "All jobs submitted - they keep running after this script exits. Check with: squeue -u \$USER"