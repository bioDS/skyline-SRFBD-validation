#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=10000M
#SBATCH --cpus-per-task=1
#SBATCH --output=runs/slurm_logs/slurm-%A_%a.out
#SBATCH --error=runs/slurm_logs/slurm-%A_%a.err
#SBATCH --job-name=sranges-validation
#SBATCH --array=1-200
#SBATCH --signal=B:TERM@120

# Usage (from the skyline-SRFBD-validation folder):
#   mkdir -p runs/slurm_logs                  # once; SLURM won't create it
#   sbatch run_inference.sh                   # highest-numbered run in runs/
#   sbatch run_inference.sh 3                 # runs/run_3
#   sbatch --array=1-50 run_inference.sh 3    # different number of trees
#
# All output from BEAST and this script is written to
# runs/run_NUM/slurm_logs/. The files in runs/slurm_logs/ only catch messages
# SLURM itself writes (time-limit or out-of-memory kills) and are deleted
# when empty.
#
# Each task starts its own private Xvfb. Xvfb picks a free display itself
# (-displayfd), so tasks sharing a node can't collide, and it is always shut
# down and its files removed when the task ends, including on time-limit
# (SLURM sends SIGTERM 120 s before the limit, see --signal above).

BASE_DIR="$SLURM_SUBMIT_DIR"
TASK=$SLURM_ARRAY_TASK_ID
FALLBACK_OUT="$BASE_DIR/runs/slurm_logs/slurm-${SLURM_ARRAY_JOB_ID}_${TASK}.out"
FALLBACK_ERR="$BASE_DIR/runs/slurm_logs/slurm-${SLURM_ARRAY_JOB_ID}_${TASK}.err"

if [ -n "$1" ]; then
    RUN_NAME="run_$1"
else
    RUN_NAME=$(ls -1 "$BASE_DIR/runs" | grep -E '^run_[0-9]+$' | sort -V | tail -n 1)
fi
RUN_DIR="$BASE_DIR/runs/$RUN_NAME"

if [ -z "$RUN_NAME" ] || [ ! -d "$RUN_DIR/inf/${TASK}" ]; then
    # still writing to runs/slurm_logs/ here, so this message is kept
    echo "Missing $RUN_DIR/inf/${TASK}" >&2
    exit 1
fi

# from here on, all output goes into the run folder
mkdir -p "$RUN_DIR/slurm_logs"
exec > "$RUN_DIR/slurm_logs/slurm-${SLURM_ARRAY_JOB_ID}_${TASK}.out" \
    2> "$RUN_DIR/slurm_logs/slurm-${SLURM_ARRAY_JOB_ID}_${TASK}.err"

echo "Run: $RUN_NAME  Task: $TASK  Node: $(hostname)  Started: $(date)"
cd "$RUN_DIR/inf/${TASK}/"

JAR_PATH="$BASE_DIR/ssr.jar"
VERSION_DIR="$BASE_DIR/Intellij_module_versions"
JAVA="$HOME/skyline/zulu17.54.21-ca-fx-jdk17.0.13-linux_x64/bin/java"

# ---------------------------------------------------------------------------
# Clear X leftovers from earlier jobs on this node: lock files whose process
# is gone, and sockets nothing is listening on. Live servers (including other
# tasks running right now) are never touched.
# ---------------------------------------------------------------------------
for f in /tmp/.X*-lock; do
    [ -e "$f" ] || continue
    pid=$(tr -d ' ' < "$f" 2>/dev/null)
    [ -n "$pid" ] && [ ! -d "/proc/$pid" ] && rm -f "$f" 2>/dev/null
done
live=$(grep -o '/tmp/\.X11-unix/X[0-9]*' /proc/net/unix | sed 's/.*X//' | sort -u)
for s in /tmp/.X11-unix/X*; do
    [ -e "$s" ] || continue
    n=${s##*/X}
    grep -qx "$n" <<< "$live" || rm -f "$s" 2>/dev/null
done

# ---------------------------------------------------------------------------
# Start a private Xvfb; it chooses a free display and reports it on fd 3.
# ---------------------------------------------------------------------------
XVFB_LOG="$RUN_DIR/slurm_logs/xvfb-${SLURM_ARRAY_JOB_ID}_${TASK}.log"
fifo=$(mktemp -u /tmp/xvfb-fifo.XXXXXX) && mkfifo "$fifo"
Xvfb -displayfd 3 -screen 0 1280x1024x24 -nolisten tcp 3>"$fifo" 2>"$XVFB_LOG" &
XVFB_PID=$!
read -r DISPLAY_NUM < "$fifo"
rm -f "$fifo"
if [ -z "$DISPLAY_NUM" ]; then
    echo "Xvfb failed to start, see $XVFB_LOG" >&2
    exit 1
fi
export DISPLAY=":$DISPLAY_NUM"
echo "Xvfb running on display $DISPLAY (pid $XVFB_PID)"

stop_xvfb() {
    kill "$XVFB_PID" 2>/dev/null
    wait "$XVFB_PID" 2>/dev/null
    # Only remove files that still belong to our server, in case another
    # task has already claimed this display number.
    lock="/tmp/.X${DISPLAY_NUM}-lock"
    if [ -f "$lock" ] && [ "$(tr -d ' ' < "$lock")" = "$XVFB_PID" ]; then
        rm -f "$lock" "/tmp/.X11-unix/X${DISPLAY_NUM}"
    fi
}
trap stop_xvfb EXIT

# On SIGTERM (time limit warning or scancel), stop BEAST, then carry on to
# the normal end of the script so logging and cleanup still happen.
on_signal() {
    echo "Signal received at $(date), stopping BEAST" >&2
    kill -TERM "$JAVA_PID" 2>/dev/null
}
trap on_signal TERM INT HUP

# ---------------------------------------------------------------------------
# Run BEAST in the background and wait, so signals are handled immediately.
# ---------------------------------------------------------------------------
"$JAVA" \
    -Djava.awt.headless=true -Djava.library.path=/lib/x86_64-linux-gnu/ \
    --add-modules javafx.fxml,javafx.base \
    -jar "$JAR_PATH" -threads -1 -seed 42 \
    -version_file "$VERSION_DIR/version-ssr.xml" \
    -version_file "$VERSION_DIR/version-sr.xml" \
    -version_file "$VERSION_DIR/version-feast.xml" \
    -version_file "$VERSION_DIR/version-beast.xml" \
    -version_file "$VERSION_DIR/version-sa.xml" \
    -version_file "$VERSION_DIR/version-beastLabs.xml" \
    -version_file "$VERSION_DIR/version-mm.xml" \
    -version_file "$VERSION_DIR/version-beastfx.xml" \
    -overwrite ssRanges_inference.xml &
JAVA_PID=$!

wait "$JAVA_PID"
STATUS=$?
# If a signal interrupted the wait, keep waiting until BEAST has really exited
while kill -0 "$JAVA_PID" 2>/dev/null; do
    wait "$JAVA_PID"
    STATUS=$?
done

echo "Finished: $(date)  Exit status: $STATUS"

# remove the SLURM fallback files if nothing was written to them
for f in "$FALLBACK_OUT" "$FALLBACK_ERR"; do
    [ -f "$f" ] && [ ! -s "$f" ] && rm -f "$f"
done

# keep the Xvfb log only if something went wrong
[ "$STATUS" -eq 0 ] && rm -f "$XVFB_LOG"

exit $STATUS