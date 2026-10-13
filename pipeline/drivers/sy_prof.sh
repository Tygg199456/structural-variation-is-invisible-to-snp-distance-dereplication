#!/usr/bin/env bash
set -uo pipefail
TASK="${1:-sy_prof}"; TASKDIR="$HOME/work/SY"; R="$TASKDIR/run/$TASK"
mkdir -p "$TASKDIR/logs" "$R/steps"; date +%s > "$R/start"
trap 'echo "$? $(date +%s)" > "$R/exit"' EXIT
export PYTHONUNBUFFERED=1
cd "$TASKDIR/pipeline/plague-decay" || exit 2
echo "=== deletion profile dump start $(date) ==="
timeout -k 60 7200 <ENV_PREFIX>/dsm13/bin/python /tmp/dump_profiles.py
echo "rc=$?"
echo "=== done $(date) ==="
