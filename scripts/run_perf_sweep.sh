#!/usr/bin/env bash
# Run every blas_parallel perf_* executable in /tmp/stage-{e,q,m}/build/tests/blas_parallel/
# at OMP=1, pinned to P-core 0, and append the structured output to a TSV.
#
# Usage:
#   scripts/run_perf_sweep.sh
#
# Survives single-routine crashes (heap corruption, assertion failures, etc.)
# so the rest of the sweep continues.
set -u

OUTDIR="${OUTDIR:-bench_reports}"
mkdir -p "$OUTDIR"
JSON="$OUTDIR/perf_sweep.json"
TSV="$OUTDIR/perf_sweep.tsv"
LOG="$OUTDIR/perf_sweep.log"
: > "$JSON"
: > "$TSV"
: > "$LOG"

echo -e "target\troutine\tkey\tsize\titers\toverlay_GFs\tmigrated_GFs\tratio" >> "$TSV"

for target in e q m; do
    case "$target" in
        e) tname=kind10; tdir=/tmp/stage-e/build/tests/blas_parallel ;;
        q) tname=kind16; tdir=/tmp/stage-q/build/tests/blas_parallel ;;
        m) tname=multifloats; tdir=/tmp/stage-m/build/tests/blas_parallel ;;
    esac
    if [[ ! -d "$tdir" ]]; then
        echo "[skip] $tname: no build dir at $tdir" >&2 | tee -a "$LOG"
        continue
    fi
    for exe in "$tdir"/perf_*; do
        [[ -x "$exe" ]] || continue
        name=$(basename "$exe")
        echo "[run] $tname/$name" >&2
        # Run with a hard timeout (5 min) per routine and pipe to awk.
        # Use a temp file to decouple the subprocess from the pipe so a
        # crash in the perf executable doesn't kill the parent shell.
        TMP=$(mktemp)
        if BLAS_PERF_JSON="$JSON" OMP_NUM_THREADS=1 \
              timeout 300 taskset -c 0 "$exe" > "$TMP" 2>>"$LOG"; then
            : ok
        else
            echo "[fail] $tname/$name exit=$?" >> "$LOG"
        fi
        awk -v t="$tname" -v r="$name" '
            /^#/ {next}
            NF >= 6 {
                gsub(/x$/, "", $7);
                printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", t, $1, $2, $3, $4, $5, $6, $7;
            }' "$TMP" >> "$TSV"
        rm -f "$TMP"
    done
done

echo "wrote $TSV, $JSON; log at $LOG"
