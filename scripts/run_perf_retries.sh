#!/usr/bin/env bash
# Re-run timed-out perf executables with smaller sizes/iters and a longer
# per-routine cap. Appends results to the same TSV/JSON as the main sweep.
#
# Reads the [fail] lines from reports/perf_sweep.log and looks up the
# matching executable in /tmp/stage-{e,q,m}/build/tests/blas_parallel/.
set -u

OUTDIR="${OUTDIR:-reports}"
JSON="$OUTDIR/perf_sweep.json"
TSV="$OUTDIR/perf_sweep.tsv"
LOG="$OUTDIR/perf_sweep.log"
TIMEOUT=${TIMEOUT:-900}
SIZES=${SIZES:-64,128,256}
ITERS=${ITERS:-3}

grep "^\[fail\]" "$LOG" | sed -E 's|\[fail\] ([a-z0-9]+)/perf_([a-z0-9]+).*|\1 \2|' | sort -u | while read -r tname rname; do
    case "$tname" in
        kind10) tdir=/tmp/stage-e/build/tests/blas_parallel ;;
        kind16) tdir=/tmp/stage-q/build/tests/blas_parallel ;;
        multifloats) tdir=/tmp/stage-m/build/tests/blas_parallel ;;
    esac
    exe="$tdir/perf_$rname"
    [[ -x "$exe" ]] || { echo "[skip] $tname/$rname (no exe)" >&2; continue; }
    echo "[retry] $tname/perf_$rname (SIZES=$SIZES, ITERS=$ITERS)" >&2
    TMP=$(mktemp)
    if BLAS_PERF_JSON="$JSON" BLAS_PERF_SIZES="$SIZES" BLAS_PERF_ITERS="$ITERS" \
           OMP_NUM_THREADS=1 timeout "$TIMEOUT" taskset -c 0 "$exe" > "$TMP" 2>>"$LOG"; then
        : ok
    else
        echo "[retry-fail] $tname/$rname exit=$?" >> "$LOG"
    fi
    awk -v t="$tname" -v r="perf_$rname" '
        /^#/ {next}
        NF >= 6 {
            gsub(/x$/, "", $7);
            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", t, $1, $2, $3, $4, $5, $6, $7;
        }' "$TMP" >> "$TSV"
    rm -f "$TMP"
done

echo "retries appended to $TSV"
