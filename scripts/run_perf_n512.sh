#!/usr/bin/env bash
# Fill the N=512 gap for the 18 perf harnesses that timed out at 10 iters.
# Runs each at BLAS_PERF_SIZES=512 BLAS_PERF_ITERS=1 BLAS_PERF_WARMUP=0
# with a generous 900s per-routine cap; appends to the same TSV/JSON.
set -u

OUTDIR="${OUTDIR:-bench_reports}"
JSON="$OUTDIR/perf_sweep.json"
TSV="$OUTDIR/perf_sweep.tsv"
LOG="$OUTDIR/perf_sweep.log"
TIMEOUT=${TIMEOUT:-900}

ROUTINES=(
    "kind10:ytrsm"
    "kind16:qgemm" "kind16:qgemmtr" "kind16:qsymm" "kind16:qtrmm" "kind16:qtrsm"
    "kind16:xgemm" "kind16:xgemmtr" "kind16:xhemm" "kind16:xherk" "kind16:xsymm"
    "kind16:xsyrk" "kind16:xtrmm" "kind16:xtrsm"
    "multifloats:wgemm" "multifloats:whemm" "multifloats:wsymm" "multifloats:wtrmm"
)

for entry in "${ROUTINES[@]}"; do
    tname="${entry%%:*}"
    rname="${entry##*:}"
    case "$tname" in
        kind10) tdir=/tmp/stage-e/build/tests/blas_parallel ;;
        kind16) tdir=/tmp/stage-q/build/tests/blas_parallel ;;
        multifloats) tdir=/tmp/stage-m/build/tests/blas_parallel ;;
    esac
    exe="$tdir/perf_$rname"
    [[ -x "$exe" ]] || { echo "[skip] $tname/$rname (no exe)" >&2; continue; }
    echo "[n512] $tname/perf_$rname" >&2
    TMP=$(mktemp)
    if BLAS_PERF_JSON="$JSON" BLAS_PERF_SIZES=512 BLAS_PERF_ITERS=1 BLAS_PERF_WARMUP=0 \
           OMP_NUM_THREADS=1 timeout "$TIMEOUT" taskset -c 0 "$exe" > "$TMP" 2>>"$LOG"; then
        : ok
    else
        echo "[n512-fail] $tname/$rname exit=$?" >> "$LOG"
    fi
    awk -v t="$tname" '
        /^#/ {next}
        NF >= 6 {
            gsub(/x$/, "", $7);
            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", t, $1, $2, $3, $4, $5, $6, $7;
        }' "$TMP" >> "$TSV"
    rm -f "$TMP"
done

echo "N=512 cells appended to $TSV"
