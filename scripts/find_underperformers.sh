#!/usr/bin/env bash
# Catalog sub-parity overlay routines from a perf_sweep.tsv.
# Usage: ./scripts/find_underperformers.sh [tsv-path] [ratio-threshold]
set -eu
TSV="${1:-reports/full-omp1-fresh/perf_sweep.tsv}"
THR="${2:-0.95}"

# Group by (target, routine), show worst ratio per group + cell count.
echo "=== Worst cell per routine, ratio < $THR (iters >= 100 only) ==="
awk -F'\t' -v thr="$THR" '
NR>1 && int($5) >= 100 && $8+0 < thr {
    key = $1 "\t" $2
    if (!(key in min_ratio) || $8+0 < min_ratio[key]+0) {
        min_ratio[key] = $8
        worst_cell[key] = sprintf("%-12s N=%-5s", $3, $4)
    }
    cnt[key]++
}
END {
    n = 0
    for (k in min_ratio) {
        items[n++] = sprintf("%.3f\t%s\t%s\t(slow cells: %d)\n", min_ratio[k], k, worst_cell[k], cnt[k])
    }
    # sort by ratio (column 1)
    asort(items)
    for (i = 1; i <= n; i++) printf "%s", items[i]
}' "$TSV"

echo
echo "=== Cells <$THR with iters>=100, sorted by ratio ==="
awk -F'\t' -v thr="$THR" '
NR>1 && int($5) >= 100 && $8+0 < thr {
    printf "%.3f\t%-12s %-10s %-14s N=%-5s ov=%s mig=%s\n", $8, $1, $2, $3, $4, $6, $7
}' "$TSV" | sort -n | head -30
