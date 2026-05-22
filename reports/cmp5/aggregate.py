#!/usr/bin/env python3
"""Aggregate the 4-variant perf sweep into a comparison table.

Input columns (cmp5_raw.tsv):
    run_id run_binary omp taskset routine key size iters overlay_GFs migrated_GFs

The same C source is linked into both binaries, so for each (routine, key,
size) we get four overlay GF/s readings (one per run_id) and four
migrated GF/s readings (the migrated_ symbol timing inside each run).

Output:
    cmp5.tsv  — wide table with columns
        routine key size epopenblas-omp1 epopenblas-omp4
                  parallel-blas-omp1 parallel-blas-omp4 migrated-serial
                  mig_par_omp1 mig_par_omp4 mig_ep_omp1 mig_ep_omp4

`migrated-serial` is the migrated_GFs from the parallel-blas-omp1 run.
The mig_* columns expose all four migrated readings so we can spot-check
that the migrated_ symbol's timing is stable across runs (it should be —
`_serial` contains no OpenMP).
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path


RAW = Path(__file__).parent / "cmp5_raw.tsv"
OUT = Path(__file__).parent / "cmp5.tsv"

OVERLAYS = ["epopenblas-omp1", "epopenblas-omp4",
            "parallel-blas-omp1", "parallel-blas-omp4"]
MIG_KEYS = {
    "epopenblas-omp1":     "mig_ep_omp1",
    "epopenblas-omp4":     "mig_ep_omp4",
    "parallel-blas-omp1":  "mig_par_omp1",
    "parallel-blas-omp4":  "mig_par_omp4",
}


def main():
    overlay = defaultdict(dict)   # (routine, key, size) → {variant: gflops}
    migrated = defaultdict(dict)  # (routine, key, size) → {mig_*: gflops}
    with RAW.open() as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            key = (row["routine"], row["key"], int(row["size"]))
            run = row["run_id"]
            overlay[key][run] = float(row["overlay_GFs"])
            migrated[key][MIG_KEYS[run]] = float(row["migrated_GFs"])

    rows = []
    for k in sorted(overlay):
        routine, ks, size = k
        o = overlay[k]
        m = migrated[k]
        rows.append({
            "routine": routine,
            "key": ks,
            "size": size,
            "epopenblas-omp1":    o.get("epopenblas-omp1", ""),
            "epopenblas-omp4":    o.get("epopenblas-omp4", ""),
            "parallel-blas-omp1": o.get("parallel-blas-omp1", ""),
            "parallel-blas-omp4": o.get("parallel-blas-omp4", ""),
            "migrated-serial":    m.get("mig_par_omp1", ""),
            "mig_par_omp1": m.get("mig_par_omp1", ""),
            "mig_par_omp4": m.get("mig_par_omp4", ""),
            "mig_ep_omp1":  m.get("mig_ep_omp1",  ""),
            "mig_ep_omp4":  m.get("mig_ep_omp4",  ""),
        })

    cols = ["routine", "key", "size",
            "epopenblas-omp1", "epopenblas-omp4",
            "parallel-blas-omp1", "parallel-blas-omp4",
            "migrated-serial",
            "mig_par_omp1", "mig_par_omp4", "mig_ep_omp1", "mig_ep_omp4"]
    with OUT.open("w") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({c: (f"{r[c]:.4f}" if isinstance(r[c], float) else r[c])
                        for c in cols})
    print(f"wrote {OUT}  ({len(rows)} rows)")


if __name__ == "__main__":
    main()
