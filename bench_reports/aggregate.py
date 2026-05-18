#!/usr/bin/env python3
"""Aggregate full-omp1 bench outputs into a single markdown report.

Reads bench_reports/full-omp1/{prec}-bench_{routine}.txt files and emits
a per-routine table with all (transpose × size) measurements.
"""
import re
import sys
from pathlib import Path
from collections import defaultdict

ROOT = Path(__file__).parent / "full-omp1"

# Output table per routine: precision, routine, combo, size, GFLOP/s overlay,
# GFLOP/s migrated, speedup.

PREC_NAME = {"k10": "kind10", "k16": "kind16", "m": "multifloats"}
PREC_ORDER = ["k10", "k16", "m"]

ROUTINE_ORDER = ["gemm", "trsm", "trmm", "syrk", "herk", "symm", "hemm"]

def routine_key(name):
    # bench_egemm -> 'gemm'
    for r in ROUTINE_ORDER:
        if r in name:
            return r
    return "?"

def parse_file(p):
    """Yield (combo, size, ov, mig, speedup) tuples from a bench .txt."""
    lines = p.read_text().splitlines()
    # First non-empty line is the header; subsequent lines are rows.
    has_trans = False
    for ln in lines:
        if "GFLOP/s" in ln and "speedup" in ln:
            has_trans = "trans" in ln
            continue
        if ln.startswith("wrote") or not ln.strip():
            continue
        if has_trans:
            # "NN       64          0.0556             0.0528       1.05x"
            m = re.match(r"^\s*([NTC]{2})\s+(\d+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.]+)x\s*$", ln)
            if m:
                yield (m.group(1), int(m.group(2)), float(m.group(3)), float(m.group(4)), float(m.group(5)))
        else:
            # "       64          0.7146             0.4848       1.47x"
            m = re.match(r"^\s*(\d+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.]+)x\s*$", ln)
            if m:
                yield (None, int(m.group(1)), float(m.group(2)), float(m.group(3)), float(m.group(4)))

# Collect all results: {routine: [(prec, prefix, combo, size, ov, mig, sp)]}
results = defaultdict(list)

for f in sorted(ROOT.glob("*-bench_*.txt")):
    fname = f.name
    m = re.match(r"^(k10|k16|m)-bench_([a-z]+)\.txt$", fname)
    if not m:
        continue
    prec, prefix = m.group(1), m.group(2)
    r = routine_key(prefix)
    for combo, size, ov, mig, sp in parse_file(f):
        results[r].append((prec, prefix, combo, size, ov, mig, sp))

# Emit markdown
out = ["# Full OMP=1 BLAS overlay vs migrated benchmark",
       "",
       "Per-routine speedup at OMP=1. Run with iters=3, warmup=1, sizes=64/128/256/512.",
       "",
       "Speedup column: `t_migrated / t_overlay` (>1 = overlay wins).",
       ""]

for r in ROUTINE_ORDER:
    if r not in results:
        continue
    rows = results[r]
    has_combo = any(c is not None for _, _, c, _, _, _, _ in rows)
    out.append(f"## {r}")
    out.append("")
    if has_combo:
        hdr = "| prec | routine | trans | size | overlay GF | migrated GF | speedup |"
        sep = "|------|---------|-------|------|-----------|-------------|---------|"
    else:
        hdr = "| prec | routine | size | overlay GF | migrated GF | speedup |"
        sep = "|------|---------|------|-----------|-------------|---------|"
    out.append(hdr)
    out.append(sep)
    # Group by (prec, prefix) for readability
    rows_sorted = sorted(rows, key=lambda x: (PREC_ORDER.index(x[0]), x[1], x[2] or "", x[3]))
    for prec, prefix, combo, size, ov, mig, sp in rows_sorted:
        if has_combo:
            out.append(f"| {PREC_NAME[prec]} | {prefix} | {combo or '—'} | {size} | {ov:.4f} | {mig:.4f} | {sp:.2f}× |")
        else:
            out.append(f"| {PREC_NAME[prec]} | {prefix} | {size} | {ov:.4f} | {mig:.4f} | {sp:.2f}× |")
    out.append("")

Path(__file__).parent.joinpath("full-omp1-report.md").write_text("\n".join(out))
print(f"Wrote {Path(__file__).parent / 'full-omp1-report.md'}")
print(f"Total rows: {sum(len(v) for v in results.values())}")
