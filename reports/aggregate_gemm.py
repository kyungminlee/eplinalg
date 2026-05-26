#!/usr/bin/env python3
"""Comprehensive gemm-only bench table — every (precision, routine,
TRANSA, TRANSB, size) row with parallel-blas / migrated-Fortran GFLOPS and speedup.

Scope: parallel-blas overlay vs migrated Fortran reference at OMP=1.
"""
import re
from pathlib import Path
from collections import defaultdict

ROOT = Path(__file__).parent / "gemm-only"

PREC_NAME = {"k10": "kind10", "k16": "kind16", "m": "multifloats"}
PREC_ORDER = ["k10", "k16", "m"]

rows = []  # (prec, prefix, ta, tb, size, ov, mig, sp)
for f in sorted(ROOT.glob("*-bench_*.txt")):
    m = re.match(r"^(k10|k16|m)-bench_([a-z]+)\.txt$", f.name)
    if not m: continue
    prec, prefix = m.group(1), m.group(2)
    for ln in f.read_text().splitlines():
        if "GFLOP" in ln or "wrote" in ln or not ln.strip(): continue
        # gemm row: "NN       64          0.0556             0.0528       1.05x"
        mt = re.match(r"^\s*([NTC])([NTC])\s+(\d+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.]+)x\s*$", ln)
        if mt:
            rows.append((prec, prefix, mt.group(1), mt.group(2), int(mt.group(3)),
                         float(mt.group(4)), float(mt.group(5)), float(mt.group(6))))

# Group by (prec, prefix) for separate tables.
by_routine = defaultdict(list)
for r in rows:
    by_routine[(r[0], r[1])].append(r)

# Sort each routine's rows: TA, TB, size
out = ["# Comprehensive GEMM bench — parallel-blas vs migrated, OMP=1, all (TA, TB, size) combos",
       "",
       "Sizes vary by precision (kind16 omits s=512 to keep run time reasonable).",
       "Each row: parallel-blas GFLOP/s vs migrated-Fortran GFLOP/s, speedup = `t_migrated / t_parallel-blas` (>1 = parallel-blas wins).",
       "",
       "Scope: parallel-blas overlay only.",
       ""]

for prec in PREC_ORDER:
    for (p, prefix), rs in sorted(by_routine.items()):
        if p != prec: continue
        out.append(f"## {PREC_NAME[prec]} — `{prefix}`")
        out.append("")
        out.append("| TA | TB | size | parallel-blas GF/s | migrated GF/s | speedup |")
        out.append("|----|----|------|-------------------|---------------|---------|")
        for _, _, ta, tb, size, ov, mig, sp in sorted(rs, key=lambda r: (r[2], r[3], r[4])):
            marker = " 🟢" if sp >= 1.0 else " 🔴"
            out.append(f"| {ta} | {tb} | {size} | {ov:.4f} | {mig:.4f} | {sp:.2f}×{marker} |")
        out.append("")

Path(__file__).parent.joinpath("gemm-only-report.md").write_text("\n".join(out))
print(f"Wrote {len(rows)} rows across {len(by_routine)} routines")
