#!/usr/bin/env python3
"""Compact summary: one row per (precision, routine) showing speedup
range across all (combo × size) measurements."""
import re
from pathlib import Path
from collections import defaultdict
from statistics import median

ROOT = Path(__file__).parent / "full-omp1"
PREC_NAME = {"k10": "kind10", "k16": "kind16", "m": "multifloats"}
PREC_ORDER = ["k10", "k16", "m"]

speedups = defaultdict(list)  # (prec, prefix) -> list of (combo, size, sp)
for f in sorted(ROOT.glob("*-bench_*.txt")):
    m = re.match(r"^(k10|k16|m)-bench_([a-z]+)\.txt$", f.name)
    if not m:
        continue
    prec, prefix = m.group(1), m.group(2)
    for ln in f.read_text().splitlines():
        if "GFLOP/s" in ln or "wrote" in ln or not ln.strip():
            continue
        mt = re.match(r"^\s*([NTC]{2})?\s*(\d+)\s+[\d.eE+-]+\s+[\d.eE+-]+\s+([\d.]+)x\s*$", ln)
        if mt:
            speedups[(prec, prefix)].append((mt.group(1) or "—", int(mt.group(2)), float(mt.group(3))))

out = ["# Compact summary: OMP=1 speedup ranges",
       "",
       "Per (precision, routine), speedup distribution across all measured (combo × size) cases.",
       "",
       "| prec | routine | combos × sizes | min spd | median | max spd | <1.0× cases |",
       "|------|---------|----------------|---------|--------|---------|-------------|"]
for prec in PREC_ORDER:
    for (p, prefix), sps in sorted(speedups.items()):
        if p != prec: continue
        vals = [s for _, _, s in sps]
        n = len(vals)
        bad = sum(1 for v in vals if v < 1.0)
        out.append(f"| {PREC_NAME[prec]} | {prefix} | {n} | {min(vals):.2f}× | {median(vals):.2f}× | {max(vals):.2f}× | {bad}/{n} |")
out.append("")

# Worst cases (regressions)
worst = []
for (prec, prefix), sps in speedups.items():
    for combo, size, sp in sps:
        if sp < 1.0:
            worst.append((sp, PREC_NAME[prec], prefix, combo, size))
worst.sort()
out.append("## Worst regressions (overlay slower than migrated)")
out.append("")
out.append("| speedup | prec | routine | combo | size |")
out.append("|---------|------|---------|-------|------|")
for sp, prec, prefix, combo, size in worst[:30]:
    out.append(f"| {sp:.2f}× | {prec} | {prefix} | {combo} | {size} |")

Path(__file__).parent.joinpath("full-omp1-summary.md").write_text("\n".join(out))
print(f"Wrote summary; {sum(len(v) for v in speedups.values())} total measurements.")
print(f"Routines: {len(speedups)}")
