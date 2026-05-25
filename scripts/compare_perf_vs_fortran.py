#!/usr/bin/env python3
"""Compare C-harness perf sweep against prior Fortran bench numbers.

Both harnesses measure the same comparison — parallel-blas overlay GF/s
vs migrated Fortran reference GF/s — so the ratios should agree. They
sometimes don't (e.g. Addendum 14 — iTLB churn from oversized parallel-blas
overlay symbols artificially inflates the gap in the Fortran bench).

Reads:
  reports/perf_sweep.tsv     — newer C harness (kernel-isolated)
  reports/full-omp1/*.txt    — older Fortran bench harness
  reports/l2/*.txt           — older Fortran bench harness (L2 subset)
  reports/l3-other/*.txt     — older Fortran bench harness (L3 misc)
  reports/gemm-only/*.txt    — older Fortran bench harness (gemm-only)

Reports cells where the two harnesses' parallel-blas / migrated ratio
diverges by more than 5%, with focus on sub-1.0× cells that flip to parity.

Usage: scripts/compare_perf_vs_fortran.py [--out reports/perf_vs_fortran.md]
"""
from __future__ import annotations
import argparse
import csv
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
BENCH_DIR = REPO / 'reports'

TARGET_PREFIX = {'k10': 'kind10', 'k16': 'kind16', 'mf': 'multifloats'}

FNAME_RE = re.compile(r'^(k10|k16|mf)-bench_([a-z0-9]+)\.txt$')

def parse_fortran_bench(path: Path) -> list[dict]:
    """Parse a Fortran bench .txt file.

    Format: header line "trans  size  subject GFLOP/s  migrated GFLOP/s  speedup"
            (older .txt captures may say "overlay GFLOP/s" instead of
             "subject GFLOP/s" — both mean the same thing; the header
             line is skipped by parsing logic either way)
            then rows: "<key>  <size>  <subject_GFs>  <mig_GFs>  <ratio>x"
    Some files use different first-column names (uplo, side+uplo, etc.) —
    we just take the first whitespace-delimited token as the key.
    """
    rows = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split()
        if len(parts) < 4:
            continue
        # Look for an "Nx" ratio in the last column.
        last = parts[-1]
        if not last.endswith('x'):
            continue
        try:
            ratio = float(last[:-1])
        except ValueError:
            continue
        # Try to back out: key, size, ov, mg, ratio
        try:
            size = int(parts[1])
        except ValueError:
            continue
        try:
            ov = float(parts[2])
            mg = float(parts[3])
        except ValueError:
            continue
        rows.append({'key': parts[0], 'size': size,
                     'parallel_blas': ov, 'migrated': mg, 'ratio': ratio})
    return rows

def load_fortran_data():
    """Walk fortran bench dirs and yield (target, routine, [rows])."""
    data = defaultdict(list)
    for sub in ('full-omp1', 'l2', 'l3-other', 'gemm-only'):
        d = BENCH_DIR / sub
        if not d.is_dir():
            continue
        for f in sorted(d.glob('*.txt')):
            m = FNAME_RE.match(f.name)
            if not m:
                continue
            tgt = TARGET_PREFIX[m.group(1)]
            rname = m.group(2)
            rows = parse_fortran_bench(f)
            # bench_<rname> → parallel-blas overlay routine name, e.g. egemm
            # Sometimes the prefix already matches (k10-bench_egemm → egemm).
            data[(tgt, rname)].extend(rows)
    return data

def load_c_harness_data(tsv: Path):
    """Return {(target, routine): [rows]} from the perf sweep TSV."""
    data = defaultdict(list)
    with tsv.open() as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for r in reader:
            try:
                size = int(r['size'])
                ov = float(r['parallel_blas_GFs'])
                mg = float(r['migrated_GFs'])
                ratio = float(r['ratio'])
            except (KeyError, ValueError):
                continue
            data[(r['target'], r['routine'])].append({
                'key': r['key'], 'size': size,
                'parallel_blas': ov, 'migrated': mg, 'ratio': ratio,
            })
    return data

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tsv', default='reports/perf_sweep.tsv')
    ap.add_argument('--out', default='reports/perf_vs_fortran.md')
    args = ap.parse_args()

    fortran = load_fortran_data()
    c_data = load_c_harness_data(Path(args.tsv))

    out = []
    out.append('# parallel-blas / migrated ratio — C harness vs Fortran harness (OMP=1)')
    out.append('')
    out.append('Both harnesses measure the same comparison — parallel-blas overlay GF/s '
               'vs migrated Fortran reference GF/s. Per-routine: median ratio across all '
               '(key, size) cells, reported by both harnesses.')
    out.append('')
    out.append('| target | routine | Fortran-harness median | C-harness median | F #cells | C #cells | Δ |')
    out.append('|--------|---------|----------------------:|----------------:|---------:|---------:|--:|')

    keys_union = sorted(set(fortran) | set(c_data))
    for (tgt, rname) in keys_union:
        f_rows = fortran.get((tgt, rname), [])
        c_rows = c_data.get((tgt, rname), [])
        if not f_rows and not c_rows:
            continue
        f_med = statistics.median(r['ratio'] for r in f_rows) if f_rows else None
        c_med = statistics.median(r['ratio'] for r in c_rows) if c_rows else None
        mark = ''
        if f_med is not None and c_med is not None:
            delta = c_med - f_med
            if abs(delta) > 0.05:
                mark = ' ⚠' if abs(delta) > 0.15 else ' ·'
            delta_str = f'{delta:+.3f}{mark}'
        else:
            delta_str = '—'
        out.append(
            f'| {tgt} | `{rname}` | '
            f'{(f_med and f"{f_med:.3f}") or "—"} | '
            f'{(c_med and f"{c_med:.3f}") or "—"} | '
            f'{len(f_rows)} | {len(c_rows)} | {delta_str} |'
        )

    out.append('')
    out.append('## Cells with biggest divergence (|Δ ratio| > 0.10)')
    out.append('')
    out.append('Same (target, routine, key, size) tuple in both datasets. Both ratios are parallel-blas / migrated.')
    out.append('')
    out.append('| target | routine | key | size | Fortran-harness ratio | C-harness ratio | Δ |')
    out.append('|--------|---------|-----|------|---------------------:|---------------:|--:|')
    divergence = []
    for (tgt, rname), f_rows in fortran.items():
        c_rows = c_data.get((tgt, rname), [])
        if not c_rows:
            continue
        c_lookup = {(r['key'], r['size']): r['ratio'] for r in c_rows}
        for r in f_rows:
            key = (r['key'], r['size'])
            if key in c_lookup:
                d = c_lookup[key] - r['ratio']
                if abs(d) > 0.10:
                    divergence.append((tgt, rname, *key, r['ratio'], c_lookup[key], d))
    divergence.sort(key=lambda x: -abs(x[-1]))
    for tgt, rname, k, sz, fr, cr, d in divergence[:40]:
        out.append(f'| {tgt} | `{rname}` | `{k}` | {sz} | {fr:.3f} | {cr:.3f} | {d:+.3f} |')
    out.append(f'')
    out.append(f'(Top 40 of {len(divergence)} divergent cells.)')

    Path(args.out).write_text('\n'.join(out) + '\n')
    print(f'wrote {args.out}: {len(keys_union)} routine pairs')

if __name__ == '__main__':
    main()
