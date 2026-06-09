#!/usr/bin/env python3
"""Aggregate per-test precision reports into a single Markdown summary.

Each test program writes <routine>.<target>.json into the reports
directory at run time. This script reads all of them and emits two
tables:

  1. Per-routine agreement (rows = routines, cols = targets), showing
     the worst-case digits-of-agreement across that routine's cases.
  2. Per-case detail under each routine (collapsible header per
     routine), useful for finding which problem shape regressed.

Usage:
    python scripts/precision_report.py <reports-dir> [--out <file>]

Reads every *.json in <reports-dir>; writes Markdown to stdout (or
<file> if --out is given).
"""

import argparse
import json
from collections import defaultdict
from pathlib import Path
from typing import Iterable


# Stable, sensible target column order — matches the build's prefix
# scheme. Anything unexpected is appended at the end alphabetically.
TARGET_ORDER = ['kind4', 'kind8', 'kind10', 'kind16', 'multifloats']


def load_reports(reports_dir: Path) -> list[dict]:
    """Read every *.json report and tag each with its (routine, target)."""
    out = []
    for p in sorted(reports_dir.glob('*.json')):
        try:
            data = json.loads(p.read_text())
        except json.JSONDecodeError as exc:
            print(f'# WARNING: failed to parse {p}: {exc}')
            continue
        if 'routine' not in data or 'target' not in data:
            print(f'# WARNING: {p} missing routine/target keys')
            continue
        if not isinstance(data.get('cases'), list):
            print(f'# WARNING: {p} missing or non-list "cases" — skipping')
            continue
        out.append(data)
    return out


def order_targets(targets: Iterable[str]) -> list[str]:
    """Sort with TARGET_ORDER first, then alphabetical for unknowns."""
    seen = set(targets)
    head = [t for t in TARGET_ORDER if t in seen]
    tail = sorted(seen - set(TARGET_ORDER))
    return head + tail


def fmt_digits(d: float | None) -> str:
    if d is None:
        return '—'
    if d >= 90:
        return 'exact'  # 99.00 sentinel from the Fortran writer
    return f'{d:.2f}'


def fmt_err(e: float | None) -> str:
    if e is None:
        return '—'
    if e == 0.0:
        return '0'
    return f'{e:.2e}'


def render_summary(reports: list[dict]) -> str:
    """Per-routine worst-case digits across targets."""
    by_routine: dict[str, dict[str, dict]] = defaultdict(dict)
    for r in reports:
        by_routine[r['routine']][r['target']] = r

    targets = order_targets({r['target'] for r in reports})
    if not targets:
        return '_(no reports)_\n'

    lines = []
    lines.append('## Per-routine summary')
    lines.append('')
    lines.append(
        'Worst-case agreement (in decimal digits) across all problem '
        'shapes per routine, per target. "exact" means bit-identical '
        'agreement (max_rel_err = 0) — typical for kind16 since both '
        'reference and target are REAL(KIND=16). A `✗` marks a failed '
        'case (max_rel_err exceeded its tolerance).'
    )
    lines.append('')
    header = '| routine |' + ''.join(f' {t} |' for t in targets)
    sep    = '|---------|' + ''.join('---|' for _ in targets)
    lines.append(header)
    lines.append(sep)

    for routine in sorted(by_routine):
        cells = []
        for t in targets:
            r = by_routine[routine].get(t)
            if not r or not r['cases']:
                cells.append('—')
                continue
            worst_case = max(r['cases'], key=lambda c: c['max_rel_err'])
            mark = ' ✗' if not worst_case['passed'] else ''
            cells.append(fmt_digits(worst_case.get('digits')) + mark)
        lines.append(f'| {routine} |' + ''.join(f' {c} |' for c in cells))

    lines.append('')
    return '\n'.join(lines)


def render_detail(reports: list[dict]) -> str:
    """Per-case detail tables, one section per (routine, target)."""
    by_routine: dict[str, dict[str, dict]] = defaultdict(dict)
    for r in reports:
        by_routine[r['routine']][r['target']] = r

    targets = order_targets({r['target'] for r in reports})

    lines = ['## Per-case detail', '']
    for routine in sorted(by_routine):
        lines.append(f'### {routine}')
        lines.append('')
        lines.append('| target | case | max_rel_err | tolerance | digits | pass |')
        lines.append('|--------|------|-------------|-----------|--------|------|')
        for t in targets:
            r = by_routine[routine].get(t)
            if not r:
                continue
            for c in r['cases']:
                pass_mark = '✓' if c['passed'] else '✗'
                lines.append(
                    f'| {t} | {c["case"]} | {fmt_err(c["max_rel_err"])} '
                    f'| {fmt_err(c["tolerance"])} '
                    f'| {fmt_digits(c.get("digits"))} | {pass_mark} |'
                )
        lines.append('')
    return '\n'.join(lines)


def render_header(reports: list[dict]) -> str:
    targets = order_targets({r['target'] for r in reports})
    routines = sorted({r['routine'] for r in reports})
    n_cases = sum(len(r['cases']) for r in reports)
    n_failed = sum(
        1 for r in reports for c in r['cases'] if not c['passed']
    )
    status = 'PASS' if n_failed == 0 else f'{n_failed} FAILED CASES'

    lines = ['# BLAS differential precision report', '']
    lines.append(f'**Status:** {status}')
    lines.append('')
    lines.append(
        f'**Coverage:** {len(routines)} routines × {len(targets)} targets '
        f'= {len(reports)} report files, {n_cases} test cases.'
    )
    lines.append(f'**Targets:** {", ".join(targets)}')
    lines.append('')
    lines.append(
        'Each test runs the same operation through the migrated '
        'extended-precision routine and through a REAL(KIND=16) '
        'reference (Netlib BLAS compiled with `-freal-8-real-16`). '
        'Results are compared in REAL(KIND=16); the table reports the '
        'worst-case relative error converted to decimal digits of '
        'agreement.'
    )
    lines.append('')
    return '\n'.join(lines)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('reports_dir', type=Path,
                    help='Directory containing *.json reports')
    ap.add_argument('--out', type=Path, default=None,
                    help='Output file (default: stdout)')
    args = ap.parse_args()

    if not args.reports_dir.is_dir():
        ap.error(f'{args.reports_dir} is not a directory')

    reports = load_reports(args.reports_dir)
    if not reports:
        ap.error(f'no JSON reports found in {args.reports_dir}')

    md = '\n'.join([
        render_header(reports),
        render_summary(reports),
        render_detail(reports),
    ])

    if args.out:
        args.out.write_text(md)
        print(f'wrote {args.out}')
    else:
        print(md)

    n_failed = sum(
        1 for r in reports for c in r['cases'] if not c['passed']
    )
    return 1 if n_failed else 0


if __name__ == '__main__':
    raise SystemExit(main())
