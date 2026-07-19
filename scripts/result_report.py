#!/usr/bin/env python3
"""Generate test/integration/RESULT.md — one giant table per library listing
every migrated double-precision entry and its differential precision
on each target.

Inputs:
  --procedures doc/user/api/procedures.md   canonical migrated-entry list
  --reports DIR [DIR ...]          per-target precision_reports dirs

Output: a Markdown file (or stdout) with sections BLAS / LAPACK /
PBLAS / ScaLAPACK; each row is the double-prefix entry name with one
column per target showing worst-case digits-of-agreement (or "—"
when no test driver exists for that target).
"""

import argparse
import re
from collections import defaultdict
from pathlib import Path

from precision_report import (TARGET_ORDER, fmt_digits, load_reports,
                              order_targets)


LIBRARIES = ['BLAS', 'LAPACK', 'PBLAS', 'ScaLAPACK']
SECTION_RE = re.compile(r'^###\s+(BLAS|LAPACK|PBLAS|ScaLAPACK)\s+—')
ROW_RE = re.compile(r'^\|\s*([^|]+?)\s*\|\s*([^|]+?)\s*\|\s*([^|]+?)\s*\|')

# LAPACK / ScaLAPACK auxiliary heuristic, aligned with the LAPACK
# Users' Guide split (drivers + computational routines vs. auxiliaries):
#   * `xLA*` infix after the precision prefix is always auxiliary.
#   * Names ending in a digit are unblocked variants or recursive
#     helpers (`xGETF2`, `xPOTF2`, `xGEQR2`, `xGEHD2`, `xGEQRT2`,
#     `xPOTRF2`, `xGEQRT3`, …) — auxiliaries by convention. The one
#     widely-documented exception is `qp3` (rank-revealing QR), which
#     is a top-level computational routine.
#   * ScaLAPACK b-prefix entries are back-transformation helpers
#     (PROCEDURES.md) — auxiliary.
AUX_LA_RE = re.compile(r'^[bp]?[diz]+la')
TRAILING_DIGIT_RE = re.compile(r'[0-9]$')
DIGIT_DRIVER_SUFFIXES = ('qp3',)


def is_user_facing(library: str, name: str) -> bool:
    if library in ('BLAS', 'PBLAS'):
        return True
    if library == 'ScaLAPACK' and name.startswith('b'):
        return False
    if AUX_LA_RE.match(name):
        return False
    if TRAILING_DIGIT_RE.search(name):
        return any(name.endswith(s) for s in DIGIT_DRIVER_SUFFIXES)
    return True


def parse_procedures(md_path: Path) -> dict[str, list[str]]:
    """Return {library: [double_entry_name, ...]} preserving uniqueness."""
    by_lib: dict[str, list[str]] = {lib: [] for lib in LIBRARIES}
    seen: dict[str, set[str]] = {lib: set() for lib in LIBRARIES}
    current: str | None = None
    for raw in md_path.read_text().splitlines():
        m = SECTION_RE.match(raw)
        if m:
            current = m.group(1)
            continue
        if current is None or not raw.startswith('|'):
            continue
        m = ROW_RE.match(raw)
        if not m:
            continue
        stem, _single, double = m.group(1), m.group(2), m.group(3)
        # Skip header / separator rows
        if stem in ('stem', '------') or set(stem) <= {'-', ' '}:
            continue
        name = double.strip().strip('`').strip()
        if not name or name == '—':
            continue
        if name not in seen[current]:
            seen[current].add(name)
            by_lib[current].append(name)
    for lib in by_lib:
        by_lib[lib].sort()
    return by_lib


def worst_digits(report: dict) -> tuple[str | None, bool]:
    """Return (formatted-cell, all-passed). cell is fmt_digits()."""
    cases = report.get('cases') or []
    if not cases:
        return None, True
    worst = max(cases, key=lambda c: c['max_rel_err'])
    return fmt_digits(worst.get('digits')), all(c['passed'] for c in cases)


def render_library(name: str, entries: list[str],
                   by_routine: dict[str, dict[str, dict]],
                   targets: list[str]) -> str:
    lines = [f'## {name}', '']
    n_total = len(entries)
    n_tested = sum(1 for e in entries if e in by_routine)
    n_user = sum(1 for e in entries if is_user_facing(name, e))
    lines.append(f'_{n_tested} of {n_total} migrated entries have a '
                 f'dedicated test driver. {n_user} are user-facing '
                 f'(✓ in the first column); the rest are internal '
                 f'auxiliaries called only from user-facing drivers._')
    lines.append('')
    header = '| user | entry |' + ''.join(f' {t} |' for t in targets)
    sep    = '|------|-------|' + ''.join('--------|' for _ in targets)
    lines.append(header)
    lines.append(sep)
    for entry in entries:
        flag = '✓' if is_user_facing(name, entry) else ''
        cells = []
        per_t = by_routine.get(entry, {})
        for t in targets:
            r = per_t.get(t)
            if not r:
                cells.append('—')
                continue
            cell, ok = worst_digits(r)
            if cell is None:
                cells.append('—')
            else:
                cells.append(cell + (' ✗' if not ok else ''))
        lines.append(f'| {flag} | {entry} |' + ''.join(f' {c} |' for c in cells))
    lines.append('')
    return '\n'.join(lines)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--procedures', type=Path, required=True,
                    help='Path to doc/user/api/procedures.md')
    ap.add_argument('--reports', type=Path, nargs='+', required=True,
                    help='One or more precision_reports directories')
    ap.add_argument('--out', type=Path, default=None,
                    help='Output file (default: stdout)')
    args = ap.parse_args()

    by_lib = parse_procedures(args.procedures)

    reports: list[dict] = []
    for d in args.reports:
        if not d.is_dir():
            ap.error(f'{d} is not a directory')
        reports.extend(load_reports(d))

    by_routine: dict[str, dict[str, dict]] = defaultdict(dict)
    for r in reports:
        by_routine[r['routine']][r['target']] = r

    targets = order_targets({r['target'] for r in reports}) \
        if reports else list(TARGET_ORDER)

    out: list[str] = []
    out.append('# Per-routine precision results')
    out.append('')
    out.append(
        'One row per migrated double-precision entry, one column per '
        'target. Cells show worst-case digits-of-agreement against the '
        'quad-precision Netlib reference across all tested problem '
        'shapes (`exact` = bit-identical, `—` = no test driver for '
        'that target). kind16 is the differential reference precision, '
        'so `exact` there is expected. The first column marks '
        'user-facing entries (✓) vs. internal LAPACK / ScaLAPACK '
        'auxiliaries (`xLA*` / `pxLA*` / `bx*` — exercised only '
        'transitively through the drivers above).'
    )
    out.append('')
    out.append('Generated by `scripts/result_report.py` from '
               '`doc/user/api/procedures.md` plus the per-target '
               '`precision_reports/*.json` files. See `test/integration/REPORT.md` '
               'for status / methodology / coverage narrative.')
    out.append('')

    for lib in LIBRARIES:
        out.append(render_library(lib, by_lib[lib], by_routine, targets))

    md = '\n'.join(out)
    if args.out:
        args.out.write_text(md)
        print(f'wrote {args.out}')
    else:
        print(md)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
