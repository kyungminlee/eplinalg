"""Extend the MUMPS keep-kind manifest with call-site entries.

The base sweep (``mumps_sweep_keep_kind.py``) protects DP declarations
in paired s/d and c/z files. This complementary sweep handles the
other half of the problem: lines that *call* a shared (copy_files)
routine with DP arguments.

Why: shared modules like ``mumps_load.F`` / ``mumps_comm_buffer_common.F``
declare DP parameters and are copied verbatim (so they keep their
DP interface across all target builds). Their callers, however, are
migrated arithmetic-specific files. If a caller uses ``dble(x)`` to
build the DP argument, the normal migrator rewrites it to
``REAL(x, KIND=16)`` — which mismatches the DP (REAL(8)) parameter.

Fix: append a manifest entry for every line (including continuation
lines) that forms part of a ``CALL <shared_routine>(...)`` statement
or a reference to a shared function. The keep-kind sentinel then
shields ``dble(`` / ``dcmplx(`` on that line from the intrinsic
rewrite. ``DOUBLE PRECISION`` declarations on the same line are also
shielded as a side-effect (harmless — those lines tend not to contain
DP declarations).

Usage:
    python mumps_sweep_keep_kind_calls.py [--append] [--verbose]

By default, prints additions to stdout. With ``--append``, merges them
into ``recipes/mumps/keep-kind.manifest`` (de-duplicating against the
existing entries).
"""

import argparse
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SRC_DIR = REPO / 'external' / 'MUMPS_5.9.0' / 'src'
MANIFEST = REPO / 'recipes' / 'mumps' / 'keep-kind.manifest'

# copy_files stems (must match recipes/mumps.yaml). Anything listed
# here contributes its SUBROUTINE / FUNCTION names to the "DP-stable"
# set; callers of those names get their lines marked keep-kind.
def _load_copy_files() -> set[str]:
    """Read the mumps.yaml recipe and return its copy_files stems."""
    try:
        import yaml  # type: ignore
    except ImportError:
        yaml = None
    recipe = REPO / 'recipes' / 'mumps.yaml'
    if yaml is None or not recipe.exists():
        # Fallback for environments without PyYAML: return nothing so
        # the sweep exits cleanly. The main migration is the real source
        # of truth.
        return set()
    data = yaml.safe_load(recipe.read_text())
    return set(str(s).upper() for s in (data.get('copy_files') or []))


COPY_FILES = _load_copy_files()

_SUB_RE = re.compile(
    r'^\s+(?:RECURSIVE\s+|PURE\s+|ELEMENTAL\s+)*'
    r'(?:SUBROUTINE|FUNCTION)\s+([A-Z_][A-Z0-9_]*)\s*\(',
    re.IGNORECASE,
)


def collect_shared_routines() -> set[str]:
    """Parse each copy_files source file and extract its top-level
    SUBROUTINE / FUNCTION names."""
    out: set[str] = set()
    for stem in COPY_FILES:
        # Try .F first (all MUMPS Fortran sources use .F).
        paths = [SRC_DIR / f'{stem.lower()}.F']
        for p in paths:
            if not p.exists():
                continue
            for raw in p.read_text(errors='replace').splitlines():
                if not raw or raw[0] in ('C', 'c', '*', '!'):
                    continue
                m = _SUB_RE.match(raw)
                if m:
                    out.add(m.group(1).upper())
    return out


def continuation_ranges(lines: list[str]) -> list[tuple[int, int]]:
    """Return a list of (start, end) 1-based line ranges where each
    range corresponds to a single logical statement (starting a
    statement, plus any continuation lines marked with ``&`` in col 6)."""
    ranges: list[tuple[int, int]] = []
    i = 0
    n = len(lines)
    while i < n:
        start = i
        # Advance through following continuation lines: in fixed-form
        # Fortran, a continuation line has a non-blank, non-zero char in
        # column 6 (1-based). The first line of a statement has a blank
        # there.
        end = i
        while end + 1 < n:
            nxt = lines[end + 1]
            if (len(nxt) >= 6 and nxt[5] not in (' ', '0', '\t')
                    and (not nxt or nxt[0] not in ('C', 'c', '*'))):
                end += 1
            else:
                break
        ranges.append((start + 1, end + 1))  # 1-based
        i = end + 1
    return ranges


def find_call_sites(shared: set[str]) -> set[tuple[str, int]]:
    """Scan every .F source under SRC_DIR and return ``(basename, lineno)``
    pairs for each line that is part of a CALL / bare-reference to a
    shared routine."""
    out: set[tuple[str, int]] = set()
    # Word-boundary regex matching any of the shared routine names.
    if not shared:
        return out
    name_alt = '|'.join(sorted(re.escape(n) for n in shared))
    ref_re = re.compile(r'\b(?:' + name_alt + r')\b', re.IGNORECASE)

    for src in sorted(SRC_DIR.glob('*.F')):
        text = src.read_text(errors='replace')
        lines = text.splitlines()
        ranges = continuation_ranges(lines)
        for (s, e) in ranges:
            stmt = '\n'.join(lines[s - 1:e])
            # Exclude the SUBROUTINE/FUNCTION definition itself.
            first = lines[s - 1]
            if _SUB_RE.match(first):
                continue
            # Exclude comment-only / preprocessor lines.
            if first and first[0] in ('C', 'c', '*', '!', '#'):
                continue
            if not ref_re.search(stmt):
                continue
            # Only mark lines that actually contain tokens we would
            # want to shield. Lines that reference a shared routine
            # but carry no ``dble(`` / ``dcmplx(`` are noise — the
            # migrator wouldn't touch them anyway.
            for ln in range(s, e + 1):
                body = lines[ln - 1]
                if re.search(r'\bdble\s*\(', body, re.IGNORECASE) \
                        or re.search(r'\bdcmplx\s*\(', body,
                                     re.IGNORECASE):
                    out.add((src.name, ln))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--append', action='store_true',
                    help='Merge into the existing manifest (default: stdout only)')
    ap.add_argument('--verbose', action='store_true')
    args = ap.parse_args()

    shared = collect_shared_routines()
    if args.verbose:
        print(f'# {len(shared)} shared routines:')
        for n in sorted(shared):
            print(f'#   {n}')

    sites = find_call_sites(shared)
    entries = sorted(
        f'external/MUMPS_5.9.0/src/{name}:{ln}'
        for (name, ln) in sites
    )

    if not args.append:
        for e in entries:
            print(e)
        return 0

    existing = set()
    if MANIFEST.exists():
        for line in MANIFEST.read_text().splitlines():
            t = line.strip()
            if t and not t.startswith('#'):
                existing.add(t)

    added = [e for e in entries if e not in existing]
    if not added:
        print(f'# No new entries to add (manifest already has {len(existing)})')
        return 0

    with MANIFEST.open('a') as f:
        f.write('\n# Call-site entries (added by mumps_sweep_keep_kind_calls.py)\n')
        for e in added:
            f.write(e + '\n')
    print(f'# Added {len(added)} entries to {MANIFEST.relative_to(REPO)}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
