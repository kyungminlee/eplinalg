"""CLI entry point for the general-purpose migration pipeline.

Usage (from the src/ directory):
    uv run python -m migrator migrate  ../recipes/blas.yaml output/ --target kind16
    uv run python -m migrator verify   output/
    uv run python -m migrator build    ../recipes/blas.yaml output/ --target kind16
    uv run python -m migrator run      ../recipes/blas.yaml work/ --target kind16
    uv run python -m migrator stage    /tmp/staging --target multifloats
"""

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path

from .cmake_gen import _generate_cmake
from .pipeline import (
    run_divergence_report, run_migration,
)
from .prepare import prepare_recipe, run_prepare, verify_patches
from .prefix_classifier import classify_symbols
from .symbol_scanner import scan_symbols
from .cli_common import _get_target_mode, _parser_args
from .fortran.lex import is_comment_line
from .staging import cmd_stage



def cmd_verify_patches(args):
    """Symmetric-patch CI check.

    Exits non-zero with one error per patch whose hunks touch a
    precision-prefixed file without touching all four siblings.
    Allow-list patches that should stay asymmetric in the recipe's
    ``asymmetric_patches:`` field.
    """
    errors = verify_patches(args.recipe, project_root=args.project_root)
    if errors:
        for e in errors:
            print(e, file=sys.stderr)
        return 1
    print(f'{args.recipe.name}: all patches symmetric')
    return 0


def cmd_prepare(args):
    """Stage upstream sources and apply the recipe's patch list.

    Output goes to ``<project_root>/build/staged-sources/<library>/`` and
    is idempotent: a ``.prepared.stamp`` file inside the staged tree
    short-circuits when no listed patch is newer than the stamp. Pass
    ``--rebuild`` to wipe and re-stage.
    """
    staged_root = run_prepare(
        recipe_path=args.recipe,
        project_root=args.project_root,
        rebuild=args.rebuild,
    )
    print(f'Staged: {staged_root}')


def cmd_migrate(args):
    """Run the migration step."""
    parser, parser_cmd = _parser_args(args)
    target = _get_target_mode(args)
    run_migration(
        recipe_path=args.recipe,
        output_dir=args.output_dir,
        target_mode=target,
        dry_run=args.dry_run,
        project_root=args.project_root,
        parser=parser,
        parser_cmd=parser_cmd,
    )


def cmd_diverge(args):
    """Report every co-family pair whose migrated text differs."""
    parser, parser_cmd = _parser_args(args)
    target = _get_target_mode(args)
    report = run_divergence_report(
        recipe_path=args.recipe,
        target_mode=target,
        project_root=args.project_root,
        parser=parser,
        parser_cmd=parser_cmd,
        apply_whitelist=not getattr(args, 'no_whitelist', False),
    )
    total = len(report)
    # Optional filtering on diff content.
    try:
        if args.grep:
            pat = re.compile(args.grep, re.IGNORECASE)
            report = [r for r in report if any(pat.search(l) for l in r['diff'])]
        if args.exclude:
            pat = re.compile(args.exclude, re.IGNORECASE)
            report = [r for r in report if not any(pat.search(l) for l in r['diff'])]
    except re.error as exc:
        print(f'error: invalid regex: {exc}', file=sys.stderr)
        return 2

    for entry in report:
        header = (f'### {entry["other"]} vs {entry["canonical"]}'
                  f' → {entry["target"]} (+{len(entry["diff"])})')
        print(header)
        diff = entry['diff'] if args.full else entry['diff'][:args.context]
        for line in diff:
            print(line[:args.max_width])
        if not args.full and len(entry['diff']) > args.context:
            print(f'  ...{len(entry["diff"]) - args.context} more')
        print()

    shown = len(report)
    if args.grep or args.exclude:
        print(f'{shown} shown / {total} divergent pairs')
    else:
        print(f'{total} divergent pairs')
    return 1 if total else 0


def _is_fixed_form_comment(line: str) -> bool:
    """A fixed-form line is a comment if its first character is C/c/*/!
    OR if its first non-whitespace character is ``!`` (the inline-comment
    marker can also start a whole-line comment when it stands alone).

    Deliberately a superset of :func:`migrator.fortran.lex.is_comment_line`
    (which only tests the column-1 marker): the column-overflow check in
    ``cmd_check`` must also ignore empty lines and indented ``!`` comments,
    which are code-irrelevant but can exceed column 72.
    """
    if not line:
        return True
    if is_comment_line(line):
        return True
    return line.lstrip().startswith('!')


def _is_free_form_comment(line: str) -> bool:
    return not line or line.lstrip().startswith('!')


# Files to exclude from the *standard*-precision sibling archive that
# the embedded Fortran template (and ``cmake/CMakeLists.txt``) builds
# alongside the migrated one. Same EXCLUDE_REGEX semantics as the
# shared ``add_standard_fortran_library`` calls — stems matched
# case-insensitively against the file basename minus its extension.
#
# blas:   ``dsdot`` / ``sdsdot`` are cross-precision (mixed S+D) and
#         not migratable; the migrator also drops them.
# lapack: cross-precision routines (``dsgesv`` / ``zcposv`` / …) and
#         the migrator-introduced extended-precision helpers
#         (``la_constants_ey`` / ``la_xisnan_qx`` / …) — those are
#         migrated content, not upstream content, and don't belong
#         in the standard archive.
_REF_EXCLUDE_STEMS: dict[str, set[str]] = {
    'blas': {'dsdot', 'sdsdot'},
    'lapack': {
        'dsgesv', 'zcgesv', 'dsposv', 'zcposv', 'dsgels', 'zcgels',
        'dlag2s', 'slag2d', 'zlag2c', 'clag2z', 'dlat2s', 'zlat2c',
        'la_constants_ey', 'la_constants_qx', 'la_constants_mw',
        'la_xisnan_ey', 'la_xisnan_qx', 'la_xisnan_mw',
    },
}

# Libraries for which the embedded template emits a standard-precision
# sibling archive. Matches the ``add_standard_fortran_library`` /
# ``add_standard_c_library`` set in ``cmake/CMakeLists.txt``. Other
# Fortran recipes (mumps, scalapack_tools, xblas) intentionally don't
# get a sibling: mumps's sources need extra include directories the
# embedded template doesn't wire; scalapack_tools' three helpers are
# already inside the std scalapack archive; xblas is C-only and the C
# template doesn't ship a sibling.
_REF_LIBRARIES: set[str] = {
    'blas', 'lapack', 'ptzblas', 'pbblas', 'scalapack',
}


def _collect_ref_sources(config) -> list[Path]:
    """Collect upstream Fortran sources for the standard-precision
    sibling archive built alongside the migrated one.

    Globs ``config.source_dir`` for ``.f``/``.F``/``.f90``/``.F90``
    files, pulls in any ``extra_migrate_files`` rooted under
    ``external/`` (so LAPACK picks up ``INSTALL/dlamch.f``), then
    applies the per-library exclude list above.

    Returns ``[]`` for non-Fortran recipes — the C templates don't
    build a sibling archive (yet) — and for recipes outside
    ``_REF_LIBRARIES`` (mumps in particular has Fortran sources that
    INCLUDE per-arithmetic headers from a sibling directory the
    embedded template doesn't wire onto the std-archive's include
    path). Callers treat an empty list as "no standard archive
    emitted".
    """
    if config.language != 'fortran':
        return []
    if config.library not in _REF_LIBRARIES:
        return []
    fortran_exts = {'.f', '.f90'}  # case-insensitive match below
    sources: list[Path] = []
    if config.source_dir.is_dir():
        for f in config.source_dir.iterdir():
            if f.is_file() and f.suffix.lower() in fortran_exts:
                sources.append(f)
    for p in config.extra_migrate_files:
        p = Path(p) if not isinstance(p, Path) else p
        if 'external' in p.parts and p.suffix.lower() in fortran_exts and p.is_file():
            sources.append(p)
    excl = _REF_EXCLUDE_STEMS.get(config.library, set())
    sources = [p for p in sources if p.stem.lower() not in excl]
    # Dedupe + stable order. Resolve to absolute paths so the
    # generated CMakeLists.txt finds them regardless of where the
    # build directory lives.
    seen: dict[tuple, Path] = {}
    for f in sources:
        f = f.resolve()
        try:
            st = f.stat()
            key = (st.st_dev, st.st_ino)
        except OSError:
            key = ('missing', str(f))
        seen.setdefault(key, f)
    return sorted(seen.values())




def cmd_verify(args):
    """Verify migrated output: residual types, literals, column width."""
    out_dir = args.output_dir
    src_dir = out_dir / 'src'
    if not src_dir.is_dir():
        # Fall back to flat layout
        src_dir = out_dir
    errors = 0

    target_mode = _get_target_mode(args)

    f_files = sorted(list(src_dir.glob('*.f')) + list(src_dir.glob('*.F')))
    f90_files = sorted(list(src_dir.glob('*.f90')) + list(src_dir.glob('*.F90')))
    all_files = f_files + f90_files

    print(f'Verifying {len(all_files)} files in {src_dir}')
    print()

    # Read each file once and cache the splitlines result. Without this
    # cache cmd_verify did 3-4 separate read_text passes over every
    # source file in src_dir, dominating wall-time on big libraries.
    file_lines: dict[Path, list[str]] = {
        f: f.read_text(errors='replace').splitlines() for f in all_files
    }

    # Determine if the target uses module-based constructors
    is_constructor_based = target_mode.real_constructor is not None

    # Check residual precision types in code lines
    print('Residual precision types (code lines):')
    residuals = 0
    for f in f_files:
        for i, line in enumerate(file_lines[f], 1):
            if _is_fixed_form_comment(line):
                continue
            if re.search(r'DOUBLE\s+PRECISION|COMPLEX\*16|COMPLEX\*8|DOUBLE\s+COMPLEX|REAL\*[48]',
                         line, re.IGNORECASE):
                print(f'  {f.name}:{i}: {line.strip()}')
                residuals += 1
    for f in f90_files:
        for i, line in enumerate(file_lines[f], 1):
            if _is_free_form_comment(line):
                continue
            if re.search(r'kind\s*\(\s*1\.[de]0\s*\)', line, re.IGNORECASE):
                print(f'  {f.name}:{i}: {line.strip()}')
                residuals += 1
            if re.search(r'\bdouble\s+precision\b', line, re.IGNORECASE):
                print(f'  {f.name}:{i}: {line.strip()}')
                residuals += 1
    if residuals == 0:
        print('  OK')
    else:
        errors += residuals
    print()

    # Build constructor-stripping patterns from target_mode
    def _strip_constructors(line: str) -> str:
        """Remove target-type constructor wrappers from a line."""
        if not target_mode.real_constructor:
            return line
        ctor = re.escape(target_mode.real_constructor)
        line = re.sub(rf"{ctor}\(\s*'[^']*'\s*\)", '', line, flags=re.IGNORECASE)
        line = re.sub(rf'{ctor}\(\s*[^)]*\s*\)', '', line, flags=re.IGNORECASE)
        if target_mode.complex_constructor:
            cctor = re.escape(target_mode.complex_constructor)
            line = re.sub(rf"{cctor}\(\s*'[^']*'\s*\)", '', line, flags=re.IGNORECASE)
            line = re.sub(rf'{cctor}\(\s*[^)]*\s*\)', '', line, flags=re.IGNORECASE)
        return line

    # Check residual D-exponent literals in code lines
    print('Residual D-exponent literals (code lines):')
    d_lits = 0
    for f in f_files:
        for i, line in enumerate(file_lines[f], 1):
            if _is_fixed_form_comment(line):
                continue
            cleaned_line = _strip_constructors(line)
            if re.search(r'[0-9]\.[0-9]*[Dd][+-]?[0-9]', cleaned_line):
                print(f'  {f.name}:{i}: {line.strip()}')
                d_lits += 1
    for f in f90_files:
        for i, line in enumerate(file_lines[f], 1):
            if _is_free_form_comment(line):
                continue
            cleaned_line = _strip_constructors(line)
            # In constructor mode, also reject _wp suffixed literals
            if is_constructor_based and re.search(r'\d+\.\d*_wp', cleaned_line, re.IGNORECASE):
                print(f'  {f.name}:{i}: {line.strip()}')
                d_lits += 1
    if d_lits == 0:
        print('  OK')
    else:
        errors += d_lits
    print()

    if is_constructor_based:
        # Constructor-based target: residual unconverted FP PARAMETER /
        # DATA statements (those that mention a value-shaped numeric and
        # are not commented out).
        print('Residual FP PARAMETER/DATA (code lines):')
        leftover = 0
        for f in all_files:
            is_fixed = f.suffix.lower() == '.f'
            for i, line in enumerate(file_lines[f], 1):
                if is_fixed and _is_fixed_form_comment(line):
                    continue
                if not is_fixed and _is_free_form_comment(line):
                    continue
                m = re.match(r'\s+(PARAMETER|DATA)\b', line, re.IGNORECASE)
                if not m:
                    continue
                cleaned = _strip_constructors(line)
                if re.search(r'\d+\.\d*[DdEe][+-]?\d+|\d*\.\d+', cleaned):
                    print(f'  {f.name}:{i}: {line.strip()}')
                    leftover += 1
        if leftover == 0:
            print('  OK')
        else:
            errors += leftover
        print()

    # Check column width for fixed-form code lines
    print('Column overflow (code lines > 72 chars):')
    overflows = 0
    for f in f_files:
        for i, line in enumerate(file_lines[f], 1):
            if _is_fixed_form_comment(line):
                continue
            if len(line) > 72:
                print(f'  {f.name}:{i}: {len(line)} chars')
                overflows += 1
    if overflows == 0:
        print('  OK')
    else:
        errors += overflows
    print()

    if errors:
        print(f'FAILED: {errors} issue(s)')
        sys.exit(1)
    else:
        print('PASSED')


def cmd_build(args):
    """Generate CMake project and build static libraries."""
    output_dir = args.output_dir
    target_mode = _get_target_mode(args)
    src_dir = output_dir / 'src'
    if not src_dir.is_dir():
        src_dir = output_dir

    config = prepare_recipe(args.recipe, args.project_root)
    lib_name = config.library

    # Classify source files into common vs precision-specific
    symbols = scan_symbols(config.source_dir, config.language,
                           config.extensions,
                           extra_c_return_types=tuple(config.c_return_types))
    classification = classify_symbols(symbols)
    independent = classification.independent

    if config.language == 'c':
        files = sorted(list(src_dir.glob('*.c')))
    else:
        # Honor the recipe's extensions list (normalized to lowercase in
        # load_recipe) so libraries that use .F (MUMPS) or any
        # non-default extension are picked up. Case-insensitive match on
        # the actual filename suffix.
        allowed = {e.lower() for e in config.extensions}
        files = sorted(
            p for p in src_dir.iterdir()
            if p.is_file() and p.suffix.lower() in allowed
        )
    # copy_files are routed to the precision library unconditionally.
    # This single-target build has no per-target file special-casing,
    # and some copy_files entries are target-specific verbatim copies
    # whose bodies USE precision modules (LAPACK's LA_XISNAN_QX USEs
    # LA_CONSTANTS_QX) — a common-lib copy would create a forbidden
    # common -> precision module dependency. NOTE this deliberately
    # differs from ``cmd_stage`` (staging.py), which routes those
    # LA_* pairs via ``_la_own``/``_la_foreign`` first and can then
    # safely send the remaining, genuinely shareable copy_files
    # entries (MUMPS common modules, PBLAS Fortran helpers) to the
    # shared common archive.
    #
    # The LA_CONSTANTS/LA_XISNAN pair filter below IS shared with
    # staging.py: only the target's own suffix pair (_ey/_qx/_mw) may
    # be compiled. Foreign pairs are not just dead weight — the *_mw
    # pair does ``use multifloats``, which doesn't exist in a
    # kind10/kind16 build, so compiling it would break the build.
    _la_suffixes = ('_ey', '_qx', '_mw')
    _own_suffix = target_mode.la_constants_suffix.lower()
    _la_foreign = {
        f'la_{base}{s}'
        for base in ('constants', 'xisnan')
        for s in _la_suffixes
        if s != _own_suffix
    }
    common_files, precision_files = [], []
    for f in files:
        if f.stem in _la_foreign:
            continue
        rel = f.relative_to(output_dir)
        stem = f.stem.upper()
        # ``force_common`` pins a stem to the family-independent archive
        # regardless of scanner assignment (mirror of ``copy_files``,
        # which forces PRECISION). Highest priority. See config.py.
        if stem in config.force_common:
            common_files.append(str(rel))
        elif stem in config.copy_files:
            precision_files.append(str(rel))
        elif stem in independent:
            common_files.append(str(rel))
        else:
            precision_files.append(str(rel))

    ref_sources = _collect_ref_sources(config)

    print(f'Generating CMake project in {output_dir}/')
    print(f'  Common:    {len(common_files)} files')
    print(f'  Precision: {len(precision_files)} files')
    if ref_sources:
        print(f'  Reference: {len(ref_sources)} files (standard-precision sibling)')

    proj_root = (args.project_root or args.recipe.resolve().parent.parent)
    _generate_cmake(output_dir, lib_name, target_mode, common_files, precision_files,
                    language=config.language, project_root=proj_root,
                    ref_sources=ref_sources)

    # Configure and build
    build_dir = output_dir / '_build'
    cmake_cmd = 'cmake'

    # Configure
    configure_args = [
        cmake_cmd, '-S', str(output_dir), '-B', str(build_dir),
        '-DCMAKE_BUILD_TYPE=Release',
    ]
    if config.language != 'c' and args.fc:
        configure_args.append(f'-DCMAKE_Fortran_COMPILER={args.fc}')

    build_dir.mkdir(parents=True, exist_ok=True)
    configure_log = build_dir / 'configure.log'
    build_log = build_dir / 'build.log'

    print('\nConfiguring...')
    with configure_log.open('w') as logf:
        r = subprocess.run(configure_args, stdout=logf,
                           stderr=subprocess.STDOUT, text=True)
    if r.returncode != 0:
        # Tail the log to stderr so the user sees the actual cause.
        log_text = configure_log.read_text(errors='replace')
        tail = log_text.splitlines()[-50:]
        print('\n'.join(tail), file=sys.stderr)
        print(f'\nConfigure failed; full log: {configure_log}', file=sys.stderr)
        sys.exit(1)

    # Build (parallel)
    jobs = os.cpu_count() or 4
    print(f'Building ({jobs} parallel jobs)...')
    with build_log.open('w') as logf:
        r = subprocess.run(
            [cmake_cmd, '--build', str(build_dir), '-j', str(jobs)],
            stdout=logf, stderr=subprocess.STDOUT, text=True,
        )
    if r.returncode != 0:
        log_text = build_log.read_text(errors='replace')
        lines = log_text.splitlines()
        # Surface lines mentioning errors plus the last 50 lines for context.
        err_lines = [l for l in lines
                     if ': error' in l.lower() or 'Error:' in l]
        for l in err_lines[:30]:
            print(f'  {l}')
        print('\n--- last 50 lines of build.log ---', file=sys.stderr)
        print('\n'.join(lines[-50:]), file=sys.stderr)
        print(f'\nBuild failed ({len(err_lines)} error line(s) detected); '
              f'full log: {build_log}', file=sys.stderr)
        sys.exit(1)

    # Report results
    pmap = target_mode.prefix_map
    pair_pfx = f"{pmap['R'].lower()}{pmap['C'].lower()}"
    precision_lib_name = f'lib{pair_pfx}{lib_name}.a'
    common_lib_name = f'lib{lib_name}_common.a'
    ref_lib_name = f'lib{lib_name}.a'

    print('\nBuild succeeded:')
    for name in [common_lib_name, ref_lib_name, precision_lib_name]:
        matches = list(build_dir.rglob(name))
        if matches:
            p = matches[0]
            size = p.stat().st_size
            print(f'  {name}: {size // 1024}K')

    print(f'\nLibraries in {build_dir}/')


def cmd_run(args):
    """Run the full pipeline: migrate → diverge → verify → build."""
    work_dir = args.work_dir
    output_dir = work_dir / 'output'
    src_dir = output_dir / 'src'
    src_dir.mkdir(parents=True, exist_ok=True)

    print('=' * 60)
    print('  Step 1: Migrate')
    print('=' * 60)
    args.output_dir = src_dir
    args.dry_run = False
    cmd_migrate(args)

    print()
    print('=' * 60)
    print('  Step 2: Divergence')
    print('=' * 60)
    args.output_dir = src_dir
    if not hasattr(args, 'grep'):
        args.grep = None
    if not hasattr(args, 'exclude'):
        args.exclude = None
    if not hasattr(args, 'context'):
        args.context = 8
    if not hasattr(args, 'full'):
        args.full = False
    if not hasattr(args, 'max_width'):
        args.max_width = 200
    rc_diverge = cmd_diverge(args) or 0

    print()
    print('=' * 60)
    print('  Step 3: Verify')
    print('=' * 60)
    args.output_dir = output_dir
    rc_verify = 0
    try:
        cmd_verify(args)
    except SystemExit as e:
        if e.code:
            rc_verify = int(e.code) if isinstance(e.code, int) else 1
            print('Verify failed, continuing...')

    print()
    print('=' * 60)
    print('  Step 4: Build (CMake)')
    print('=' * 60)
    args.output_dir = output_dir
    rc_build = cmd_build(args) or 0

    return rc_build or rc_verify or rc_diverge










def _add_parser_args(p):
    """Add --parser and --parser-cmd arguments to a subparser."""
    p.add_argument(
        '--parser', default=None,
        choices=['flang', 'gfortran'],
        help='Parse tree backend for Fortran migration '
             '(default: regex-only, no compiler)')
    p.add_argument(
        '--parser-cmd', default=None,
        help='Explicit path to the parser compiler binary '
             '(overrides PATH lookup)')

def _add_target_args(p):
    p.add_argument('--target', type=str, default='kind16',
                   help='Target name (e.g. "multifloats", "kind16") or path to a target .yaml file')

def build_parser():
    parser = argparse.ArgumentParser(
        prog='migrator',
        description='General-purpose type migration for numerical libraries'
    )
    sub = parser.add_subparsers(dest='command', required=True)

    # --- prepare ---
    p = sub.add_parser(
        'prepare',
        help='Stage upstream sources for a recipe and apply its patch list',
    )
    p.add_argument('recipe', type=Path, help='Recipe YAML file')
    p.add_argument('--project-root', type=Path, default=None)
    p.add_argument('--rebuild', action='store_true',
                   help='Wipe and re-stage even if the cache stamp is fresh')
    p.set_defaults(func=cmd_prepare)

    # --- verify-patches ---
    p = sub.add_parser(
        'verify-patches',
        help='CI check: every patch that touches a precision-prefixed file '
             'must touch all four siblings (or be listed in '
             'asymmetric_patches:)',
    )
    p.add_argument('recipe', type=Path, help='Recipe YAML file')
    p.add_argument('--project-root', type=Path, default=None)
    p.set_defaults(func=cmd_verify_patches)

    # --- migrate ---
    p = sub.add_parser('migrate', help='Migrate source files')
    p.add_argument('recipe', type=Path, help='Recipe YAML file')
    p.add_argument('output_dir', type=Path, help='Output directory')
    _add_target_args(p)
    p.add_argument('--dry-run', action='store_true')
    p.add_argument('--project-root', type=Path, default=None)
    _add_parser_args(p)
    p.set_defaults(func=cmd_migrate)

    # --- diverge ---
    p = sub.add_parser('diverge',
                       help='Report co-family pairs with differing output')
    p.add_argument('recipe', type=Path, help='Recipe YAML file')
    _add_target_args(p)
    p.add_argument('--project-root', type=Path, default=None)
    p.add_argument('--grep', default=None,
                   help='Regex: only show entries with diff matching')
    p.add_argument('--exclude', default=None,
                   help='Regex: drop entries whose diff matches')
    p.add_argument('--context', type=int, default=8,
                   help='Max diff lines per entry (default 8)')
    p.add_argument('--full', action='store_true',
                   help='Print full diff per entry (ignores --context)')
    p.add_argument('--max-width', type=int, default=200,
                   help='Truncate each diff line to this many chars')
    p.add_argument('--no-whitelist', action='store_true',
                   help='Bypass expected_divergences / defer_all_divergences '
                        'whitelist for this run')
    _add_parser_args(p)
    p.set_defaults(func=cmd_diverge)

    # --- verify ---
    p = sub.add_parser('verify', help='Verify migrated output')
    p.add_argument('output_dir', type=Path)
    _add_target_args(p)
    p.set_defaults(func=cmd_verify)

    # --- build ---
    p = sub.add_parser('build', help='Generate CMake project and build')
    p.add_argument('recipe', type=Path, help='Recipe YAML file')
    p.add_argument('output_dir', type=Path)
    _add_target_args(p)
    p.add_argument('--fc', default='gfortran')
    p.add_argument('--project-root', type=Path, default=None)
    p.set_defaults(func=cmd_build)

    # --- run (full pipeline) ---
    p = sub.add_parser('run', help='Run full pipeline')
    p.add_argument('recipe', type=Path, help='Recipe YAML file')
    p.add_argument('work_dir', type=Path, help='Working directory')
    _add_target_args(p)
    p.add_argument('--fc', default='gfortran')
    p.add_argument('--project-root', type=Path, default=None)
    _add_parser_args(p)
    p.set_defaults(func=cmd_run)

    # --- stage (migrate all + unified cmake) ---
    p = sub.add_parser('stage',
                       help='Migrate all libraries into a staging directory '
                            'with a unified CMake build')
    p.add_argument('staging_dir', type=Path,
                   help='Output staging directory')
    _add_target_args(p)
    p.add_argument('--project-root', type=Path, default=None)
    p.add_argument('--libraries', nargs='+', default=None,
                   help='Subset of libraries to migrate (default: all)')
    _add_parser_args(p)
    p.set_defaults(func=cmd_stage)

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()
    rc = args.func(args)
    if isinstance(rc, int) and rc != 0:
        sys.exit(rc)


if __name__ == '__main__':
    main()
