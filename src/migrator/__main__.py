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
import shutil
import subprocess
import sys
from pathlib import Path

from .pipeline import (
    run_divergence_report, run_migration,
)
from .prepare import prepare_recipe, run_prepare, verify_patches
from .prefix_classifier import classify_symbols
from .symbol_scanner import scan_symbols
from .target_mode import load_target

# BLACS-style dual-entry-point detector used by ``cmd_stage`` to
# identify C sources that switch their public symbol via the
# ``INTFACE == C_CALL`` / ``CallFromC`` macros. Hoisted to module scope
# so it's compiled once instead of per-library in the cmd_stage loop.
_DUAL_ENTRY_C_RE = re.compile(
    r'#\s*if\s*\(?\s*INTFACE\s*==\s*C_CALL\b'
    r'|#\s*ifdef\s+CallFromC\b'
    r'|#\s*if\s+defined\s*\(\s*CallFromC\s*\)',
)


def _get_target_mode(args):
    """Construct TargetMode based on CLI arguments."""
    target_str = getattr(args, 'target', None) or 'kind16'
    return load_target(target_str)

def _parser_args(args):
    """Extract parser/parser_cmd from CLI args."""
    parser = getattr(args, 'parser', None)
    parser_cmd = getattr(args, 'parser_cmd', None)
    return parser, parser_cmd


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
    marker can also start a whole-line comment when it stands alone)."""
    if not line:
        return True
    if line[0] in ('C', 'c', '*', '!'):
        return True
    stripped = line.lstrip()
    return stripped.startswith('!')


def _is_free_form_comment(line: str) -> bool:
    return not line or line.lstrip().startswith('!')


def _patch_libseq_mpi_f(path: Path) -> None:
    """Extend libseq's ``MUMPS_COPY`` with MPI_REAL16 / MPI_COMPLEX32
    cases so reductions on REAL(KIND=16) / COMPLEX(KIND=16) buffers
    dispatch correctly under our libmpiseq variant. Upstream's
    ``MUMPS_COPY`` only knows the standard MPI datatypes; the migrated
    qmumps archive passes MPI_REAL16 (Intel MPI = 1275072555) for
    kind16 reductions.

    Patches the staged copy at ``_mpiseq_src/mpi.f``; upstream's
    ``external/MUMPS_5.8.2/libseq/mpi.f`` stays read-only. BLACS /
    ScaLAPACK forwarders inside the same file are deliberately KEPT
    — libmpiseq stands in for those archives in the ``_seq`` test
    link, and the real BLACS / ScaLAPACK archives aren't linked there
    so there's no duplicate-symbol collision.
    """
    src = path.read_text()

    # Extend MUMPS_COPY's dispatch with MPI_REAL16 / MPI_COMPLEX32
    # cases, and append matching MUMPS_COPY_* helpers. Anchor the
    # insertion before the existing ``ELSE\n IERR=1`` fallthrough.
    fallthrough = '      ELSE\n        IERR=1\n        RETURN\n      END IF'
    extra_dispatch = (
        # kind16: REAL(16) / COMPLEX(16)
        '      ELSE IF ( DATATYPE .EQ. MPI_REAL16 ) THEN\n'
        '      CALL MUMPS_COPY_REAL16( SENDBUF, RECVBUF, CNT, SS, RS )\n'
        '      ELSE IF ( DATATYPE .EQ. MPI_COMPLEX32 ) THEN\n'
        '      CALL MUMPS_COPY_COMPLEX32( SENDBUF, RECVBUF, CNT, SS, RS )\n'
        # kind10: 80-bit extended real / complex map to MPI's long
        # double tokens (no MPI_REAL10 in standard MPI).
        '      ELSE IF ( DATATYPE .EQ. MPI_LONG_DOUBLE ) THEN\n'
        '      CALL MUMPS_COPY_REAL10( SENDBUF, RECVBUF, CNT, SS, RS )\n'
        '      ELSE IF ( DATATYPE .EQ. MPI_C_LONG_DOUBLE_COMPLEX ) THEN\n'
        '      CALL MUMPS_COPY_COMPLEX20( SENDBUF, RECVBUF, CNT, SS, RS )\n'
        # multifloats: cmake/mpiseq_c_stubs.c encodes derived-type
        # sentinels as 0x10000000 | total_bytes. float64x2 → 16-byte
        # element (sentinel 268435472); complex64x2 → 32-byte
        # (268435488). MPI_Type_c2f / MPI_Op_c2f in Intel mpi.h are
        # passthrough casts, so the Fortran handle is the same value.
        '      ELSE IF ( DATATYPE .EQ. 268435472 ) THEN\n'
        '      CALL MUMPS_COPY_FLOAT64X2( SENDBUF, RECVBUF, CNT, SS, RS )\n'
        '      ELSE IF ( DATATYPE .EQ. 268435488 ) THEN\n'
        '      CALL MUMPS_COPY_COMPLEX64X2( SENDBUF, RECVBUF, CNT, SS, RS )\n'
    )
    if 'MPI_REAL16' not in src and fallthrough in src:
        src = src.replace(fallthrough, extra_dispatch + fallthrough, 1)

    extra_helpers = """
      SUBROUTINE MUMPS_COPY_REAL16( S, R, N, SS, RS )
      IMPLICIT NONE
      INTEGER N, SS, RS
      REAL(KIND=16) S(N),R(N)
      INTEGER I
      DO I = 1, N
        R(I+RS) = S(I+SS)
      END DO
      RETURN
      END SUBROUTINE MUMPS_COPY_REAL16
      SUBROUTINE MUMPS_COPY_COMPLEX32( S, R, N, SS, RS )
      IMPLICIT NONE
      INTEGER N, SS, RS
      COMPLEX(KIND=16) S(N),R(N)
      INTEGER I
      DO I = 1, N
        R(I+RS) = S(I+SS)
      END DO
      RETURN
      END SUBROUTINE MUMPS_COPY_COMPLEX32
      SUBROUTINE MUMPS_COPY_REAL10( S, R, N, SS, RS )
      IMPLICIT NONE
      INTEGER N, SS, RS
      REAL(KIND=10) S(N),R(N)
      INTEGER I
      DO I = 1, N
        R(I+RS) = S(I+SS)
      END DO
      RETURN
      END SUBROUTINE MUMPS_COPY_REAL10
      SUBROUTINE MUMPS_COPY_COMPLEX20( S, R, N, SS, RS )
      IMPLICIT NONE
      INTEGER N, SS, RS
      COMPLEX(KIND=10) S(N),R(N)
      INTEGER I
      DO I = 1, N
        R(I+RS) = S(I+SS)
      END DO
      RETURN
      END SUBROUTINE MUMPS_COPY_COMPLEX20
      SUBROUTINE MUMPS_COPY_FLOAT64X2( S, R, N, SS, RS )
      IMPLICIT NONE
      INTEGER N, SS, RS
      REAL(KIND=8) S(2*N),R(2*N)
      INTEGER I
      DO I = 1, 2*N
        R(I+2*RS) = S(I+2*SS)
      END DO
      RETURN
      END SUBROUTINE MUMPS_COPY_FLOAT64X2
      SUBROUTINE MUMPS_COPY_COMPLEX64X2( S, R, N, SS, RS )
      IMPLICIT NONE
      INTEGER N, SS, RS
      REAL(KIND=8) S(4*N),R(4*N)
      INTEGER I
      DO I = 1, 4*N
        R(I+4*RS) = S(I+4*SS)
      END DO
      RETURN
      END SUBROUTINE MUMPS_COPY_COMPLEX64X2
"""
    if 'SUBROUTINE MUMPS_COPY_REAL16' not in src:
        src = src.rstrip() + '\n' + extra_helpers

    path.write_text(src)


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
#         (``la_constants_ep`` / ``la_xisnan_mf`` / …) — those are
#         migrated content, not upstream content, and don't belong
#         in the standard archive.
_REF_EXCLUDE_STEMS: dict[str, set[str]] = {
    'blas': {'dsdot', 'sdsdot'},
    'lapack': {
        'dsgesv', 'zcgesv', 'dsposv', 'zcposv', 'dsgels', 'zcgels',
        'dlag2s', 'slag2d', 'zlag2c', 'clag2z', 'dlat2s', 'zlat2c',
        'la_constants_ep', 'la_constants_mf',
        'la_xisnan_ep', 'la_xisnan_mf',
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


def _collect_source_files(src_dir: Path, language: str) -> list[Path]:
    """Discover migrated source files in ``src_dir`` for the given language.

    Honors all four Fortran extension cases (``.f``/``.F``/``.f90``/``.F90``).
    Dedupe uses ``(st_dev, st_ino)`` so case-insensitive filesystems (where
    ``*.f`` and ``*.F`` glob the same physical file) do not double-stage.
    """
    if language == 'c':
        patterns = ('*.c', '*.f', '*.F', '*.f90', '*.F90')
    else:
        patterns = ('*.f', '*.F', '*.f90', '*.F90')
    seen: dict[tuple, Path] = {}
    for pat in patterns:
        for f in src_dir.glob(pat):
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


def _generate_cmake(output_dir: Path, lib_name: str, target_mode,
                    common_files: list[str], precision_files: list[str],
                    language: str = 'fortran',
                    project_root: Path | None = None,
                    ref_sources: list[Path] | None = None):
    """Generate a self-contained CMakeLists.txt in the output directory."""
    pmap = target_mode.prefix_map
    real_pfx = pmap['R'].lower()
    precision_lib = f'{real_pfx}{lib_name}'
    common_lib = f'{lib_name}_common'

    common_list = '\n    '.join(sorted(common_files))
    precision_list = '\n    '.join(sorted(precision_files))
    ref_list = '\n    '.join(f'"{p}"' for p in (ref_sources or []))

    # Default path to the vendored Intel MPI headers. ``project_root``
    # is resolved at generation time, so the generated CMakeLists.txt
    # works when built from a fresh out-of-tree output directory.
    _impi_default = str(((project_root or Path.cwd())
                         / 'external' / 'impi-headers').resolve())

    if language == 'c':
        # When targeting multifloats, the migrated C sources `#include
        # "multifloats_bridge.h"`. Copy + patch the bridge header into
        # the output dir so the migrated sources find it on the include
        # path, mirroring what cmd_stage does for the shared driver.
        c_mf_link = ''
        c_mf_deps = ''
        if target_mode.module_name is not None:
            _root = project_root or Path.cwd()
            mf_local = _root / 'external' / 'multifloats-mpi'
            bridge_h_src = mf_local / 'multifloats_bridge.h'
            if bridge_h_src.is_file():
                helpers_dst = output_dir / '_helpers'
                helpers_dst.mkdir(exist_ok=True)
                staged = helpers_dst / bridge_h_src.name
                shutil.copy2(bridge_h_src, staged)
                # Guard the bridge header's `#include <mpi.h>` against
                # the C++ MPI bindings — without it, scalapack_c's
                # C-as-C++ build pulls thousands of mpicxx.h templates
                # into the migrator-injected ``extern "C" { … }`` wrap
                # and fails to link. Same patch cmd_stage applies.
                text = staged.read_text()
                if 'MPICH_SKIP_MPICXX' not in text:
                    text = text.replace(
                        '#include <mpi.h>',
                        '#define MPICH_SKIP_MPICXX 1\n'
                        '#define OMPI_SKIP_MPICXX 1\n'
                        '#include <mpi.h>',
                        1,
                    )
                    staged.write_text(text)
            c_mf_link = """
# multifloats: FetchContent (or local via -DMULTIFLOATS_DIR) so the
# migrated sources can link against ``libmultifloats.a`` (C++) and
# include ``multifloats_bridge.h`` (staged into ./_helpers/).
set(BUILD_TESTING OFF CACHE BOOL "Disable tests in fetched multifloats" FORCE)
set(MULTIFLOATS_BUILD_BENCH OFF CACHE BOOL "Disable benches in fetched multifloats" FORCE)
if(DEFINED MULTIFLOATS_DIR)
    message(STATUS "Using local multifloats: ${MULTIFLOATS_DIR}")
    add_subdirectory(${MULTIFLOATS_DIR}
        ${CMAKE_CURRENT_BINARY_DIR}/_mf EXCLUDE_FROM_ALL)
else()
    include(FetchContent)
    set(MULTIFLOATS_GIT_REPO "https://github.com/kyungminlee/multifloats.git"
        CACHE STRING "Git URL for the multifloats library")
    set(MULTIFLOATS_GIT_TAG "v0.6.0"
        CACHE STRING "Git tag/branch/commit for multifloats (>= v0.6.0)")
    FetchContent_Declare(multifloats_fetch
        GIT_REPOSITORY ${MULTIFLOATS_GIT_REPO}
        GIT_TAG        ${MULTIFLOATS_GIT_TAG}
    )
    FetchContent_Populate(multifloats_fetch)
    add_subdirectory(
        ${multifloats_fetch_SOURCE_DIR}
        ${CMAKE_CURRENT_BINARY_DIR}/_mf EXCLUDE_FROM_ALL)
endif()
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/_helpers)
"""
            c_mf_deps = f"""
if(TARGET multifloats)
    target_link_libraries({precision_lib} PUBLIC multifloats)
    if(TARGET {common_lib})
        target_link_libraries({common_lib} PUBLIC multifloats)
    endif()
endif()
# multifloats's bridge header (multifloats_bridge.h) uses C++ ``using``
# declarations to expose ``float64x2`` at file scope. Migrated .c
# bodies need a C++ translation unit for those typedefs to resolve.
# Reclassify the migrated sources to LANGUAGE CXX so they go through
# the C++ compiler while keeping their .c extension on disk.
set_source_files_properties(${{PRECISION_SOURCES}} PROPERTIES LANGUAGE CXX)
if(COMMON_SOURCES)
    set_source_files_properties(${{COMMON_SOURCES}} PROPERTIES LANGUAGE CXX)
endif()
# C-as-C++ flags (mirror add_migrated_c_library in the shared driver):
#   - cxx_std_17: multifloats requires it.
#   - -fpermissive: BLACS/PBLAS bodies have implicit ``void *`` →
#     ``char *`` casts (e.g. BI_iMPI_amn.c:12) that C accepts and C++
#     rejects. Tolerated on g++; clang++ does NOT honor this flag,
#     so a clang-based multifloats build still needs the explicit-cast
#     overrides in ``recipes/<lib>/mfc_overrides/``.
#   - -Wno-write-strings: silences `const char[]` -> `char *` literals.
#   - MPICH_SKIP_MPICXX / OMPI_SKIP_MPICXX: keep mpicxx.h's templates
#     out of the migrator-injected ``extern "C" {{ ... }}`` wrap.
target_compile_features({precision_lib} PRIVATE cxx_std_17)
target_compile_options({precision_lib} PRIVATE -fpermissive -Wno-write-strings)
target_compile_definitions({precision_lib} PRIVATE
    MPICH_SKIP_MPICXX OMPI_SKIP_MPICXX)
if(TARGET {common_lib})
    target_compile_features({common_lib} PRIVATE cxx_std_17)
    target_compile_options({common_lib} PRIVATE -fpermissive -Wno-write-strings)
    target_compile_definitions({common_lib} PRIVATE
        MPICH_SKIP_MPICXX OMPI_SKIP_MPICXX)
endif()
"""

        # Add CXX to project() languages when multifloats is in play —
        # multifloats's targets request cxx_std_17 features from the
        # embedding project. Harmless on KIND targets (small detect
        # cost, no sources compiled as C++).
        project_langs = 'C CXX' if target_mode.module_name is not None else 'C'

        cmake = f"""\
cmake_minimum_required(VERSION 3.20)
project({precision_lib} {project_langs})

# --- Compiler flags ---
set(CMAKE_C_FLAGS "${{CMAKE_C_FLAGS}} -w")

# --- MPI: default to vendored Intel MPI headers ---
# ``external/impi-headers/`` ships mpi.h and mpif.h at the Intel MPI
# ABI so the build compiles against a stable MPI surface without
# requiring every contributor to install an MPI runtime. Link-time
# libraries still come from whichever MPI runtime the user provides
# (impi-rt / OpenMPI / MPICH — headers are ABI-compatible).
# Users who want a different MPI's *headers* can override IMPI_HEADERS.
if(NOT DEFINED IMPI_HEADERS)
    set(IMPI_HEADERS "{_impi_default}"
        CACHE PATH "Path to vendored Intel MPI headers")
endif()
include_directories(${{IMPI_HEADERS}})
find_package(MPI COMPONENTS C QUIET)
{c_mf_link}
# --- Common (type-independent) library ---
set(COMMON_SOURCES
    {common_list}
)

# --- Precision-specific library ---
set(PRECISION_SOURCES
    {precision_list}
)

# Header include path
include_directories(${{CMAKE_CURRENT_SOURCE_DIR}}/src)

if(COMMON_SOURCES)
    add_library({common_lib} STATIC ${{COMMON_SOURCES}})
endif()

add_library({precision_lib} STATIC ${{PRECISION_SOURCES}})
if(TARGET {common_lib})
    target_link_libraries({precision_lib} PUBLIC {common_lib})
endif()
{c_mf_deps}
# --- Install rules ---
install(TARGETS {precision_lib} ARCHIVE DESTINATION lib)
if(TARGET {common_lib})
    install(TARGETS {common_lib} ARCHIVE DESTINATION lib)
endif()
"""
    else:
        # If multifloats, we need to link against the multifloats library
        # AND build the la_constants_mf / la_xisnan_mf helper modules
        # that the migrated source depends on for la_constants USE clauses.
        mf_link = ""
        mf_deps = ""
        if target_mode.module_name is not None:
            # Resolve absolute paths to external dependencies so the
            # generated CMakeLists.txt works from any output directory.
            # MF helpers (la_constants_mf.f90 / la_xisnan_mf.f90) live
            # under recipes/<lib>/mf_helpers/. Prefer the per-recipe
            # directory if present; otherwise fall back to the upstream
            # SRC dir which (historically) shipped the EP helpers.
            _root = project_root or Path.cwd()
            _per_recipe_mf = _root / 'recipes' / lib_name / 'mf_helpers'
            if _per_recipe_mf.is_dir():
                _helpers_default = str(_per_recipe_mf.resolve())
            else:
                _helpers_default = str((_root / 'external' / 'lapack-3.12.1' / 'SRC').resolve())
            # multifloats-mpi extras: Fortran-side MPI handle module
            # used by MUMPS (``USE multifloats_mpi_f``).
            _mf_mpi_dir = (_root / 'external' / 'multifloats-mpi').resolve()
            mf_link = f"""
# Fetch the multifloats library from GitHub (default) or use a local
# checkout via -DMULTIFLOATS_DIR=/path/to/multifloats. We add the
# multifloats *top-level* directory so its CMakeLists.txt runs — the
# src/CMakeLists.txt references ``CMAKE_SOURCE_DIR/include`` which is
# wrong outside a top-level build. Tests/benches are suppressed via
# cache variables set before the subdirectory add.
set(BUILD_TESTING OFF CACHE BOOL "Disable tests in fetched multifloats" FORCE)
set(MULTIFLOATS_BUILD_BENCH OFF CACHE BOOL "Disable benches in fetched multifloats" FORCE)
if(DEFINED MULTIFLOATS_DIR)
    message(STATUS "Using local multifloats: ${{MULTIFLOATS_DIR}}")
    add_subdirectory(${{MULTIFLOATS_DIR}}
        ${{CMAKE_CURRENT_BINARY_DIR}}/_mf EXCLUDE_FROM_ALL)
else()
    include(FetchContent)
    set(MULTIFLOATS_GIT_REPO "https://github.com/kyungminlee/multifloats.git"
        CACHE STRING "Git URL for the multifloats library")
    # multifloats v0.6.0 fixed the ${{CMAKE_SOURCE_DIR}} include-path
    # leak (upstream issue #23). Earlier tags fail at configure when
    # add_subdirectory'd; don't drop below this floor.
    set(MULTIFLOATS_GIT_TAG "v0.6.0"
        CACHE STRING "Git tag/branch/commit for multifloats (>= v0.6.0)")
    message(STATUS "Fetching multifloats from ${{MULTIFLOATS_GIT_REPO}} (${{MULTIFLOATS_GIT_TAG}})")
    FetchContent_Declare(multifloats_fetch
        GIT_REPOSITORY ${{MULTIFLOATS_GIT_REPO}}
        GIT_TAG        ${{MULTIFLOATS_GIT_TAG}}
    )
    FetchContent_Populate(multifloats_fetch)
    add_subdirectory(
        ${{multifloats_fetch_SOURCE_DIR}}
        ${{CMAKE_CURRENT_BINARY_DIR}}/_mf EXCLUDE_FROM_ALL)
endif()

# Build the la_constants_mf and la_xisnan_mf helper modules. These
# re-export multifloats's MF_* constants under the M/W-prefixed names
# that the migrated LAPACK source uses via its rewritten
# ``USE LA_CONSTANTS_MF`` clause.
set(MF_HELPERS_DIR "{_helpers_default}"
    CACHE PATH "Directory containing la_constants_mf.f90 / la_xisnan_mf.f90")
if(EXISTS "${{MF_HELPERS_DIR}}/la_constants_mf.f90")
    add_library(la_constants_mf STATIC
        "${{MF_HELPERS_DIR}}/la_constants_mf.f90")
    set_target_properties(la_constants_mf PROPERTIES
        Fortran_MODULE_DIRECTORY ${{CMAKE_CURRENT_BINARY_DIR}}/mod)
    target_include_directories(la_constants_mf PUBLIC
        $<BUILD_INTERFACE:${{CMAKE_CURRENT_BINARY_DIR}}/mod>)
    # la_constants_mf.f90 does ``use multifloats, only: real64x2`` —
    # the Fortran module lives in ``multifloatsf`` (the Fortran half),
    # not the C++ ``multifloats`` target.
    if(TARGET multifloatsf)
        target_link_libraries(la_constants_mf PUBLIC
            $<BUILD_INTERFACE:multifloatsf>)
    endif()
    if(TARGET multifloats)
        target_link_libraries(la_constants_mf PUBLIC multifloats)
    endif()
endif()
if(EXISTS "${{MF_HELPERS_DIR}}/la_xisnan_mf.f90")
    add_library(la_xisnan_mf STATIC
        "${{MF_HELPERS_DIR}}/la_xisnan_mf.f90")
    set_target_properties(la_xisnan_mf PROPERTIES
        Fortran_MODULE_DIRECTORY ${{CMAKE_CURRENT_BINARY_DIR}}/mod)
    target_include_directories(la_xisnan_mf PUBLIC
        $<BUILD_INTERFACE:${{CMAKE_CURRENT_BINARY_DIR}}/mod>)
    if(TARGET multifloatsf)
        target_link_libraries(la_xisnan_mf PUBLIC
            $<BUILD_INTERFACE:multifloatsf>)
    endif()
    if(TARGET multifloats)
        target_link_libraries(la_xisnan_mf PUBLIC multifloats)
    endif()
endif()

# multifloats_mpi_f.f90: Fortran module exposing the C-side MPI
# datatype handles (MPI_FLOAT64X2 / MPI_DD_SUM / ...) via bind(c).
# MUMPS's migrated source `USE multifloats_mpi_f` requires the .mod;
# other libraries route MPI through C and don't need this target.
set(MF_MPI_DIR "{_mf_mpi_dir}"
    CACHE PATH "Directory containing multifloats_mpi_f.f90")
if(EXISTS "${{MF_MPI_DIR}}/multifloats_mpi_f.f90")
    add_library(multifloats_mpi_f STATIC
        "${{MF_MPI_DIR}}/multifloats_mpi_f.f90")
    set_target_properties(multifloats_mpi_f PROPERTIES
        Fortran_MODULE_DIRECTORY ${{CMAKE_CURRENT_BINARY_DIR}}/mod)
    target_include_directories(multifloats_mpi_f PUBLIC
        $<BUILD_INTERFACE:${{CMAKE_CURRENT_BINARY_DIR}}/mod>)
    if(MPI_Fortran_FOUND)
        target_link_libraries(multifloats_mpi_f PUBLIC MPI::MPI_Fortran)
    endif()
endif()
"""
            mf_deps = f"""
if(TARGET multifloats)
    target_link_libraries({precision_lib} PUBLIC multifloats)
endif()
# multifloatsf is the Fortran half of multifloats — provides the
# ``multifloats.mod`` module that the migrated source's ``use
# multifloats`` clauses resolve against. Wrap in $<BUILD_INTERFACE:>
# because multifloats owns its own install/export set; we just need
# the .mod path during this build.
if(TARGET multifloatsf)
    target_link_libraries({precision_lib} PUBLIC
        $<BUILD_INTERFACE:multifloatsf>)
endif()
if(TARGET la_constants_mf)
    target_link_libraries({precision_lib} PUBLIC la_constants_mf)
endif()
if(TARGET la_xisnan_mf)
    target_link_libraries({precision_lib} PUBLIC la_xisnan_mf)
endif()
if(TARGET multifloats_mpi_f)
    target_link_libraries({precision_lib} PUBLIC multifloats_mpi_f)
endif()
"""

        cmake = f"""\
cmake_minimum_required(VERSION 3.20)
project({precision_lib} Fortran C CXX)

# --- Compiler flags ---
set(CMAKE_Fortran_FLAGS "${{CMAKE_Fortran_FLAGS}} -w")
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${{CMAKE_Fortran_FLAGS}} -std=legacy")
endif()

# Enable Fortran preprocessing for .F90 files
set(CMAKE_Fortran_PREPROCESS ON)

# --- MPI: default to vendored Intel MPI headers ---
# See note in the C template: headers come from ``external/impi-headers/``
# unconditionally; the runtime comes from whichever MPI the user links
# against at final link time. MUMPS uses ``INCLUDE 'mpif.h'`` in 231
# source files and never ``USE mpi``, so F77 headers are enough.
if(NOT DEFINED IMPI_HEADERS)
    set(IMPI_HEADERS "{_impi_default}"
        CACHE PATH "Path to vendored Intel MPI headers")
endif()
include_directories(${{IMPI_HEADERS}})
find_package(MPI COMPONENTS Fortran QUIET)

# Detect extended-precision (KIND=10 / KIND=16) support.
# Shared probe sits in cmake/DetectExtendedPrecision.cmake; copied
# next to this CMakeLists.txt at generation time so the staging
# tree stays self-contained.
include(${{CMAKE_CURRENT_SOURCE_DIR}}/DetectExtendedPrecision.cmake)
{mf_link}
# --- Standard-precision sibling archive ---
# Built from upstream Fortran sources alongside the migrated archive.
# Carries the original S/D/C/Z entry points and the precision-
# independent helpers (LSAME, XERBLA, LA_XISNAN module, ...) that
# the migrated archive's bodies reference but don't ship themselves.
# The migrated archive PUBLIC-links this so downstreams resolve both
# symbol families through one link line. Modules also flow to the
# migrated build via the shared module directory (e.g. la_xisnan.mod
# from std-precision is what la_xisnan_ep.F90 `use`s).
set(REF_SOURCES
    {ref_list}
)

if(REF_SOURCES)
    add_library({lib_name} STATIC ${{REF_SOURCES}})
    set_target_properties({lib_name} PROPERTIES
        Fortran_MODULE_DIRECTORY ${{CMAKE_CURRENT_BINARY_DIR}}/mod)
    target_include_directories({lib_name} PUBLIC
        $<BUILD_INTERFACE:${{CMAKE_CURRENT_BINARY_DIR}}/mod>)
    # The migrator can lengthen tokens (e.g. MPI_DOUBLE_COMPLEX →
    # MPI_C_LONG_DOUBLE_COMPLEX), pushing fixed-form .F lines past
    # column 72. Disable the line-length limit on the standard
    # archive too — same flag as the precision archive — so the
    # build stays consistent across the pair.
    if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
        target_compile_options({lib_name} PRIVATE
            $<$<COMPILE_LANGUAGE:Fortran>:-ffixed-line-length-none>)
    endif()
endif()

# --- Common (type-independent) library ---
set(COMMON_SOURCES
    {common_list}
)

# --- Precision-specific library ---
set(PRECISION_SOURCES
    {precision_list}
)

if(COMMON_SOURCES)
    add_library({common_lib} STATIC ${{COMMON_SOURCES}})
    set_target_properties({common_lib} PROPERTIES
        Fortran_MODULE_DIRECTORY ${{CMAKE_CURRENT_BINARY_DIR}}/mod)
    target_include_directories({common_lib} PUBLIC
        $<BUILD_INTERFACE:${{CMAKE_CURRENT_BINARY_DIR}}/mod>)
endif()

add_library({precision_lib} STATIC ${{PRECISION_SOURCES}})
set_target_properties({precision_lib} PROPERTIES
    Fortran_MODULE_DIRECTORY ${{CMAKE_CURRENT_BINARY_DIR}}/mod)
target_include_directories({precision_lib} PUBLIC
    $<BUILD_INTERFACE:${{CMAKE_CURRENT_BINARY_DIR}}/mod>
    $<BUILD_INTERFACE:${{CMAKE_CURRENT_SOURCE_DIR}}>)
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    target_compile_options({precision_lib} PRIVATE
        $<$<COMPILE_LANGUAGE:Fortran>:-ffixed-line-length-none>)
endif()
if(TARGET {common_lib})
    target_link_libraries({precision_lib} PUBLIC {common_lib})
endif()
if(TARGET {lib_name})
    target_link_libraries({precision_lib} PUBLIC {lib_name})
endif()
if(MPI_Fortran_FOUND)
    target_link_libraries({precision_lib} PUBLIC MPI::MPI_Fortran)
    if(TARGET {common_lib})
        target_link_libraries({common_lib} PUBLIC MPI::MPI_Fortran)
    endif()
    if(TARGET {lib_name})
        target_link_libraries({lib_name} PUBLIC MPI::MPI_Fortran)
    endif()
endif()
{mf_deps}

# --- Install rules ---
install(TARGETS {precision_lib} ARCHIVE DESTINATION lib)
if(TARGET {common_lib})
    install(TARGETS {common_lib} ARCHIVE DESTINATION lib)
endif()
if(TARGET {lib_name})
    install(TARGETS {lib_name} ARCHIVE DESTINATION lib)
endif()
"""
    (output_dir / 'CMakeLists.txt').write_text(cmake)

    # Ship the shared extended-precision probe alongside the generated
    # CMakeLists.txt. The Fortran template ``include(...)``s it; the C
    # template doesn't need it but the file is cheap to copy and keeps
    # the staging tree self-contained.
    if language != 'c':
        _root = project_root or Path.cwd()
        probe = _root / 'cmake' / 'DetectExtendedPrecision.cmake'
        if probe.exists():
            shutil.copy2(probe, output_dir / probe.name)


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
                           config.extensions, config.library_path,
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
    # copy_files are routed to the precision library unconditionally:
    # their stems may collide with module names the symbol scanner has
    # now classified as ``independent`` (post-MODULE-scanner), but the
    # bodies often USE precision-specific modules that live in the
    # precision lib — a common-lib copy would create a forbidden
    # common -> precision dependency. Keeping copy_files in precision
    # preserves the one-way precision -> common link direction.
    common_files, precision_files = [], []
    for f in files:
        rel = f.relative_to(output_dir)
        stem = f.stem.upper()
        if stem in config.copy_files:
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

    print(f'\nConfiguring...')
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
    real_pfx = pmap['R'].lower()
    precision_lib_name = f'lib{real_pfx}{lib_name}.a'
    common_lib_name = f'lib{lib_name}_common.a'
    ref_lib_name = f'lib{lib_name}.a'

    print(f'\nBuild succeeded:')
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


# Topologically sorted library build order for the unified CMake project.
# Each entry is (library_name, recipe_filename).
LIBRARY_ORDER = [
    ('blas',        'blas.yaml'),
    ('xblas',       'xblas.yaml'),
    ('blacs',       'blacs.yaml'),
    ('lapack',      'lapack.yaml'),
    ('ptzblas',     'ptzblas.yaml'),
    # NUMROC / ICEIL / ILCM. Hoisted out of scalapack so that pbblas
    # downstreams (which don't link libqscalapack) can still reach them.
    ('scalapack_tools', 'scalapack_tools.yaml'),
    ('pbblas',      'pbblas.yaml'),
    ('pblas',       'pblas.yaml'),
    ('scalapack',   'scalapack.yaml'),
    ('scalapack_c', 'scalapack_c.yaml'),
    ('mumps',       'mumps.yaml'),
]


BASELINE_TARGETS = ('kind4', 'kind8')


def cmd_stage(args):
    """Migrate all libraries into a structured staging directory.

    Produces a self-contained directory that can be built with:
        cmake -S <staging> -B <staging>/build && cmake --build <staging>/build -j

    For ``--target kind4`` / ``--target kind8`` the migration step is
    skipped: those are un-migrated single/double precision baselines used
    for diff'ing against the quad reference, so all that's needed is a
    staging tree with upstream sources, tests/, and a target_config.cmake
    that points the test framework at the standard archive (LIB_PREFIX="").
    """
    target_str = getattr(args, 'target', None) or 'kind16'
    if target_str in BASELINE_TARGETS:
        return _stage_baseline(args, target_str)

    staging_dir = args.staging_dir.resolve()
    target_mode = _get_target_mode(args)
    parser, parser_cmd = _parser_args(args)
    proj_root = (args.project_root or Path(__file__).resolve().parent.parent.parent)
    recipes_dir = proj_root / 'recipes'

    # Determine which libraries to stage
    if args.libraries:
        lib_set = set(args.libraries)
        libraries = [(n, r) for n, r in LIBRARY_ORDER if n in lib_set]
    else:
        libraries = list(LIBRARY_ORDER)

    staged = []
    for lib_name, recipe_file in libraries:
        recipe_path = recipes_dir / recipe_file
        if not recipe_path.exists():
            print(f'Warning: recipe {recipe_path} not found, skipping {lib_name}')
            continue

        lib_dir = staging_dir / lib_name
        src_dir = lib_dir / 'src'
        src_dir.mkdir(parents=True, exist_ok=True)

        print(f'\n{"=" * 60}')
        print(f'  Migrating: {lib_name}')
        print(f'{"=" * 60}')

        # Run migration
        run_migration(
            recipe_path=recipe_path,
            output_dir=src_dir,
            target_mode=target_mode,
            dry_run=False,
            project_root=proj_root,
            parser=parser,
            parser_cmd=parser_cmd,
        )

        # Classify files into common vs precision-specific
        config = prepare_recipe(recipe_path, proj_root)
        symbols = scan_symbols(config.source_dir, config.language,
                               config.extensions, config.library_path,
                               extra_c_return_types=tuple(config.c_return_types))
        classification = classify_symbols(symbols)
        independent = classification.independent

        # Pick up precision-independent Fortran helpers staged via
        # ``copy_files`` (e.g. PBLAS/SRC/pilaenv.f) when the recipe is C —
        # CMake's ``add_library(… STATIC …)`` handles mixed C + Fortran
        # sources because both languages are enabled at the top
        # project(). Copy-files are precision-independent by contract,
        # so they land in COMMON_SOURCES below.
        files = _collect_source_files(src_dir, config.language)

        # MF helper modules are built as separate targets; exclude from manifests.
        _mf_helpers = {'la_constants_mf', 'la_xisnan_mf'}
        # la_constants_ep provides extended/quad precision constants —
        # not needed by multifloats target, and may fail on compilers
        # that don't support REAL(KIND=16).
        if target_mode.module_name is not None:
            _mf_helpers.update({'la_constants_ep', 'la_xisnan_ep'})

        common_files, precision_files = [], []
        for f in files:
            if f.stem in _mf_helpers:
                continue
            rel = f'src/{f.name}'
            # ``copy_files`` entries are precision-independent by
            # contract (the file is staged verbatim, no prefix rename).
            # The symbol scanner may never have visited them — e.g. a
            # Fortran ``copy_files`` entry in a C recipe — so they
            # won't appear in ``independent``. Treat them as common
            # explicitly.
            if (f.stem.upper() in independent
                    or f.stem.upper() in config.copy_files):
                common_files.append(rel)
            else:
                precision_files.append(rel)

        # Identify C sources that gate their entry-point signature on
        # the ``INTFACE == C_CALL`` macro (upstream BLACS pattern). Each
        # such file exposes a Fortran-callable symbol (e.g.
        # ``blacs_gridinfo_``) in the default build and a C-callable
        # symbol (``Cblacs_gridinfo``) when compiled with
        # ``-DCallFromC``. The CMake helper compiles these sources
        # twice so the final static library ships both entry points.
        # Detection is a cheap regex scan of the staged source; any
        # library that does not use the pattern emits an empty list.
        dual_files = []
        if config.language == 'c':
            for f in files:
                if f.suffix.lower() != '.c':
                    continue
                try:
                    text = f.read_text(errors='replace')
                except OSError:
                    continue
                if _DUAL_ENTRY_C_RE.search(text):
                    dual_files.append(f'src/{f.name}')

        # Write manifest.cmake
        common_list = '\n    '.join(common_files) if common_files else ''
        precision_list = '\n    '.join(precision_files) if precision_files else ''
        dual_list = '\n    '.join(dual_files) if dual_files else ''
        manifest = f"""\
set({lib_name}_COMMON_SOURCES
    {common_list}
)

set({lib_name}_PRECISION_SOURCES
    {precision_list}
)

set({lib_name}_DUAL_INTERFACE_SOURCES
    {dual_list}
)

set({lib_name}_LANGUAGE {config.language})
"""
        (lib_dir / 'manifest.cmake').write_text(manifest)
        print(f'  Manifest: {len(common_files)} common, '
              f'{len(precision_files)} precision files')
        staged.append(lib_name)

    # Copy MF helper modules into staging so it's self-contained
    pmap = target_mode.prefix_map
    lib_prefix = pmap['R'].lower()
    lib_prefix_complex = pmap['C'].lower()
    needs_mf = target_mode.module_name is not None
    staged_list = ';'.join(staged)

    helpers_src = proj_root / 'recipes' / 'lapack' / 'mf_helpers'
    helpers_dst = staging_dir / '_helpers'
    if needs_mf:
        helpers_dst.mkdir(exist_ok=True)
        for name in ['la_constants_mf.f90', 'la_xisnan_mf.f90']:
            src = helpers_src / name
            if src.exists():
                shutil.copy2(src, helpers_dst / name)
        # Copy multifloats bridge files (C++ bridge header + MPI registration)
        mf_local = proj_root / 'external' / 'multifloats-mpi'
        bridge_h = mf_local / 'multifloats_bridge.h'
        mpi_cpp = mf_local / 'multifloats_mpi.cpp'
        if bridge_h.exists():
            shutil.copy2(bridge_h, helpers_dst / bridge_h.name)
            # Skip the C++ MPI bindings (mpicxx.h). Without this guard,
            # any migrated source compiled as C++ that transitively
            # pulls in <mpi.h> through the bridge gets thousands of
            # template declarations from mpicxx.h. Those templates
            # cannot live inside the ``extern "C" { … }`` wrap that
            # the c_migrator post-pass injects around .c bodies, so
            # link of scalapack_c (whose REDIST sources include
            # redist.h → multifloats_bridge.h → mpi.h) fails. Setting
            # MPICH_SKIP_MPICXX / OMPI_SKIP_MPICXX before the include
            # is the documented way to compile MPI clients without
            # the C++ bindings.
            staged_bridge = helpers_dst / bridge_h.name
            text = staged_bridge.read_text()
            text = text.replace(
                '#include <mpi.h>',
                '#define MPICH_SKIP_MPICXX 1\n'
                '#define OMPI_SKIP_MPICXX 1\n'
                '#include <mpi.h>',
                1,
            )
            staged_bridge.write_text(text)
        if mpi_cpp.exists():
            shutil.copy2(mpi_cpp, helpers_dst / mpi_cpp.name)
        # multifloats_mpi_f.f90: Fortran module exposing the C-side
        # MPI_FLOAT64X2 / MPI_DD_SUM / etc. handles via bind(c). MUMPS
        # and any other library that calls MPI from Fortran directly
        # USEs this module so the rewritten ``MPI_FLOAT64X2`` token is
        # a known integer at compile time. (BLACS/PBLAS go through C
        # and use the extern "C" handles from the bridge header
        # instead.)
        mpi_f90 = mf_local / 'multifloats_mpi_f.f90'
        if mpi_f90.exists():
            shutil.copy2(mpi_f90, helpers_dst / mpi_f90.name)

    target_config = f"""\
# Generated by: python -m migrator stage --target {target_mode.name}
set(TARGET_NAME "{target_mode.name}")
set(LIB_PREFIX "{lib_prefix}")
set(LIB_PREFIX_COMPLEX "{lib_prefix_complex}")
set(NEEDS_MULTIFLOATS {'TRUE' if needs_mf else 'FALSE'})
set(C_AS_CXX {'TRUE' if needs_mf else 'FALSE'})
set(MF_HELPERS_DIR "${{CMAKE_CURRENT_SOURCE_DIR}}/_helpers")
set(STAGED_LIBRARIES {staged_list})
"""
    (staging_dir / 'target_config.cmake').write_text(target_config)

    # Copy CMake files to staging directory. ``CMakePresets.json`` rides
    # along so users can `cmake --preset=linux-impi` from the staged
    # tree without having to re-discover Intel MPI's wrapper paths.
    cmake_dir = proj_root / 'cmake'
    for cmake_file in ['CMakeLists.txt', 'FortranCompiler.cmake',
                       'DetectExtendedPrecision.cmake',
                       'CMakePresets.json',
                       'mpiseq_qx_stubs.f',
                       'mpiseq_mw_stubs.f90',
                       'mpiseq_c_stubs.c']:
        src = cmake_dir / cmake_file
        if src.exists():
            shutil.copy2(src, staging_dir / cmake_file)
        else:
            print(f'Warning: {src} not found')

    # Plant the refblas_quad / reflapack_quad symbol-rename helper
    # alongside the other build-time scripts so tests/blas/refblas and
    # tests/lapack/reflapack can locate it via find_file in the staging
    # tree (see those CMakeLists for the search-path list).
    scripts_src = proj_root / 'scripts' / 'refquad_rename_archive.sh'
    if scripts_src.exists():
        scripts_dst = staging_dir / 'scripts'
        scripts_dst.mkdir(exist_ok=True)
        shutil.copy2(scripts_src, scripts_dst / scripts_src.name)

    # Copy tests/ subtree so the unified CMakeLists.txt can pick it up
    # via add_subdirectory(tests) when BUILD_TESTING=ON.
    tests_src = proj_root / 'tests'
    if tests_src.is_dir():
        tests_dst = staging_dir / 'tests'
        if tests_dst.exists():
            shutil.rmtree(tests_dst)
        shutil.copytree(tests_src, tests_dst)

    # Copy parallel_blas/ overlay sources (opt-in via PARALLEL_BLAS=ON
    # at CMake time). The unified CMakeLists picks it up via
    # add_subdirectory(parallel_blas) when the directory is present.
    parallel_blas_src = proj_root / 'src' / 'parallel_blas'
    if parallel_blas_src.is_dir():
        parallel_blas_dst = staging_dir / 'parallel_blas'
        if parallel_blas_dst.exists():
            shutil.rmtree(parallel_blas_dst)
        shutil.copytree(parallel_blas_src, parallel_blas_dst)

    # Copy epopenblas/ overlay sources — separate OpenBLAS-D/Z-port-at-kind10
    # overlay archive (NOT wired into the public composite). The unified
    # CMakeLists picks it up via add_subdirectory(epopenblas) when the
    # directory is present.
    epopenblas_src = proj_root / 'src' / 'epopenblas'
    if epopenblas_src.is_dir():
        epopenblas_dst = staging_dir / 'epopenblas'
        if epopenblas_dst.exists():
            shutil.rmtree(epopenblas_dst)
        shutil.copytree(epopenblas_src, epopenblas_dst)

    # Copy vendored Netlib BLAS source for the differential precision
    # tests' refblas_quad reference library (compiled with gfortran's
    # -freal-8-real-16 to promote KIND=8 entities to KIND=16 in-place).
    # Tests fall back to system -lblas if this directory is absent.
    netlib_blas_src = proj_root / 'external' / 'lapack-3.12.1' / 'BLAS' / 'SRC'
    if netlib_blas_src.is_dir():
        refblas_dst = staging_dir / '_refblas_src'
        if refblas_dst.exists():
            shutil.rmtree(refblas_dst)
        shutil.copytree(netlib_blas_src, refblas_dst)

    # Same recipe for LAPACK: vendored Netlib SRC/ promoted to quad
    # precision gives tests/lapack/reflapack/ a KIND=16 reference to
    # compare the migrated qlapack/elapack/ddlapack against. The
    # INSTALL/ directory provides dlamch.f / droundup_lwork.f, which
    # LAPACK SRC routines call but which aren't in SRC/ itself — copy
    # them into _reflapack_src/ alongside the SRC contents so a single
    # glob compiles the full reference.
    netlib_lapack_src = proj_root / 'external' / 'lapack-3.12.1' / 'SRC'
    if netlib_lapack_src.is_dir():
        reflapack_dst = staging_dir / '_reflapack_src'
        if reflapack_dst.exists():
            shutil.rmtree(reflapack_dst)
        shutil.copytree(netlib_lapack_src, reflapack_dst)
        install_src = proj_root / 'external' / 'lapack-3.12.1' / 'INSTALL'
        for fname in ('dlamch.f', 'droundup_lwork.f'):
            src = install_src / fname
            if src.is_file():
                shutil.copy2(src, reflapack_dst / fname)

    # Stage standard-precision source directories for the std archives
    # built alongside each migrated extension. The CMakeLists.txt
    # invokes add_standard_fortran_library / add_standard_c_library
    # against these directories. Sibling to _refblas_src/_reflapack_src
    # but used for production link deps, not just tests.
    # _pblas_src/ includes PTOOLS/ as a child subdirectory (matching
    # the upstream layout) so PTOOLS sources' ``#include "../pblas.h"``
    # resolves without an explicit include-path remap. Same shape for
    # _scalapack_src/ which contains REDIST/SRC and shares the
    # ``../some_header.h`` convention. PBLAS's internal subdirectories
    # (PBBLAS, PTZBLAS) are NOT included under _pblas_src — those are
    # owned by the separate ptzblas / pbblas std archives.
    # _scalapack_src/ is the upstream SRC/ tree; scalapack_c REDIST
    # routines live alongside SRC under REDIST/SRC, but we stage them
    # together inside _scalapack_src/REDIST/ so REDIST sources'
    # ``#include "../redist.h"`` (or the matching SRC headers)
    # resolve relative to _scalapack_src/.
    _std_dirs = [
        ('_blacs_src',     'scalapack-2.2.3/BLACS/SRC'),
        ('_pblas_src',     'scalapack-2.2.3/PBLAS/SRC'),
        ('_ptzblas_src',   'scalapack-2.2.3/PBLAS/SRC/PTZBLAS'),
        ('_pbblas_src',    'scalapack-2.2.3/PBLAS/SRC/PBBLAS'),
        ('_scalapack_src', 'scalapack-2.2.3/SRC'),
        ('_scalapack_tools_src', 'scalapack-2.2.3/TOOLS'),
        ('_scalapack_redist_src', 'scalapack-2.2.3/REDIST/SRC'),
        # MUMPS sequential MPI stub. Lets cmake build a single-process
        # ``libmpiseq`` archive alongside the migrated qmumps; tests can
        # link it instead of MPI::MPI_Fortran for plain (no mpiexec)
        # executables. Stubs print a "should not be called" error if a
        # collective/comm primitive that requires multi-rank coordination
        # is invoked, so libseq is NPROCS=1-only by construction.
        ('_mpiseq_src',    'MUMPS_5.8.2/libseq'),
        # MUMPS upstream src/ + include/. The recipe (which is fortran-
        # only) skips every *MUMPS_C / MUMPS_C_TYPES header and every
        # *.c file, so the migrated qmumps archive ships without a C
        # interface. tests/mumps's C-bridge build re-uses upstream
        # mumps_c.c (compiled twice with quad-precision type overrides
        # supplied from tests/mumps/c/include/, see B2 in
        # tests/mumps/TODO.md), plus mumps_common.c, mumps_addr.c, and
        # the IO/save/thread/metis/pord/scotch helpers, all of which are
        # type-agnostic and compile verbatim. Staging the whole src/
        # tree (including the .F siblings we don't need here) is
        # cheaper than per-file plumbing and matches the convention.
        ('_mumps_upstream_src',     'MUMPS_5.8.2/src'),
        ('_mumps_upstream_include', 'MUMPS_5.8.2/include'),
    ]
    for dst_name, rel_src in _std_dirs:
        src = proj_root / 'external' / rel_src
        if not src.is_dir():
            continue
        dst = staging_dir / dst_name
        if dst.exists():
            shutil.rmtree(dst)
        shutil.copytree(src, dst)

    # libseq's mpi.f bundles BLACS/ScaLAPACK forwarders (which collide
    # with the real migrated archives) and only knows the standard MPI
    # datatypes in MUMPS_COPY (no MPI_REAL16 / MPI_COMPLEX32 needed by
    # qmumps reductions). Patch the staged copy to fix both. Upstream
    # external/ stays read-only.
    mpiseq_dst = staging_dir / '_mpiseq_src' / 'mpi.f'
    if mpiseq_dst.is_file():
        _patch_libseq_mpi_f(mpiseq_dst)

    print(f'\n{"=" * 60}')
    print(f'  Staging complete: {len(staged)} libraries')
    print(f'{"=" * 60}')
    print(f'  Target:  {target_mode.name} (prefix: {lib_prefix})')
    print(f'  Output:  {staging_dir}')
    print(f'\nTo build:')
    print(f'  cmake -S {staging_dir} -B {staging_dir}/build -DCMAKE_BUILD_TYPE=Release')
    print(f'  cmake --build {staging_dir}/build -j')


def _stage_baseline(args, target_name: str):
    """Stage an unmigrated baseline tree for kind4 / kind8.

    No per-library migration is run — kind4 / kind8 are the upstream
    S/D/C/Z entry points themselves, served by the ``add_standard_*``
    archives the unified CMake build always assembles. We just need the
    upstream source trees, the tests/ subtree, the cmake/ glue, and a
    target_config.cmake that signals "no migrated archives" via an
    empty LIB_PREFIX so the test framework links the std archive
    directly.
    """
    staging_dir = args.staging_dir.resolve()
    proj_root = (args.project_root or
                 Path(__file__).resolve().parent.parent.parent)
    staging_dir.mkdir(parents=True, exist_ok=True)

    print(f'\n{"=" * 60}')
    print(f'  Baseline staging: {target_name} (no migration)')
    print(f'{"=" * 60}')

    # target_config.cmake: empty LIB_PREFIX, no multifloats, no migrated
    # libs in STAGED_LIBRARIES. The parent CMakeLists.txt's
    # add_migrated_* helpers are no-ops without per-lib manifest.cmake
    # files (which we don't write here), so the build resolves to just
    # the standard archives.
    target_config = (
        f'# Generated by: python -m migrator stage --target {target_name}\n'
        f'# Baseline (un-migrated) target — see targets/{target_name}.yaml.\n'
        f'set(TARGET_NAME "{target_name}")\n'
        'set(LIB_PREFIX "")\n'
        'set(LIB_PREFIX_COMPLEX "")\n'
        'set(NEEDS_MULTIFLOATS FALSE)\n'
        'set(C_AS_CXX FALSE)\n'
        'set(MF_HELPERS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/_helpers")\n'
        'set(STAGED_LIBRARIES )\n'
    )
    (staging_dir / 'target_config.cmake').write_text(target_config)

    # CMake glue (top-level CMakeLists, FortranCompiler module, presets,
    # mpiseq stubs). Same files cmd_stage copies; baseline reuses them.
    cmake_dir = proj_root / 'cmake'
    for cmake_file in ['CMakeLists.txt', 'FortranCompiler.cmake',
                       'DetectExtendedPrecision.cmake',
                       'CMakePresets.json',
                       'mpiseq_qx_stubs.f',
                       'mpiseq_mw_stubs.f90',
                       'mpiseq_c_stubs.c']:
        src = cmake_dir / cmake_file
        if src.exists():
            shutil.copy2(src, staging_dir / cmake_file)

    scripts_src = proj_root / 'scripts' / 'refquad_rename_archive.sh'
    if scripts_src.exists():
        scripts_dst = staging_dir / 'scripts'
        scripts_dst.mkdir(exist_ok=True)
        shutil.copy2(scripts_src, scripts_dst / scripts_src.name)

    # tests/ subtree.
    tests_src = proj_root / 'tests'
    if tests_src.is_dir():
        tests_dst = staging_dir / 'tests'
        if tests_dst.exists():
            shutil.rmtree(tests_dst)
        shutil.copytree(tests_src, tests_dst)

    # Upstream BLAS / LAPACK / ScaLAPACK / PBLAS / BLACS / MUMPS sources.
    # For libraries with a recipe we stage from build/staged-sources/<lib>/
    # so the baseline column links against the patched archives (closing
    # the gap the kind4/kind8 column previously had — patched in migrated
    # archive but broken in baseline). Libraries without a recipe
    # (TOOLS / REDIST / libseq / MUMPS include/) come straight from
    # external/.
    recipes_dir = proj_root / 'recipes'

    def _staged_or_external(rel_src: str, recipe_name: str | None) -> Path:
        if recipe_name:
            recipe_path = recipes_dir / f'{recipe_name}.yaml'
            if recipe_path.exists():
                return run_prepare(recipe_path, project_root=proj_root)
        return proj_root / 'external' / rel_src

    def _stage_dst(dst_name: str, src: Path) -> None:
        if not src.is_dir():
            return
        dst = staging_dir / dst_name
        if dst.exists():
            shutil.rmtree(dst)
        shutil.copytree(src, dst, ignore=shutil.ignore_patterns('.prepared.stamp'))

    _stage_dst('_refblas_src',
               _staged_or_external('lapack-3.12.1/BLAS/SRC', 'blas'))

    reflapack_src = _staged_or_external('lapack-3.12.1/SRC', 'lapack')
    _stage_dst('_reflapack_src', reflapack_src)
    if reflapack_src.is_dir():
        # Pull dlamch / slamch / droundup_lwork / sroundup_lwork from
        # external/INSTALL (these aren't in any recipe's source_dir, so
        # they don't get staged; baseline needs them alongside SRC).
        install_src = proj_root / 'external' / 'lapack-3.12.1' / 'INSTALL'
        reflapack_dst = staging_dir / '_reflapack_src'
        for fname in ('dlamch.f', 'droundup_lwork.f', 'slamch.f',
                      'sroundup_lwork.f'):
            src = install_src / fname
            if src.is_file():
                shutil.copy2(src, reflapack_dst / fname)

    _std_dirs: list[tuple[str, str, str | None]] = [
        ('_blacs_src',            'scalapack-2.2.3/BLACS/SRC',         'blacs'),
        ('_pblas_src',            'scalapack-2.2.3/PBLAS/SRC',         'pblas'),
        ('_ptzblas_src',          'scalapack-2.2.3/PBLAS/SRC/PTZBLAS', 'ptzblas'),
        ('_pbblas_src',           'scalapack-2.2.3/PBLAS/SRC/PBBLAS',  'pbblas'),
        ('_scalapack_src',        'scalapack-2.2.3/SRC',               'scalapack'),
        ('_scalapack_tools_src',  'scalapack-2.2.3/TOOLS',             None),
        ('_scalapack_redist_src', 'scalapack-2.2.3/REDIST/SRC',        None),
        ('_mpiseq_src',           'MUMPS_5.8.2/libseq',                None),
        ('_mumps_upstream_src',   'MUMPS_5.8.2/src',                   'mumps'),
        ('_mumps_upstream_include', 'MUMPS_5.8.2/include',             None),
    ]
    for dst_name, rel_src, recipe_name in _std_dirs:
        _stage_dst(dst_name, _staged_or_external(rel_src, recipe_name))

    print(f'  Output:  {staging_dir}')
    print(f'\nTo build:')
    print(f'  cmake -S {staging_dir} -B {staging_dir}/build -DCMAKE_BUILD_TYPE=Release')
    print(f'  cmake --build {staging_dir}/build -j')
    print(f'  ctest --test-dir {staging_dir}/build')


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

def main():
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

    args = parser.parse_args()
    rc = args.func(args)
    if isinstance(rc, int) and rc != 0:
        sys.exit(rc)


if __name__ == '__main__':
    main()
