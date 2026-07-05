"""``stage`` subcommand: materialize a self-contained migrated build tree.

Enumerates and migrates each library's sources into a staging directory,
patches the LibSeq MPI shim, classifies symbols, and lays down the CMake
manifest a downstream build consumes. Extracted verbatim from ``__main__``.
"""
import re
import shutil
import sys
from pathlib import Path

from .cli_common import _get_target_mode, _parser_args
from .libseq_patch import _patch_libseq_mpi_f
from .pipeline import run_migration
from .prefix_classifier import classify_symbols
from .prepare import prepare_recipe, run_prepare
from .symbol_scanner import scan_symbols


# BLACS-style dual-entry-point detector used by ``cmd_stage`` to
# identify C sources that switch their public symbol via the
# ``INTFACE == C_CALL`` / ``CallFromC`` macros. Hoisted to module scope
# so it's compiled once instead of per-library in the cmd_stage loop.
_DUAL_ENTRY_C_RE = re.compile(
    r'#\s*if\s*\(?\s*INTFACE\s*==\s*C_CALL\b'
    r'|#\s*ifdef\s+CallFromC\b'
    r'|#\s*if\s+defined\s*\(\s*CallFromC\s*\)',
)


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


# Topologically sorted library build order for the unified CMake project.
# Each entry is (library_name, recipe_filename).
LIBRARY_ORDER = [
    ('blas',        'blas.yaml'),
    ('xblas',       'xblas.yaml'),
    ('blacs',       'blacs.yaml'),
    ('lapack',      'lapack.yaml'),
    ('ptzblas',     'ptzblas.yaml'),
    # Leaf symbol-provider for NUMROC / ICEIL / ILCM: kept so pbblas.yaml
    # keeps its cycle-free numroc symbol context. Its migrated output is
    # no longer built — the three helpers are folded into scalapack_common
    # (migrated by scalapack.yaml). Still staged (harmless, unused manifest).
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
        valid = {n for n, _ in LIBRARY_ORDER}
        unknown = lib_set - valid
        if unknown:
            sys.exit(
                f'error: unknown library name(s) in --libraries: '
                f'{sorted(unknown)}. Valid: {sorted(valid)}'
            )
        libraries = [(n, r) for n, r in LIBRARY_ORDER if n in lib_set]
    else:
        libraries = list(LIBRARY_ORDER)

    # The real / complex family letters (e/y, q/x, m/w). Used both to
    # build the per-library archive-filename prefix below and to stamp
    # LIB_PREFIX / LIB_PREFIX_COMPLEX into target_config.cmake.
    pmap = target_mode.prefix_map
    lib_prefix = pmap['R'].lower()
    lib_prefix_complex = pmap['C'].lower()

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

        # Files that the precision-rename map RENAMES are precision-
        # specific by construction — they were given a precision prefix
        # (e.g. ``disnan.f`` → ``qisnan.f``, family ``#RISNAN``, d→q).
        # Their renamed stem can COINCIDE with a name the classifier
        # marked ``independent``: the split ``la_xisnan_ey.F90`` /
        # ``la_xisnan_qx.F90`` modules define procedures ``EISNAN`` /
        # ``QISNAN`` / ``ELAISNAN`` / ``QLAISNAN``, which land in
        # ``independent`` because Q/E/X/Y are not S/D/C/Z.
        # Without this guard the migrated plain-F77 isnan externals
        # (``qisnan_``/``eisnan_`` — whose bodies differ per target) get
        # mis-filed into the shared ``lapack_common`` archive and then
        # collide across targets in a combined release tree (same
        # filename, different content). A rename target is never a
        # legitimate common member, so route it to PRECISION regardless
        # of the polluted independent set.
        rename_targets = {
            t.upper()
            for t in classification.build_rename_map(target_mode).values()
        }

        # Pick up precision-independent Fortran helpers staged via
        # ``copy_files`` (e.g. PBLAS/SRC/pilaenv.f) when the recipe is C —
        # CMake's ``add_library(… STATIC …)`` handles mixed C + Fortran
        # sources because both languages are enabled at the top
        # project(). Copy-files are precision-independent by contract,
        # so they land in COMMON_SOURCES below.
        files = _collect_source_files(src_dir, config.language)

        # Per-target LA_CONSTANTS / LA_XISNAN precision helpers. Each
        # extended target owns exactly one prefix-pair module —
        #   kind10 → *_ey   kind16 → *_qx   multifloats → *_mw
        # (suffix taken from the target's ``la_constants_suffix``). Only
        # the target's own pair belongs in this staging tree; the other
        # targets' pairs, which sit in the shared LAPACK SRC dir, are
        # excluded. multifloats' *_mw pair is not in SRC at all — it
        # ships as separate helper libraries from
        # recipes/lapack/mf_helpers/ (staged below) — so this exclusion
        # is also what keeps the *_ey/*_qx SRC files out of a
        # multifloats build (formerly the ``module_name is not None``
        # special-case for the shared *_ep pair).
        _la_suffixes = ('_ey', '_qx', '_mw')
        _own_suffix = target_mode.la_constants_suffix.lower()
        _la_own = {f'la_constants{_own_suffix}', f'la_xisnan{_own_suffix}'}
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
            rel = f'src/{f.name}'
            stem = f.stem.upper()
            # ``force_common`` is the explicit recipe override that pins a
            # stem to the family-independent ``_common`` archive no matter
            # what the symbol scanner decided. It takes priority over every
            # other route (including the PRECISION defaults below) so a
            # recipe author can rescue files the scanner mis-assigned to a
            # non-target family — e.g. the integer BLACS entry points and
            # the type-agnostic driver, which are family-independent (one
            # copy serves e/y, q/x, m/w) and must not leak into the
            # prefixed ``eyblacs`` archive. See ``recipes/blacs.yaml``.
            # Match with OR without the trailing Fortran-underscore so a
            # recipe can list the logical name (``igamx2d``) exactly as
            # ``skip_files`` does (the C migrator strips it — see
            # c_migrator.py), keeping the two levers' conventions aligned.
            stem_nodeco = stem[:-1] if stem.endswith('_') else stem
            if stem in config.force_common or stem_nodeco in config.force_common:
                common_files.append(rel)
            # The target's own LA_CONSTANTS/LA_XISNAN pair is single-
            # precision (only this target's E/Y or Q/X block), so it is
            # precision-specific: route it to PRECISION so it lands in
            # the prefixed archive (libelapack / libqlapack), never the
            # shared lapack_common. A precision-rename target (e.g. the
            # migrated qisnan.f, whose renamed stem may collide with an
            # ``independent`` module-procedure name) is likewise
            # precision-specific by construction — see ``rename_targets``.
            elif f.stem in _la_own or stem in rename_targets:
                precision_files.append(rel)
            # ``copy_files`` entries are precision-independent by
            # contract (the file is staged verbatim, no prefix rename).
            # The symbol scanner may never have visited them — e.g. a
            # Fortran ``copy_files`` entry in a C recipe — so they
            # won't appear in ``independent``. Treat them as common
            # explicitly.
            elif stem in independent or stem in config.copy_files:
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

        # Content-driven archive-filename prefix. The classifier's family
        # slots carry an 'R' (real) / 'C' (complex) type tag; the emitted
        # archive is named after what it actually contains, independent of
        # any filename convention:
        #   real families only    → e / q / m
        #   complex families only → y / x / w
        #   both                  → ey / qx / mw   (e.g. libeylapack.a)
        #   neither (precision-independent, e.g. scalapack_tools) → no prefix
        # This is the archive FILE name only; the CMake target name and the
        # exported linker symbols stay on LIB_PREFIX, so nothing links-breaks.
        _tags = {
            slot.type_tag
            for fam in classification.families
            for slot in fam.slots
        }
        archive_prefix = (
            (lib_prefix if 'R' in _tags else '')
            + (lib_prefix_complex if 'C' in _tags else '')
        )

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

set({lib_name}_ARCHIVE_PREFIX {archive_prefix})
"""
        (lib_dir / 'manifest.cmake').write_text(manifest)
        print(f'  Manifest: {len(common_files)} common, '
              f'{len(precision_files)} precision files')
        staged.append(lib_name)

    # Copy MF helper modules into staging so it's self-contained
    needs_mf = target_mode.module_name is not None
    # kind16 (REAL(KIND=16) / __float128) keeps the *standard* MPI datatypes
    # (MPI_REAL16 / MPI_COMPLEX32) but routes every reduction through custom
    # combine ops, because Intel MPI's builtin MPI_SUM/MAX/MIN have no
    # 16-byte-real kernel. Those ops live in a small plain-C helper
    # (quad_mpi.c) plus a Fortran shim (quad_mpi_f.f90) -- the quad analogue
    # of the multifloats bridge, but far simpler (no FetchContent, no C++,
    # no custom datatypes). A KIND-based target with a custom MPI module and
    # no multifloats module_name is exactly kind16.
    needs_quad_mpi = (
        target_mode.c_mpi_module is not None and not needs_mf
    )

    # Re-stage with a subset of --libraries must not shrink the unified
    # build's library list. Read STAGED_LIBRARIES from any prior
    # target_config.cmake, keep entries whose lib_dir still exists on
    # disk, and union them with this run's freshly-staged set so the
    # rewritten config reflects everything currently present in the
    # staging tree.
    prior_staged: list[str] = []
    prior_config = staging_dir / 'target_config.cmake'
    if prior_config.exists():
        m = re.search(
            r'^\s*set\s*\(\s*STAGED_LIBRARIES\s+([^)]*)\)',
            prior_config.read_text(),
            re.MULTILINE,
        )
        if m:
            for tok in m.group(1).replace(';', ' ').split():
                tok = tok.strip().strip('"')
                if tok and (staging_dir / tok).is_dir():
                    prior_staged.append(tok)
    merged: list[str] = []
    for n in prior_staged + staged:
        if n not in merged:
            merged.append(n)
    staged_list = ';'.join(merged)

    helpers_src = proj_root / 'recipes' / 'lapack' / 'mf_helpers'
    helpers_dst = staging_dir / '_helpers'
    if needs_mf:
        helpers_dst.mkdir(exist_ok=True)
        for name in ['la_constants_mw.f90', 'la_xisnan_mw.f90']:
            src = helpers_src / name
            if src.exists():
                shutil.copy2(src, helpers_dst / name)
        # Copy the first-party multifloats-mpi bridge library wholesale.
        # It is a standalone library under runtime/ with its own
        # CMakeLists.txt (add_subdirectory'd from cmake/CMakeLists.txt);
        # staging just plants the whole directory so the relative
        # add_subdirectory(runtime/multifloats-mpi) resolves in the staged
        # tree. The MPICH_SKIP_MPICXX / OMPI_SKIP_MPICXX guard that used to
        # be spliced in here is now baked into multifloats_bridge.h itself.
        mf_local = proj_root / 'runtime' / 'multifloats-mpi'
        if mf_local.is_dir():
            mf_dst = staging_dir / 'runtime' / 'multifloats-mpi'
            if mf_dst.exists():
                shutil.rmtree(mf_dst)
            shutil.copytree(mf_local, mf_dst)

    if needs_quad_mpi:
        # Copy the first-party quad-mpi bridge library wholesale (custom
        # MPI reduce ops on the standard MPI_REAL16 / MPI_COMPLEX32). Like
        # multifloats-mpi it is a standalone library under runtime/ with
        # its own CMakeLists.txt, add_subdirectory'd from the parent build.
        quad_local = proj_root / 'runtime' / 'quad-mpi'
        if quad_local.is_dir():
            quad_dst = staging_dir / 'runtime' / 'quad-mpi'
            if quad_dst.exists():
                shutil.rmtree(quad_dst)
            shutil.copytree(quad_local, quad_dst)

    target_config = f"""\
# Generated by: python -m migrator stage --target {target_mode.name}
set(TARGET_NAME "{target_mode.name}")
set(LIB_PREFIX "{lib_prefix}")
set(LIB_PREFIX_COMPLEX "{lib_prefix_complex}")
set(NEEDS_MULTIFLOATS {'TRUE' if needs_mf else 'FALSE'})
set(NEEDS_QUAD_MPI {'TRUE' if needs_quad_mpi else 'FALSE'})
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
    # INSTALL/ directory provides {s,d}lamch.f / {s,d}roundup_lwork.f,
    # which LAPACK SRC routines call but which aren't in SRC/ itself —
    # copy them into _reflapack_src/ alongside the SRC contents so a
    # single glob compiles the full reference. Both the single- and
    # double-precision variants are needed: the genuine standard ``lapack``
    # archive carries every arithmetic (the genuine single-precision
    # solvers smumps/cmumps call slamch_ / sroundup_lwork_, the double
    # ones dlamch_ / droundup_lwork_). These verbatim upstream files are
    # unpromoted here — only the migrated LAPACK archive elides
    # roundup_lwork via _strip_roundup_lwork; the standard one keeps it.
    netlib_lapack_src = proj_root / 'external' / 'lapack-3.12.1' / 'SRC'
    if netlib_lapack_src.is_dir():
        reflapack_dst = staging_dir / '_reflapack_src'
        if reflapack_dst.exists():
            shutil.rmtree(reflapack_dst)
        shutil.copytree(netlib_lapack_src, reflapack_dst)
        install_src = proj_root / 'external' / 'lapack-3.12.1' / 'INSTALL'
        for fname in ('slamch.f', 'dlamch.f',
                      'sroundup_lwork.f', 'droundup_lwork.f'):
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
        # PORD nested-dissection ordering — ships in-tree with MUMPS and
        # is self-contained standard C (no MPI / external dep). Staging
        # its algorithm sources (PORD/lib) + headers (PORD/include) lets
        # cmake build ``libpord`` and define ``-Dpord`` so ICNTL(7)=4
        # works; without it mumps_pord.c compiles as an inert stub.
        # Precision-agnostic (permutes the integer adjacency graph), so a
        # single build serves every migrated arithmetic.
        ('_mumps_pord_src',         'MUMPS_5.8.2/PORD/lib'),
        ('_mumps_pord_include',     'MUMPS_5.8.2/PORD/include'),
        # METIS 5.1.0 nested-dissection / k-way ordering — vendored under
        # external/metis-5.1.0 with every public API symbol privately
        # namespaced (METIS_<X> → METIS_MUMPS_<X>, internal libmetis__ →
        # libmetis_MUMPS_) so this copy can never clash with a system
        # METIS at link time; the MUMPS caller sites were renamed to
        # match. Staging GKlib + libmetis sources and the public header
        # lets cmake build ``libmetis`` and define ``-Dmetis`` so
        # ICNTL(7)=5 works; without it the mumps_metis*.c compile as inert
        # stubs. Integer-graph only, so a single 32-bit-idx build serves
        # every migrated arithmetic.
        ('_mumps_metis_gklib',      'metis-5.1.0/GKlib'),
        ('_mumps_metis_lib',        'metis-5.1.0/libmetis'),
        ('_mumps_metis_include',    'metis-5.1.0/include'),
        # Scotch 7.0.4 sequential ordering (ICNTL(7)=3) — vendored under
        # external/scotch-7.0.4, built with -DSCOTCH_NAME_SUFFIX=_mumps so
        # every public SCOTCH_* and internal _SCOTCH* symbol carries a
        # _mumps suffix and this copy can never clash with a system Scotch
        # at link time; the MUMPS caller sites resolve to the suffixed
        # names through scotch_rename_mumps.h. The bison/flex parser and
        # scotch.h/scotchf.h are pre-generated and vendored, so the build
        # needs no bison/flex. Staging libscotch + esmumps sources and the
        # generated headers lets cmake build ``scotch``/``scotcherr``/
        # ``esmumps`` and define ``-Dscotch`` so ICNTL(7)=3 works; without
        # it the mumps_scotch*.c compile as inert stubs. Integer-graph
        # only, so a single build serves every migrated arithmetic.
        ('_mumps_scotch_libsrc',    'scotch-7.0.4/libscotch'),
        ('_mumps_scotch_esmumps',   'scotch-7.0.4/esmumps'),
        ('_mumps_scotch_include',   'scotch-7.0.4/include'),
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
    print('\nTo build:')
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
        ('_mumps_pord_src',       'MUMPS_5.8.2/PORD/lib',              None),
        ('_mumps_pord_include',   'MUMPS_5.8.2/PORD/include',          None),
        ('_mumps_metis_gklib',    'metis-5.1.0/GKlib',                 None),
        ('_mumps_metis_lib',      'metis-5.1.0/libmetis',              None),
        ('_mumps_metis_include',  'metis-5.1.0/include',               None),
        ('_mumps_scotch_libsrc',  'scotch-7.0.4/libscotch',            None),
        ('_mumps_scotch_esmumps', 'scotch-7.0.4/esmumps',              None),
        ('_mumps_scotch_include', 'scotch-7.0.4/include',              None),
    ]
    for dst_name, rel_src, recipe_name in _std_dirs:
        _stage_dst(dst_name, _staged_or_external(rel_src, recipe_name))

    print(f'  Output:  {staging_dir}')
    print('\nTo build:')
    print(f'  cmake -S {staging_dir} -B {staging_dir}/build -DCMAKE_BUILD_TYPE=Release')
    print(f'  cmake --build {staging_dir}/build -j')
    print(f'  ctest --test-dir {staging_dir}/build')
