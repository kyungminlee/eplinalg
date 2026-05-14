"""Recipe configuration loader.

Recipes are YAML files that describe a library's source layout and
migration parameters. The engine reads a recipe and performs the
appropriate migration.
"""

import sys
from dataclasses import dataclass, field
from pathlib import Path

try:
    import yaml
except ImportError:
    yaml = None  # type: ignore[assignment]


# Top-level keys recognized by load_recipe(). Unknown keys are warned on
# load so that a typo (``copy-files`` vs ``copy_files``) doesn't silently
# default to empty and produce a partially-migrated tree.
_KNOWN_RECIPE_KEYS: frozenset[str] = frozenset({
    'library', 'language', 'source_dir', 'extensions',
    'symbols', 'prefix',
    'skip_files', 'copy_files', 'prefer_source',
    'module_renames', 'extra_renames',
    'copy_all_originals', 'patches',
    'depends', 'extra_symbol_dirs',
    'extra_migrate_files', 'extra_c_dirs', 'extra_fortran_dirs',
    'keep_kind_manifest',
    'c_return_types', 'c_type_aliases', 'c_pointer_cast_aliases',
    'header_patches', 'overrides',
    'expected_divergences', 'defer_all_divergences',
    'asymmetric_patches', 'one_sided_cleanup',
})


@dataclass
class RecipeConfig:
    """Configuration for a single library migration."""
    library: str
    language: str                     # "fortran" or "c"
    source_dir: Path
    extensions: list[str]
    library_path: Path | None = None
    skip_files: set[str] = field(default_factory=set)
    copy_files: set[str] = field(default_factory=set)  # Copy unchanged (multi-precision utilities)
    # Source stems (uppercase, no extension) whose migrated output
    # should win as canonical, overriding the default D/Z-first
    # preference. Used to route around upstream bugs that live only
    # in the D or Z half of a precision pair (e.g. ScaLAPACK's
    # PZUNGQL / PZUNML2 call PB_TOPGET where they should call
    # PB_TOPSET; PCUNGQL / PCUNML2 have the correct restore).
    prefer_source: set[str] = field(default_factory=set)
    copy_all_originals: bool = False  # For C: copy all files, then add clones
    patches: list[str] = field(default_factory=list)
    depends: list[Path] = field(default_factory=list)  # Dependency recipe paths
    extra_symbol_dirs: list[Path] = field(default_factory=list)  # Extra dirs to scan for symbols
    # Individual files (outside ``source_dir``) to migrate in the same
    # pass as ``source_dir``. Each entry is a project-root-relative
    # path to a single ``.f``/``.f90``/``.F90``/``.c``/``.h`` file.
    # Used to pull in targeted leaf sources from shared directories
    # whose other contents belong to a different library — e.g.
    # LAPACK migrates ``INSTALL/dlamch.f`` without swallowing the
    # timer variants and test programs that live alongside it
    # (``INSTALL/droundup_lwork.f`` is deliberately not migrated;
    # the engine's ``_strip_roundup_lwork`` post-pass elides every
    # call site — see ``recipes/lapack.yaml``); PTZBLAS pulls in
    # ``TOOLS/zzdotc.f`` and ``TOOLS/zzdotu.f`` without claiming
    # the rest of ScaLAPACK's TOOLS/.
    extra_migrate_files: list[Path] = field(default_factory=list)
    # Additional C source directories to *migrate* (not just scan) in
    # the same generic-rename-map pass as ``source_dir``. Used by PBLAS
    # to pull in the PTOOLS/ helper sources alongside the SRC/ entry
    # points so the cloned ddgemm.c entry points have a real
    # PB_Cddtypeset implementation to call. Files are flat-copied into
    # ``output_dir`` (no subdirectory mirroring) so include resolution
    # stays simple.
    extra_c_dirs: list[Path] = field(default_factory=list)
    # Additional Fortran source directories whose files are migrated in
    # the same pass as ``source_dir``. Used by MUMPS to pull in the
    # per-arithmetic header files under ``external/MUMPS_5.8.2/include/``
    # (``dmumps_struc.h`` etc.), which are Fortran content despite the
    # ``.h`` extension and must be migrated so the ``INCLUDE`` statements
    # in ``dmumps_struc_def.F`` resolve against the renamed target file.
    # Files are flat-copied into ``output_dir`` (no subdir mirroring).
    extra_fortran_dirs: list[Path] = field(default_factory=list)
    # Additional C return types to recognize when scanning for function
    # definitions, as regex fragments (e.g. ``r'PBTYP_T\s*\*'``). Used
    # only when ``language == 'c'``; the default set in
    # ``symbol_scanner._C_DEFAULT_RETURN_TYPES`` is always included.
    c_return_types: list[str] = field(default_factory=list)
    # Extra library-specific C typedef renames applied after the
    # built-in double/float/SCOMPLEX/DCOMPLEX substitutions. Each entry
    # has ``from`` (list of source identifiers) and ``to`` (target
    # identifier). Both the ``to`` field and the inserted text in
    # ``header_patches`` support template substitution with the C
    # migrator's template_vars (``{REAL_TYPE}``, ``{COMPLEX_TYPE}``,
    # ``{C_REAL_TYPE}``, ``{RP}``, ``{CP}``, ``{RPU}``, ``{CPU}``).
    c_type_aliases: list[dict] = field(default_factory=list)
    # Pointer-cast aliases applied to copy-original C sources (those
    # that are precision-independent dispatchers, e.g. PB_Cconjg.c uses
    # ``((double*)CALPHA)[REAL_PART] = …`` switched on TYPE->type).
    # Each entry has ``from`` (list of full cast strings like
    # ``(double*)``) and ``to`` (the replacement, with template
    # substitution). Distinct from ``c_type_aliases`` because the
    # broad ``double → REAL_TYPE`` substitution would clobber the
    # SCPLX/DCPLX dispatch labels in those originals; pointer-cast
    # rewriting is needed for the kind-correct stride but the bare
    # ``double``/``float`` keywords must stay.
    c_pointer_cast_aliases: list[dict] = field(default_factory=list)
    # Insert new content into migrated headers after a literal anchor
    # line. Each entry: ``{'file': <relative path under source_dir>,
    # 'after': <anchor line>, 'insert': <text>}``. Used to define
    # library-specific extended-precision typedefs referenced by
    # c_type_aliases targets.
    header_patches: list[dict] = field(default_factory=list)
    # Target-gated verbatim file overrides. Structure:
    #
    #     overrides:
    #       <target_name>:
    #         src_dir: <path relative to recipe file>
    #         files:
    #           - <filename>
    #           - ...
    #
    # For the active target, each listed file is copied verbatim from
    # ``<recipe_dir>/<src_dir>/<filename>`` to ``<output_dir>/<filename>``
    # after the main C migration has produced clones and header patches,
    # so the override wins. Used for hand-written replacement kernels
    # that cannot be produced by regex substitution (e.g. BI_*vv* for
    # multifloats double-double arithmetic).
    overrides: dict = field(default_factory=dict)
    # Directory containing the recipe file, used to resolve paths in
    # ``overrides`` and similar recipe-relative references.
    recipe_dir: Path | None = None
    # Per-line "keep-kind" manifest: for each source filename (basename,
    # not stem — e.g. ``dana_aux.F``), the set of 1-based line numbers
    # whose ``DOUBLE PRECISION`` declarations must NOT be promoted.
    # Used by MUMPS, where ``DOUBLE PRECISION`` overloads "working
    # precision" and "arithmetic-agnostic DP" (timing, flop counters,
    # MPI_WTIME buffers, Scotch ABI). Generated by
    # ``scripts/mumps_sweep_keep_kind.py``. See ``recipes/README.md``.
    keep_kind_lines: dict[str, frozenset[int]] = field(default_factory=dict)
    # Post-migration module name rewrites applied to migrated Fortran
    # files only (copy_files / skip_files are untouched). Each entry
    # maps ``OLD_MODULE`` → ``NEW_MODULE``; both ``USE OLD_MODULE`` and
    # ``USE OLD_MODULE, ONLY: ...`` are rewritten. Used by MUMPS to
    # redirect migrated callers from the upstream ``MUMPS_MEMORY_MOD``
    # (kept verbatim via copy_files) to a hand-written
    # ``MUMPS_MEMORY_MOD_EP`` that adds extended-precision reallocators
    # without collapsing the original S/D/C/Z generic interface.
    module_renames: dict[str, str] = field(default_factory=dict)
    # Recipe-level forced rename entries, appended to the classifier's
    # rename map after family discovery. Used for precision-prefixed
    # symbols whose S/C sibling does not exist in the upstream source
    # (so the prefix classifier cannot pair them — ScaLAPACK's
    # ``pdlaiectb_/pdlaiectl_`` are the canonical example: the IEEE
    # big/little-endian Sturm-sequence helpers exist only for double
    # precision because the bit-shift sign trick they rely on is a
    # 64-bit-double-only hack). Each entry maps an upstream upper-cased
    # identifier to a target template that may reference {RP}/{CP}/
    # {RPU}/{CPU} via target template_vars. Applied to migrated output
    # (clones + caller bodies) but NOT to copy-original sources, so
    # the upstream (un-migrated) entry points keep their original
    # symbol names and link cleanly alongside the renamed clones.
    extra_renames: dict[str, str] = field(default_factory=dict)
    # Convergence-report whitelist. Each stem (uppercased, no extension)
    # names the canonical (D/Z) member of a co-family pair whose
    # divergence is expected — typically because the two upstream halves
    # genuinely differ (BLAS sdot line-swap, srotmg constants, MINRGP
    # tuning split between S and D, etc.) or because a patch covers only
    # one half by design. Pairs whose canonical stem appears here are
    # filtered out of the convergence report and do NOT cause CI to fail.
    # See doc/UPSTREAM_BUGS.md for individual entries.
    expected_divergences: set[str] = field(default_factory=set)
    # Coarse-grained whitelist: when True, every divergence in this
    # library is filtered out. Used for libraries where convergence is
    # currently dominated by migrator-internal asymmetries (PBLAS K&R
    # re-emergence, MUMPS kind-promotion, scalapack_c TYPE rename gap)
    # tracked separately. The fix is migrator-side; once the migrator
    # gap closes, switch the recipe back to enumerated
    # ``expected_divergences``.
    defer_all_divergences: bool = False
    # Patches under ``recipes/<lib>/patches/`` that touch only one half
    # of a co-family pair because the upstream sibling carries a
    # genuinely different bug shape (or no analogous bug). Listed here,
    # the symmetric-patch CI check (``migrator verify-patches``) skips
    # them; otherwise it fails when a precision-prefixed file is touched
    # without its siblings. Use this field for "real bug, S/C may need
    # its own future patch with a different fix" — periodically review
    # entries to see whether the sibling situation has changed.
    asymmetric_patches: list[str] = field(default_factory=list)
    # Patches that close a D↔S or Z↔C *cosmetic* asymmetry by stripping
    # upstream dead code (unused PARAMETER blocks, redundant CMPLX
    # casts, dead INTRINSIC entries, literal-style mismatches). The
    # sibling half is already in the post-patch shape — no S/C patch
    # is or will be needed. Same CI semantics as ``asymmetric_patches``
    # (skipped by symmetric check), but the separate field communicates
    # intent so reviewers don't waste time hunting for missing siblings.
    one_sided_cleanup: list[str] = field(default_factory=list)


def load_recipe(recipe_path: Path,
                project_root: Path | None = None) -> RecipeConfig:
    """Load a recipe YAML file.

    Args:
        recipe_path: Path to the .yaml recipe file.
        project_root: Project root for resolving relative paths.
            Defaults to recipe_path's grandparent.
    """
    if yaml is None:
        raise ImportError(
            'PyYAML is required to load recipe configs. '
            'Install it with: pip install pyyaml'
        )

    if project_root is None:
        # recipes/*.yaml → project root is two levels up
        project_root = recipe_path.parent.parent

    with open(recipe_path) as f:
        data = yaml.safe_load(f)

    if not isinstance(data, dict):
        raise ValueError(
            f'{recipe_path}: top-level YAML must be a mapping, '
            f'got {type(data).__name__}'
        )

    for required in ('library', 'language', 'source_dir'):
        if required not in data:
            raise KeyError(
                f'{recipe_path}: missing required recipe key {required!r}'
            )

    unknown = sorted(set(data) - _KNOWN_RECIPE_KEYS)
    if unknown:
        print(
            f'  warning: {recipe_path.name}: unknown recipe key(s) '
            f'{unknown!r}; ignored. (Known keys: {sorted(_KNOWN_RECIPE_KEYS)})',
            file=sys.stderr,
        )

    source_dir = project_root / data['source_dir']

    library_path = None
    symbols_cfg = data.get('symbols') or {}
    if symbols_cfg.get('library_path'):
        library_path = project_root / symbols_cfg['library_path']

    skip = set(s.upper() for s in data.get('skip_files', []))
    copy = set(s.upper() for s in data.get('copy_files', []))
    prefer = set(s.upper() for s in data.get('prefer_source', []))

    # Resolve dependency recipe paths relative to the recipe directory
    depends_raw = data.get('depends', [])
    depends = [recipe_path.parent / d for d in depends_raw]

    # Resolve extra symbol directories relative to the project root
    extra_dirs_raw = data.get('extra_symbol_dirs', [])
    extra_symbol_dirs = [project_root / d for d in extra_dirs_raw]

    # Load the keep-kind manifest if specified. The manifest is a plain
    # text file with one ``<path>:<lineno>`` entry per line; the path is
    # ignored (we key by basename) and the lineno is 1-based.
    keep_kind_lines: dict[str, set[int]] = {}
    manifest_rel = data.get('keep_kind_manifest')
    if manifest_rel:
        manifest_path = recipe_path.parent / manifest_rel
        for entry in manifest_path.read_text().splitlines():
            entry = entry.strip()
            if not entry or entry.startswith('#'):
                continue
            path_str, _, lineno_str = entry.rpartition(':')
            if not path_str or not lineno_str:
                continue
            basename = Path(path_str).name
            keep_kind_lines.setdefault(basename, set()).add(int(lineno_str))
    keep_kind_frozen = {k: frozenset(v) for k, v in keep_kind_lines.items()}

    return RecipeConfig(
        library=data['library'],
        language=data['language'],
        source_dir=source_dir,
        extensions=[e.lower() for e in data.get('extensions', ['.f', '.f90'])],
        library_path=library_path,
        skip_files=skip,
        copy_files=copy,
        prefer_source=prefer,
        copy_all_originals=data.get('copy_all_originals', False),
        patches=data.get('patches', []),
        depends=depends,
        extra_symbol_dirs=extra_symbol_dirs,
        extra_migrate_files=[project_root / p
                             for p in (data.get('extra_migrate_files') or [])],
        c_return_types=list(data.get('c_return_types', [])),
        c_type_aliases=list(data.get('c_type_aliases', [])),
        c_pointer_cast_aliases=list(data.get('c_pointer_cast_aliases', [])),
        header_patches=list(data.get('header_patches', [])),
        overrides=dict(data.get('overrides') or {}),
        recipe_dir=recipe_path.parent,
        extra_c_dirs=[project_root / d
                      for d in (data.get('extra_c_dirs') or [])],
        extra_fortran_dirs=[project_root / d
                            for d in (data.get('extra_fortran_dirs') or [])],
        keep_kind_lines=keep_kind_frozen,
        module_renames={
            str(k).upper(): str(v).upper()
            for k, v in (data.get('module_renames') or {}).items()
        },
        extra_renames={
            str(k).upper(): str(v)
            for k, v in (data.get('extra_renames') or {}).items()
        },
        expected_divergences={
            str(s).upper() for s in (data.get('expected_divergences') or [])
        },
        defer_all_divergences=bool(data.get('defer_all_divergences', False)),
        asymmetric_patches=list(data.get('asymmetric_patches') or []),
        one_sided_cleanup=list(data.get('one_sided_cleanup') or []),
    )
