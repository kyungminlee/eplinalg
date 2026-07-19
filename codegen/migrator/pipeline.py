"""Migration pipeline — orchestrates the full migration for any library.

Usage:
    from migrator.pipeline import run_migration
    run_migration(recipe_path, output_dir, target_mode=16)
"""

import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from .config import RecipeConfig, load_recipe
from .prepare import prepare_recipe


def _build_module_rename_pairs(
    config: RecipeConfig,
) -> list[tuple[re.Pattern[str], str]]:
    """Compile the recipe's ``module_renames:`` map into regex pairs.

    The pattern matches ``USE <OLD_MODULE>`` (case-insensitive, both
    bare ``USE`` and ``USE …, ONLY: …`` forms) and rewrites the module
    name to the new one. Same shape used by ``run_fortran_migration``
    when emitting the canonical text — surfacing this helper lets
    convergence re-migration apply the identical rule, keeping the
    on-disk canonical and the re-migrated sibling on equal footing.
    """
    pairs: list[tuple[re.Pattern[str], str]] = []
    for old_mod, new_mod in (config.module_renames or {}).items():
        pat = re.compile(r'(?i)(\bUSE\s+)' + re.escape(old_mod) + r'\b')
        pairs.append((pat, r'\g<1>' + new_mod))
    return pairs


def _apply_module_renames(
    text: str,
    pairs: list[tuple[re.Pattern[str], str]],
) -> str:
    for pat, repl in pairs:
        text = pat.sub(repl, text)
    return text


def _apply_call_arg_casts(
    text: str,
    casts: list[tuple[str, str]],
) -> str:
    """Wrap literal call-argument expressions in a cast, post-migration.

    Each ``(find, wrap)`` rewrites every occurrence of the literal
    ``find`` to ``wrap(find)``. Runs after the intrinsic rewriter (so
    the injected ``dble(...)`` is not itself promoted to
    ``REAL(..., KIND=N)``) and only on migrated output — copy_files
    verbatim sources are never touched. Idempotent guard: a ``find``
    already immediately preceded by ``wrap(`` is left alone so a
    re-run does not double-wrap. See
    :attr:`RecipeConfig.call_arg_casts`.
    """
    for find, wrap in casts:
        already = f'{wrap}({find}'
        if already in text:
            continue
        text = text.replace(find, f'{wrap}({find})')
    return text


from .symbol_scanner import scan_symbols
from .prefix_classifier import classify_symbols
from .divergence import _canonicalize_for_compare, _strip_fortran_comments
from .fortran.decls import scan_type_component_names
from .fortran_migrator import (
    migrate_file_to_string,
    target_filename,
)
from .c_migrator import migrate_c_directory
from .templates import (
    build_sub_vars as _build_sub_vars,
    expand_template as _expand_template,
)
from .privatize import privatize_tree
from .target_mode import TargetMode

from tqdm import tqdm


def _apply_extra_renames(rename_map: dict[str, str],
                         config: RecipeConfig,
                         target_mode: TargetMode) -> dict[str, str]:
    """Append recipe-declared ``extra_renames`` to ``rename_map``.

    Targets may use ``{RP}/{CP}/{RPU}/{CPU}`` template substitutions.
    Returns ``rename_map`` mutated in place. See
    :class:`RecipeConfig.extra_renames` for the canonical use case
    (precision-prefixed orphan symbols whose S/C sibling does not
    exist upstream so the classifier cannot pair them).
    """
    if not config.extra_renames:
        return rename_map
    template_vars = _build_sub_vars(target_mode)
    for src_upper, tgt_template in config.extra_renames.items():
        rename_map[src_upper] = _expand_template(tgt_template, template_vars).upper()
    return rename_map


def _build_component_oracle(
    paths, target_mode: TargetMode,
) -> tuple[frozenset[str], frozenset[str]]:
    """Global derived-type component-type oracle for formatted-output
    narrowing.

    Harvests real/complex component names from every ``TYPE...END TYPE``
    block across ``paths`` (the whole staged set), unions them, and drops
    the real∩complex intersection. That intersection is the family of
    ambiguous data-array components (``A``, ``S``, ``RHS``, ``SCHUR`` — real
    in the d/s struct, complex in the c/z struct) which are never emitted
    through a single real edit descriptor and so must not be narrowed. What
    survives in ``comp_real`` are the statistics that are real in *every*
    arithmetic (``CNTL``, ``RINFO``, ``RINFOG``, ``DKEEP``) — exactly the
    fields printed as ``id%DKEEP(160)`` in the diagnostics, which a per-file
    oracle cannot see because the struct is defined in a ``USE``d module.

    Component names are only ever consulted for ``%``-qualified references
    (see :func:`io_narrow._narrow_items`), so a bare local of the same name
    is unaffected. Returns empty sets for kind-based targets (real64x2 /
    cmplx64x2 narrowing does not apply there)."""
    if target_mode is None or target_mode.is_kind_based:
        return frozenset(), frozenset()
    comp_real: set[str] = set()
    comp_complex: set[str] = set()
    for p in paths:
        try:
            text = p.read_text(errors='replace')
        except OSError:
            continue
        r, c = scan_type_component_names(text)
        comp_real |= r
        comp_complex |= c
    ambiguous = comp_real & comp_complex
    return frozenset(comp_real - ambiguous), frozenset(comp_complex - ambiguous)


def classify_recipe_symbols(config: RecipeConfig):
    """Scan the recipe's own sources and classify precision families."""
    symbols = scan_symbols(config.source_dir, config.language,
                           config.extensions,
                           extra_c_return_types=tuple(config.c_return_types))
    return classify_symbols(symbols)


def _canonical_rank(stem: str, prefer: frozenset[str]) -> int:
    # 0 = highest priority (recipe-pinned), 1 = D/Z default,
    # 2 = S/C fallback. Sorting ascending makes the preferred
    # source the first writer for each target.
    u = stem.upper()
    if u in prefer:
        return 0
    if u and u[0] in ('D', 'Z'):
        return 1
    return 2


def _migrate_parallel(paths, rename_map: dict[str, str],
                      target_mode: TargetMode,
                      parser: str | None, parser_cmd: str | None,
                      config: RecipeConfig,
                      ) -> dict[Path, tuple[str, str] | None]:
    """Migrate ``paths`` in a process pool; return per-path results.

    Each worker runs parser + regex substitution, which is the dominant
    cost for large libraries like LAPACK. A crash on one file is
    reported as a warning and recorded as ``None`` so the caller can
    skip it without losing the rest of the batch.
    """
    workers = max(1, (os.cpu_count() or 4))
    comp_real, comp_complex = _build_component_oracle(paths, target_mode)
    results: dict[Path, tuple[str, str] | None] = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = {
            ex.submit(migrate_file_to_string, p, rename_map, target_mode,
                      parser, parser_cmd,
                      config.keep_kind_lines.get(p.name),
                      comp_real, comp_complex): p
            for p in paths
        }
        for fut in tqdm(as_completed(futures), total=len(futures),
                        desc='  Migrating', unit='file',
                        mininterval=1.0, miniters=10):
            p = futures[fut]
            try:
                results[p] = fut.result()
            except Exception as exc:
                print(f'  warning: migration crashed on {p.name}: '
                      f'{type(exc).__name__}: {exc}', file=sys.stderr)
                results[p] = None
    return results


def _postprocess_migrated(migrated: str, src_path: Path,
                          config: RecipeConfig,
                          module_rename_pairs) -> str:
    """Post-migration text fixups: module renames + call-arg casts."""
    migrated = _apply_module_renames(migrated, module_rename_pairs)
    casts = config.call_arg_casts.get(src_path.name)
    if casts:
        migrated = _apply_call_arg_casts(migrated, casts)
    return migrated


def run_fortran_migration(config: RecipeConfig, rename_map: dict[str, str],
                          output_dir: Path, target_mode: TargetMode,
                          dry_run: bool = False,
                          classification=None,
                          parser: str | None = None,
                          parser_cmd: str | None = None) -> dict:
    """Run Fortran migration pipeline."""
    # Identify precision-independent symbols
    if classification is None:
        classification = classify_recipe_symbols(config)
    independent = classification.independent

    migrated_count = 0
    copied_count = 0
    skipped: list[str] = []
    divergences: list[str] = []

    # Collect eligible source files. In addition to files with the
    # recipe's declared extensions, include any file whose stem is listed
    # in copy_files (regardless of extension) — this lets libraries
    # carry shared, type-independent headers (e.g. MUMPS's
    # mumps_tags.h / mumps_headers.h) through without having to widen
    # the extensions list and accidentally migrate them.
    #
    # extra_fortran_dirs contributes additional directories whose files
    # are migrated alongside source_dir. MUMPS uses this for
    # per-arithmetic headers (``dmumps_struc.h`` → ``qmumps_struc.h``)
    # that live in a parallel ``include/`` directory but are Fortran
    # content that must go through the full migration pipeline.
    src_dirs = [config.source_dir] + list(
        getattr(config, 'extra_fortran_dirs', []) or []
    )
    src_files = sorted(
        p for d in src_dirs for p in d.iterdir()
        if p.suffix.lower() in config.extensions
        or p.stem.upper() in config.copy_files
    )
    # extra_migrate_files pulls in specific leaf files from outside
    # source_dir / extra_fortran_dirs — used to target individual
    # helpers that live in a shared directory whose other contents
    # belong to a different library (LAPACK's INSTALL/dlamch.f,
    # PTZBLAS's TOOLS/zzdotc.f). The symbol scanner already sees these
    # via extra_symbol_dirs, so the rename map covers the referenced
    # names; this line only widens the emit loop to include them.
    extra_migrate = list(getattr(config, 'extra_migrate_files', []) or [])
    src_files = sorted(set(src_files).union(
        p for p in extra_migrate if p.is_file()
    ))
    # Convergence buffer: first writer of each output name stores its
    # text; subsequent writers must agree or we record a divergence.
    # D/Z sources are preferred as the canonical text: when a pair
    # (SGEMM, DGEMM) both target QGEMM, DGEMM's migrated body is kept
    # and SGEMM's is only consulted for the equality check. The
    # ``prefer_source`` recipe field flips this for individual stems
    # (e.g. to route around a bug that only exists in the D/Z half).
    prefer = config.prefer_source
    src_files.sort(key=lambda p: (_canonical_rank(p.stem, prefer), p.name))

    canonical_text: dict[str, str] = {}
    canonical_normalized: dict[str, str] = {}
    canonical_source: dict[str, str] = {}

    # Build module-rename regex pairs once. Applied post-migration to
    # every migrated file (copy_files are deliberately untouched so the
    # verbatim upstream module keeps its original name).
    module_rename_pairs = _build_module_rename_pairs(config)

    # Partition: copy/skip/dry-run decisions stay in the main process;
    # only the Flang-bound migration is dispatched to workers.
    to_migrate: list[Path] = []
    for src_path in src_files:
        stem_upper = src_path.stem.upper()
        if stem_upper in config.skip_files:
            skipped.append(src_path.name)
            continue
        is_copy = stem_upper in config.copy_files or stem_upper in independent
        if dry_run:
            if is_copy:
                tqdm.write(f'  {src_path.name} (copy)')
            else:
                out_name = target_filename(src_path.name, rename_map, target_mode)
                tqdm.write(f'  {src_path.name} → {out_name}')
            continue
        if is_copy:
            (output_dir / src_path.name).write_text(
                src_path.read_text(errors='replace'))
            copied_count += 1
            continue
        to_migrate.append(src_path)

    if dry_run:
        return {
            'migrated': 0, 'copied': copied_count,
            'skipped': skipped, 'divergences': divergences,
        }

    # Parallel migration. We reduce results in canonical-first order so
    # D/Z output is the one written to disk.
    results = _migrate_parallel(to_migrate, rename_map, target_mode,
                                parser, parser_cmd, config)

    # Reduce in deterministic D/Z-first order.
    for src_path in to_migrate:
        result = results.get(src_path)
        if result is None:
            skipped.append(src_path.name)
            continue
        out_name, migrated = result
        migrated = _postprocess_migrated(migrated, src_path, config,
                                         module_rename_pairs)
        normalized = _canonicalize_for_compare(
            _strip_fortran_comments(migrated, src_path.suffix)
        )

        prior = canonical_normalized.get(out_name)
        if prior is None:
            # First (canonical) writer for this target name
            (output_dir / out_name).write_text(migrated)
            canonical_text[out_name] = migrated
            canonical_normalized[out_name] = normalized
            canonical_source[out_name] = src_path.name
            migrated_count += 1
        elif prior == normalized:
            # Convergence: co-family member produced identical code
            # (comments may differ and are ignored).
            pass
        else:
            # Divergence: co-family members disagree. Keep the canonical
            # (D/Z) version on disk and record the mismatch.
            divergences.append(
                f'{src_path.name} vs {canonical_source[out_name]} '
                f'→ {out_name}'
            )

    return {
        'migrated': migrated_count,
        'copied': copied_count,
        'skipped': skipped,
        'divergences': divergences,
    }


def _scan_extra_dirs(extra_dirs: list[Path],
                     extra_c_return_types: tuple[str, ...]) -> set[str]:
    """Scan ``extra_symbol_dirs`` plus one level of subdirectories for
    Fortran and C symbols."""
    out: set[str] = set()
    for extra_dir in extra_dirs:
        if not extra_dir.is_dir():
            continue
        for lang, exts in [('fortran', ['.f', '.f90', '.F90']),
                           ('c', ['.c'])]:
            out |= scan_symbols(
                extra_dir, lang, exts,
                extra_c_return_types=extra_c_return_types)
        for sub in sorted(extra_dir.iterdir()):
            if sub.is_dir():
                for lang, exts in [('fortran', ['.f', '.f90', '.F90']),
                                   ('c', ['.c'])]:
                    out |= scan_symbols(
                        sub, lang, exts,
                        extra_c_return_types=extra_c_return_types)
    return out


def _collect_all_symbols(config: RecipeConfig,
                         project_root: Path | None) -> set[str]:
    """Scan the recipe plus every transitive dependency — including
    each dependency's ``extra_symbol_dirs`` — so the rename map
    covers kernels a dependency exposes indirectly (e.g. LAPACK's
    INSTALL/slamch.f reached via scalapack → lapack)."""
    own_c_types = tuple(config.c_return_types)
    symbols = scan_symbols(
        config.source_dir, config.language,
        config.extensions,
        extra_c_return_types=own_c_types,
    )
    symbols |= _scan_extra_dirs(config.extra_symbol_dirs, own_c_types)

    seen: set[Path] = set()
    queue = list(config.depends)
    while queue:
        dep_path = queue.pop(0)
        key = dep_path.resolve()
        if key in seen:
            continue
        seen.add(key)
        dep_cfg = load_recipe(dep_path, project_root)
        dep_c_types = tuple(dep_cfg.c_return_types)
        symbols |= scan_symbols(
            dep_cfg.source_dir, dep_cfg.language,
            dep_cfg.extensions,
            extra_c_return_types=dep_c_types,
        )
        symbols |= _scan_extra_dirs(dep_cfg.extra_symbol_dirs, dep_c_types)
        queue.extend(dep_cfg.depends)
    return symbols


def run_c_migration(config: RecipeConfig, output_dir: Path,
                    target_mode: TargetMode, dry_run: bool = False,
                    classification=None,
                    rename_map: dict[str, str] | None = None) -> dict:
    """Run C migration pipeline (clone-and-substitute).

    When `classification` and `rename_map` are supplied, the generic
    rename-map-driven migration is used (for ScaLAPACK-style libraries
    like PBLAS). Otherwise the BLACS-specific path is used.
    """
    if dry_run:
        print('  WARNING: --dry-run is not supported for C migration; '
              'no C files were analyzed. Rerun without --dry-run for '
              'the real migration.')
        return {'cloned': [], 'template_vars': {}}

    overrides = _resolve_overrides(config, target_mode)

    result = migrate_c_directory(
        config.source_dir, output_dir, target_mode,
        classification=classification,
        rename_map=rename_map,
        c_type_aliases=config.c_type_aliases,
        c_pointer_cast_aliases=config.c_pointer_cast_aliases,
        header_patches=config.header_patches,
        overrides=overrides,
        extra_c_dirs=config.extra_c_dirs,
        skip_files=config.skip_files,
        copy_files=config.copy_files,
    )
    return result


def _apply_fortran_overrides(config: RecipeConfig,
                             target_mode: TargetMode,
                             output_dir: Path) -> None:
    """Copy target-gated override files into the Fortran output dir.

    Used for hand-written shared modules that cannot be produced by
    migration (e.g. MUMPS's extended-precision reallocator module).
    """
    overrides = _resolve_overrides(config, target_mode)
    for src_path, dst_name in overrides:
        (output_dir / dst_name).write_text(src_path.read_text())
        print(f'  Override:  {dst_name} (from {src_path.name})')


def _resolve_overrides(config: RecipeConfig,
                       target_mode: TargetMode) -> list[tuple[Path, str]]:
    """Select the active target's override entries from the recipe.

    The ``overrides`` recipe field is a dict keyed by target name. For
    the currently active target, return a list of ``(src_path, dst_name)``
    pairs with ``src_path`` resolved against the recipe directory.
    Returns an empty list if the recipe has no overrides for this target.
    """
    target_entry = (config.overrides or {}).get(target_mode.name)
    if not target_entry:
        return []
    src_dir_rel = target_entry.get('src_dir', '')
    files = target_entry.get('files', [])
    if config.recipe_dir is None:
        raise RuntimeError(
            'RecipeConfig.recipe_dir is not set; cannot resolve overrides'
        )
    src_dir = (config.recipe_dir / src_dir_rel).resolve()
    return [(src_dir / fname, fname) for fname in files]


def _print_divergences(result: dict, limit: int = 10) -> None:
    """Print the divergence summary block of a migration result."""
    if not result.get('divergences'):
        return
    print(f'  Divergences: {len(result["divergences"])} '
          f'(co-family members produced differing output)')
    for d in result['divergences'][:limit]:
        print(f'    {d}')
    if len(result['divergences']) > limit:
        print(f'    ... ({len(result["divergences"]) - limit} more)')


def run_migration(recipe_path: Path, output_dir: Path,
                  target_mode=None, dry_run: bool = False,
                  project_root: Path | None = None,
                  parser: str | None = None,
                  parser_cmd: str | None = None) -> dict:
    """Run the full migration pipeline for a library.

    Args:
        recipe_path: Path to the YAML recipe file.
        output_dir: Where to write migrated files.
        target_mode: 10 or 16.
        dry_run: If True, show what would be done without writing.
        project_root: Project root for resolving relative paths.
        parser: Parse tree backend: ``'flang'``, ``'gfortran'``, or
            ``None`` (regex-only).
        parser_cmd: Explicit path to the compiler binary.

    Returns:
        Summary dict with migration statistics.
    """
    config = prepare_recipe(recipe_path, project_root)
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f'Library:     {config.library}')
    print(f'Language:    {config.language}')
    print(f'Source:      {config.source_dir}')
    print(f'Target:      {target_mode.name}')
    print()

    # Scan symbols and classify precision families
    print('Scanning symbols...')
    own_symbols = scan_symbols(
        config.source_dir, config.language,
        config.extensions,
        extra_c_return_types=tuple(config.c_return_types),
    )
    own_count = len(own_symbols)
    symbols = _collect_all_symbols(config, project_root)

    classification = classify_symbols(symbols)
    rename_map = classification.build_rename_map(target_mode)
    rename_map = _apply_extra_renames(rename_map, config, target_mode)

    # NOTE: target_mode.known_constants (ZERO/ONE/...) are handled
    # per-file by strip_known_constants_from_decls +
    # replace_known_constants. We intentionally do NOT add them to the
    # global rename_map: doing so would also rewrite the LHS aliases
    # of LAPACK ``USE LA_CONSTANTS, ONLY: zero=>dzero`` clauses (and
    # rename the lowercase ``zero`` reference everywhere) which breaks
    # the alias binding.


    print(f'  {own_count} own symbols + {len(symbols) - own_count} from dependencies')
    print(f'  {len(classification.families)} precision families')
    print(f'  {len(classification.independent)} independent symbols')
    print(f'  {len(rename_map)} renames computed')

    # Dispatch to language-specific migrator
    print(f'\nMigrating to {target_mode.name}...')

    if config.language == 'fortran':
        result = run_fortran_migration(
            config, rename_map, output_dir, target_mode, dry_run,
            classification=classification,
            parser=parser, parser_cmd=parser_cmd,
        )
        if not dry_run:
            _apply_fortran_overrides(config, target_mode, output_dir)
        if not dry_run:
            print(f'\n  Migrated: {result["migrated"]} files')
            print(f'  Copied:   {result["copied"]} files (precision-independent)')
            if result['skipped']:
                print(f'  Skipped:  {len(result["skipped"])} files')
            _print_divergences(result)
    elif config.language == 'c':
        # BLACS keeps the legacy hardcoded-pattern migrator (it carries
        # MPI typedef patches, Bdef.h rewrites, MPI_REAL16 check
        # generation, and BLACS-specific Cd*/BI_d* routine patterns
        # that have no analogue in other C libraries). Every other C
        # recipe (PBLAS / ScaLAPACK_C / XBLAS / future libraries) uses
        # the generic rename-map-driven cloner — the prefix classifier
        # discovers slot positions empirically, so the cloner is
        # naming-convention-agnostic.
        if config.library == 'blacs':
            result = run_c_migration(config, output_dir, target_mode, dry_run)
        else:
            result = run_c_migration(
                config, output_dir, target_mode, dry_run,
                classification=classification, rename_map=rename_map,
            )
        if not dry_run:
            print(f'\n  Cloned: {len(result["cloned"])} files')
            _print_divergences(result)
    else:
        raise ValueError(f'Unsupported language: {config.language}')

    # Symbol privatization (task 44): recipes that opt in via
    # ``privatize_symbols:`` get every manifest name renamed
    # ``name`` → ``ep_name`` across the full migrated output —
    # strictly after clones, header patches, extern-"C" wrapping and
    # override copies, so nothing regenerates a pristine name behind
    # the pass. Only migrated (extended-precision) stagings ever reach
    # this point; baseline stagings copy sources verbatim.
    if config.privatize_symbols and not dry_run:
        # The header-split step restores each split original (PBblacs.h,
        # pblas.h, …) to its pristine Netlib bytes and points every
        # migrated TU at the ``<pair>``-prefixed sibling; the pass must
        # leave those originals byte-identical and rename only the
        # siblings, which carry the extended stack's public C surface.
        split = result.get('split_headers') if isinstance(result, dict) \
            else None
        pristine_originals = frozenset(
            output_dir / name for name in (split or {}))
        n_privatized = privatize_tree(output_dir, config.privatize_symbols,
                                      skip=pristine_originals)
        print(f'  Privatized: {n_privatized} files rewritten with '
              f'ep_-prefixed symbols '
              f'({len(config.privatize_symbols)} manifest names)')

    return result
