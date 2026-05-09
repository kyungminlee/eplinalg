"""Stage upstream sources for a recipe and apply its patch list.

The prepared tree lives at ``<project_root>/build/staged-sources/<library>/``
and mirrors the contents of the recipe's ``source_dir``. Patches listed in
``recipe.patches`` are unified diffs (typically produced by ``git diff
--no-index``) applied with ``git apply`` from inside the staged tree.

Idempotency: a ``.prepared.stamp`` file inside the staged tree records when
preparation last completed. If the stamp is newer than every listed patch
file, the stage is skipped. Pass ``rebuild=True`` to force a clean re-stage.

This module is the input side of the pipeline reshape described in
``doc/refactor-20260509.md``. Phase A wires it as a no-op-by-default CLI
command; subsequent phases route migration to read from the staged tree.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from .config import RecipeConfig, load_recipe


STAMP_NAME = '.prepared.stamp'


def staged_root_for(project_root: Path, library: str) -> Path:
    """Where ``library``'s staged sources live."""
    return project_root / 'build' / 'staged-sources' / library


def patch_dir_for(recipe_path: Path, library: str) -> Path:
    """Directory holding the recipe's patch files."""
    return recipe_path.parent / library / 'patches'


def _resolve_patches(recipe_path: Path,
                     config: RecipeConfig) -> list[Path]:
    pdir = patch_dir_for(recipe_path, config.library)
    resolved: list[Path] = []
    missing: list[str] = []
    for name in config.patches:
        p = pdir / name
        if not p.exists():
            missing.append(name)
        resolved.append(p)
    if missing:
        raise FileNotFoundError(
            f'{recipe_path.name}: patch file(s) not found in {pdir}: '
            f'{missing!r}'
        )
    return resolved


def _is_fresh(stamp: Path, patches: list[Path]) -> bool:
    if not stamp.exists():
        return False
    stamp_mtime = stamp.stat().st_mtime
    for p in patches:
        if p.stat().st_mtime > stamp_mtime:
            return False
    return True


def run_prepare(recipe_path: Path,
                project_root: Path | None = None,
                *,
                rebuild: bool = False) -> Path:
    """Prepare the staged-sources tree for one recipe.

    Returns the staged root path.
    """
    if project_root is None:
        project_root = recipe_path.parent.parent

    config = load_recipe(recipe_path, project_root)
    staged_root = staged_root_for(project_root, config.library)
    stamp = staged_root / STAMP_NAME

    patches = _resolve_patches(recipe_path, config)

    if not rebuild and staged_root.exists() and _is_fresh(stamp, patches):
        return staged_root

    if staged_root.exists():
        shutil.rmtree(staged_root)
    staged_root.parent.mkdir(parents=True, exist_ok=True)

    if not config.source_dir.exists():
        raise FileNotFoundError(
            f'{recipe_path.name}: source_dir does not exist: '
            f'{config.source_dir}'
        )
    shutil.copytree(config.source_dir, staged_root)

    for patch in patches:
        # ``git apply`` is run from inside the staged tree so a patch
        # whose hunks reference top-level filenames (``a/dnrm2.f90``)
        # resolves directly. The patch path is made absolute first
        # because git apply re-resolves it against ``cwd=staged_root``,
        # not the caller's CWD.
        patch_abs = str(patch.resolve())
        # ``--check`` first: refusing on apply-conflict surfaces upstream
        # drift loud (vendor bumped a file the patch covers) instead of
        # silently producing a half-applied tree.
        subprocess.run(
            ['git', 'apply', '--whitespace=nowarn', '--check', patch_abs],
            cwd=staged_root, check=True,
        )
        subprocess.run(
            ['git', 'apply', '--whitespace=nowarn', patch_abs],
            cwd=staged_root, check=True,
        )

    stamp.touch()
    return staged_root


def _precision_sibling_set(filename: str,
                           upstream_root: Path | None = None) -> set[str] | None:
    """Return the co-family for ``filename`` (existing upstream files only).

    Generates candidate sibling names by precision-prefix swap, with the
    LAPACK orthogonal/unitary aliasing (``or`` ↔ ``un``) folded in:

      ``dorbdb3.f``  →  {``sorbdb3.f``, ``dorbdb3.f``, ``cunbdb3.f``, ``zunbdb3.f``}
      ``pdormrz.f``  →  {``psormrz.f``, ``pdormrz.f``, ``pcunmrz.f``, ``pzunmrz.f``}
      ``pdgeequ.f``  →  {``psgeequ.f``, ``pdgeequ.f``, ``pcgeequ.f``, ``pzgeequ.f``}

    When ``upstream_root`` is provided, candidates that don't correspond
    to a real upstream file are filtered out — so a candidate like
    ``dunbdb3.f`` (which doesn't exist; D-half is ``dorbdb3.f``) won't
    show up as a missing sibling.

    Returns ``None`` for files whose stem doesn't begin with a known
    precision prefix (e.g. ``utils.f``).
    """
    p = Path(filename)
    stem = p.stem.lower()
    suf = p.suffix

    # Detect prefix style. Try 2-letter (PS/PD/PC/PZ) first.
    for prefix_set, pfx_len in ((('ps', 'pd', 'pc', 'pz'), 2),
                                 (('s',  'd',  'c',  'z'),  1)):
        pfx = stem[:pfx_len]
        if pfx not in prefix_set or len(stem) <= pfx_len:
            continue
        tail = stem[pfx_len:]
        precision = pfx[-1]  # 's','d','c','z'
        is_real = precision in ('s', 'd')

        candidates: set[str] = set()
        for q in prefix_set:
            q_precision = q[-1]
            q_is_real = q_precision in ('s', 'd')
            # Apply or/un aliasing when crossing the real/complex boundary.
            new_tail = tail
            if is_real and not q_is_real and tail.startswith('or'):
                new_tail = 'un' + tail[2:]
            elif not is_real and q_is_real and tail.startswith('un'):
                new_tail = 'or' + tail[2:]
            candidates.add(f'{q}{new_tail}{suf}')

        if upstream_root is not None and upstream_root.is_dir():
            existing = {n for n in candidates
                        if (upstream_root / n).is_file()}
            return existing if existing else None
        return candidates
    return None


def _patch_touched_files(patch_path: Path) -> set[str]:
    """Return the set of filenames the patch's hunks add or modify."""
    touched: set[str] = set()
    for line in patch_path.read_text().splitlines():
        if line.startswith('+++ b/'):
            rel = line[len('+++ b/'):].strip().split('\t', 1)[0]
            touched.add(Path(rel).name)
    return touched


def verify_patches(recipe_path: Path,
                   project_root: Path | None = None) -> list[str]:
    """Symmetric-patch CI check.

    Aggregates every file touched by any patch in the recipe (excluding
    those whose patch is listed under ``asymmetric_patches:``). For
    each precision-prefixed file in the aggregate, requires all four
    co-family siblings to be touched by SOME patch — otherwise a fix
    in the D/Z half doesn't propagate to the S/C re-migration path.

    Sibling hunks may live in different patches; they often do
    (per-stem patches keep diffs reviewable).

    Returns a list of error strings; empty list means clean.
    """
    if project_root is None:
        project_root = recipe_path.parent.parent
    config = load_recipe(recipe_path, project_root)
    pdir = patch_dir_for(recipe_path, config.library)
    if not pdir.is_dir():
        return []

    skip = set(config.asymmetric_patches) | set(config.one_sided_cleanup)
    aggregate: set[str] = set()
    file_to_patches: dict[str, list[str]] = {}
    errors: list[str] = []

    for patch_name in config.patches:
        if patch_name in skip:
            continue
        patch_path = pdir / patch_name
        if not patch_path.exists():
            errors.append(f'{patch_name}: not found in {pdir}')
            continue
        for f in _patch_touched_files(patch_path):
            aggregate.add(f)
            file_to_patches.setdefault(f, []).append(patch_name)

    upstream_root = config.source_dir
    for f in sorted(aggregate):
        sibs = _precision_sibling_set(f, upstream_root)
        if sibs is None:
            continue
        missing = sibs - aggregate
        if missing:
            errors.append(
                f'{f} (in {file_to_patches[f]}): precision-asymmetric — '
                f'co-family siblings {sorted(missing)} are not touched '
                f'by any patch. Add a patch that touches them, or list '
                f'the patch under asymmetric_patches: in the recipe.'
            )

    return errors


def prepare_recipe(recipe_path: Path,
                   project_root: Path | None = None,
                   *,
                   rebuild: bool = False) -> RecipeConfig:
    """Load a recipe with ``source_dir`` rewritten to its staged tree.

    Runs :func:`run_prepare` to ensure ``build/staged-sources/<library>/``
    exists and reflects the recipe's patch list, then loads the recipe
    and swaps ``config.source_dir`` to point at the staged tree. All
    downstream pipeline code (`scan_symbols`, `run_fortran_migration`,
    `run_c_migration`, …) reads through ``config.source_dir``, so this
    one swap routes the entire migration through the staged tree.
    """
    if project_root is None:
        project_root = recipe_path.parent.parent
    staged_root = run_prepare(recipe_path, project_root, rebuild=rebuild)
    config = load_recipe(recipe_path, project_root)
    config.source_dir = staged_root
    return config
