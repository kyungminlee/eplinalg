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
