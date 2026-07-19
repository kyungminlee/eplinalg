"""Shared CLI helpers used across the migrator subcommands.

Small argument-to-object adapters (target-mode construction, parser
selection) and path/name conventions that several ``cmd_*`` handlers and
``staging`` depend on. Kept in their own module so both ``__main__`` and
``staging`` can import them without a circular dependency. Extracted
verbatim from ``__main__``.
"""
from pathlib import Path

from .target_mode import load_target

# Every non-empty ``la_constants_suffix`` across codegen/targets/*.yaml
# (kind10 → _ey, kind16 → _qx, multifloats → _mw). This is the universe
# of per-target LA_CONSTANTS / LA_XISNAN module suffixes used to route
# the helper pair between OWN/FOREIGN sets during builds and staging.
# The staged top-level CMakeLists.txt (cmake/CMakeLists.txt, LAPACK
# EXCLUDE_REGEX) hardcodes the same universe — keep them in sync.
LA_SUFFIXES = ('_ey', '_qx', '_mw')


def la_helper_pairs(target_mode) -> tuple[set[str], set[str]]:
    """Return ``(own, foreign)`` LA_CONSTANTS / LA_XISNAN stem sets.

    ``own`` holds the active target's pair-suffixed helper stems;
    ``foreign`` holds every other target's, i.e. the stems the current
    build must not compile from the shared LAPACK SRC directory.
    """
    own_suffix = target_mode.la_constants_suffix.lower()
    own = {f'la_constants{own_suffix}', f'la_xisnan{own_suffix}'}
    foreign = {
        f'la_{base}{sfx}'
        for base in ('constants', 'xisnan')
        for sfx in LA_SUFFIXES
        if sfx != own_suffix
    }
    return own, foreign


def recipe_project_root(recipe_path: Path) -> Path:
    """Default project root inferred from a recipe path.

    Recipes live at ``<root>/codegen/recipes/<name>.yaml``, three
    levels below the repo root. The path is used as given — callers
    that want symlink resolution resolve before calling.
    """
    return recipe_path.parent.parent.parent


def migrator_project_root() -> Path:
    """Project root located from this package's position on disk
    (``<root>/codegen/migrator/`` → repo root)."""
    return Path(__file__).resolve().parent.parent.parent


def get_target_mode(args):
    """Construct TargetMode based on CLI arguments."""
    target_str = getattr(args, 'target', None) or 'kind16'
    return load_target(target_str)


def parser_args(args):
    """Extract parser/parser_cmd from CLI args."""
    parser = getattr(args, 'parser', None)
    parser_cmd = getattr(args, 'parser_cmd', None)
    return parser, parser_cmd
