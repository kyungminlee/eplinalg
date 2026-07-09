"""Post-migration symbol privatization (``ep_`` prefix).

The family-independent support archives of the extended-precision stack
(``blacs_common``, ``pblas_common``, ``scalapack_common``,
``ptzblas_common``, plus the ``scalapack_tools`` helpers) are built from
Netlib sources and export the same symbol names as MKL's proprietary
BLACS/PBLAS/ScaLAPACK internals — while disagreeing on internal struct
layout (``BLACBUFF`` is 48 bytes in the Netlib engine, 56 in MKL's).
A consumer linking both gets exactly one definition per name
process-wide, so the two engines collapse into a mixed stack that
corrupts shared state at np>=2 and SIGSEGVs in **either** link order
once the ScaLAPACK type-3 root actually executes.

Fix: the extended stack gets its own hermetic engine. Every name in the
checked-in manifest (``recipes/privatize_ep_symbols.txt``) is renamed
``name -> ep_name`` — definitions, call sites, internal headers and
override files uniformly — in the migrated staged sources of every
recipe that opts in via the ``privatize_symbols:`` key. This is the
source-level equivalent of the ``objcopy --redefine-syms`` rewrite the
task 44 Phase 0 prototype validated end-to-end.

The pass runs at the very end of :func:`pipeline.run_migration`,
strictly after clones, header patches, the extern-"C" wrap and override
copies, so nothing regenerates a pristine name behind it. Baseline
(kind4/kind8) stagings and the verbatim-staged standard-precision trees
never run the migration pipeline at all, so they are untouched by
construction.

Manifest semantics — entries are exact linker-level names:

* **C sources** (``.c``/``.h``) are rewritten with exact,
  case-sensitive whole-token matches. The manifest carries each
  decoration spelling as its own entry (``blacs_gridinit_`` and
  ``Cblacs_gridinit``; the ``reshape`` / ``reshape_`` / ``RESHAPE``
  alias triple), so a plain token substitution covers every convention
  the sources spell out, including the ``#define name_ NAME`` mangling
  blocks. ``#include`` lines are left alone — file names never change,
  so rewriting a token inside one could only break the include.
* **Fortran sources** are rewritten case-insensitively (and
  case-preservingly) on the underscore-stripped stem of every
  single-trailing-underscore entry: ``CALL DESCINIT`` becomes
  ``CALL EP_DESCINIT`` and gfortran emits ``ep_descinit_``, matching
  the renamed definition. Double-underscore aliases (``…__``, g77
  decoration) and C-only lowercase names never enter the Fortran map.
  The ``RESHAPE`` stem shares its spelling with the Fortran intrinsic,
  but it is renamed like any other entry: the ScaLAPACK band solvers
  ``CALL RESHAPE`` (the REDIST C routine — a ``CALL`` can never be the
  intrinsic, which is a function), and neither the staged trees nor
  upstream ScaLAPACK 2.2.3 contain a genuine intrinsic ``RESHAPE(``
  use (every non-CALL match is a ``SL_GRIDRESHAPE(`` substring, which
  whole-token matching never touches).

Struct *members* are never manifest entries (they are not linker
symbols); if one ever sneaks in, the staged-header ABI assertion in
``staging.py`` (``_assert_shared_struct_abi``) aborts the stage loudly.
"""

import re
from pathlib import Path

from .fortran.renames import replace_routine_names

#: Prefix applied to every manifest name.
PRIVATIZE_PREFIX = 'ep_'

_C_IDENT_RE = re.compile(r'[A-Za-z_]\w*')
_C_INCLUDE_RE = re.compile(r'\s*#\s*include\b')
_FORTRAN_CALLABLE_RE = re.compile(r'[A-Za-z]\w*_')

_C_EXTS = frozenset({'.c', '.h'})
_FORTRAN_EXTS = frozenset({'.f', '.for', '.f90', '.f95', '.f03'})


def build_c_map(names) -> dict[str, str]:
    """Exact, case-sensitive C token map: ``name -> ep_name``."""
    return {n: PRIVATIZE_PREFIX + n for n in names}


def build_fortran_map(names) -> dict[str, str]:
    """Case-insensitive Fortran map keyed on upper-cased stems.

    Only single-trailing-underscore manifest entries are
    Fortran-callable; their stem (name minus the decoration underscore)
    is what Fortran source spells. Values are upper-case —
    :func:`replace_routine_names` transfers the token's case.
    """
    fmap: dict[str, str] = {}
    for name in names:
        if not _FORTRAN_CALLABLE_RE.fullmatch(name):
            continue
        stem = name[:-1]
        if stem.endswith('_'):
            continue  # double-underscore decoration alias (g77) — C-only
        stem_u = stem.upper()
        fmap[stem_u] = PRIVATIZE_PREFIX.upper() + stem_u
    return fmap


def privatize_c_text(text: str, c_map: dict[str, str]) -> str:
    """Rename manifest tokens in C source text (case-sensitive).

    ``#include`` lines are passed through verbatim: staged file names
    are never renamed, so touching a token inside an include could only
    break resolution.
    """
    def _sub(m: re.Match) -> str:
        tok = m.group(0)
        return c_map.get(tok, tok)

    out_lines = []
    for line in text.split('\n'):
        if _C_INCLUDE_RE.match(line):
            out_lines.append(line)
        else:
            out_lines.append(_C_IDENT_RE.sub(_sub, line))
    return '\n'.join(out_lines)


def privatize_fortran_text(text: str, f_map: dict[str, str]) -> str:
    """Rename manifest stems in Fortran source text (case-preserving)."""
    return replace_routine_names(text, f_map)


def privatize_tree(root: Path, names, skip=frozenset()) -> int:
    """Apply the privatization pass to every source file under ``root``.

    Dispatches per file extension: C map for ``.c``/``.h``, Fortran map
    for ``.f``/``.F``/``.f90``/… Files whose content is unchanged are
    not rewritten. Returns the number of files modified. Idempotent:
    ``ep_``-prefixed tokens are whole tokens that match no manifest
    entry, so a second run is a no-op.

    ``skip`` is a collection of paths passed through verbatim. It holds
    the restored-pristine split originals of the header-split step
    (``migrate_c_directory``'s ``split_headers``): those files are
    contractually byte-identical to the Netlib reference — the
    precision-prefixed siblings carry the transformed (and privatized)
    surface, and no migrated TU includes the originals.
    """
    skip = {Path(p).resolve() for p in skip}
    c_map = build_c_map(names)
    f_map = build_fortran_map(names)
    changed = 0
    for path in sorted(root.rglob('*')):
        if not path.is_file() or path.resolve() in skip:
            continue
        ext = path.suffix.lower()
        if ext in _C_EXTS:
            original = path.read_text(errors='replace')
            updated = privatize_c_text(original, c_map)
        elif ext in _FORTRAN_EXTS:
            original = path.read_text(errors='replace')
            updated = privatize_fortran_text(original, f_map)
        else:
            continue
        if updated != original:
            path.write_text(updated)
            changed += 1
    return changed
