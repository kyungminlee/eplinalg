"""Shared clone-and-substitute primitives for C sources.

Type-specific C files are near-identical clones differing only in C
type names and MPI datatype constants, so cloning is mechanical text
substitution — no clang parser needed.  Both the on-disk directory
migrations (BLACS, generic/ScaLAPACK) and the in-memory single-file
harness apply these same tables and transforms, so the tested path and
the shipped path cannot diverge.
"""

import re
from pathlib import Path

from ..templates import expand_template


# Default C type substitution rules.
# Each rule is (pattern, replacement_template).
# In replacement_template:
#   {REAL_TYPE}    → the target real type name (e.g., "QREAL", "EREAL")
#   {COMPLEX_TYPE} → the target complex type name (e.g., "QCOMPLEX")
#   {MPI_REAL}     → the target MPI real type (e.g., "MPI_QREAL")
#   {MPI_COMPLEX}  → the target MPI complex type (e.g., "MPI_QCOMPLEX")
#   {RP}           → real prefix char lowercase (e.g., "q", "e")
#   {CP}           → complex prefix char lowercase (e.g., "x", "y")
#   {RPU}          → real prefix char uppercase
#   {CPU}          → complex prefix char uppercase
#   ...
#   {C_REAL_TYPE}  → underlying C type (e.g., "__float128")

# Clone-sub builders. Each table is a flat list of (pattern, replacement)
# pairs applied in order to one C source line at a time. ``\b`` boundary
# matches catch adjacent tokens in a single pass (so no "run twice" dup
# is needed) and exclude underscore-adjacent identifiers like
# ``float64x2_t``. Prefix rules (``Cd*``/``BI_d*``) use ``(?<![A-Za-z0-9_])``
# lookbehind rather than ``\b`` because the leading prefix letter is
# lowercase and may follow a word boundary that ``\b`` would mistake.
#
# {MPI_SUM_REAL}/{MPI_SUM_COMPLEX} are textual no-ops for KIND targets
# (both expand to ``MPI_SUM``); multifloats expands them to its
# user-defined ops (``MPI_MM_SUM`` / the zz variant) registered by libmfc.


def _real_clone_subs(real_keyword: str, src_letter: str) -> list[tuple[str, str]]:
    mpi_real_src = f'MPI_{real_keyword.upper()}'
    return [
        (rf'\b{real_keyword}\b', '{REAL_TYPE}'),
        (rf'\b{mpi_real_src}\b', '{MPI_REAL}'),
        (r'\bMPI_SUM\b', '{MPI_SUM_REAL}'),
        # Function name prefixes (allow uppercase suffix for BI_dMPI_* etc.).
        (rf'(?<![A-Za-z0-9_])C{src_letter}([a-z])', r'C{RP}\1'),
        (rf'(?<![A-Za-z0-9_])BI_{src_letter}([a-zA-Z])', r'BI_{RP}\1'),
    ]


def _complex_clone_subs(real_keyword: str, own_struct: str, cross_struct: str,
                        complex_letter: str, real_letter: str) -> list[tuple[str, str]]:
    mpi_real_src = f'MPI_{real_keyword.upper()}'
    mpi_complex_src = f'MPI_{real_keyword.upper()}_COMPLEX'
    return [
        (rf'\b{own_struct}\b', '{COMPLEX_TYPE}'),
        (rf'\b{cross_struct}\b', '{COMPLEX_TYPE}'),
        (rf'\b{real_keyword}\b', '{REAL_TYPE}'),
        # MPI types — order matters: <PREC>_COMPLEX before <PREC>.
        (rf'\b{mpi_complex_src}\b', '{MPI_COMPLEX}'),
        (rf'\b{mpi_real_src}\b', '{MPI_REAL}'),
        # Gated MPI_COMPLEX rule (preserved as-is: no leading boundary).
        (r'MPI_COMPLEX([^a-zA-Z_0-9])', r'{MPI_COMPLEX}\1'),
        (r'\bMPI_SUM\b', '{MPI_SUM_COMPLEX}'),
        (rf'(?<![A-Za-z0-9_])C{complex_letter}([a-z])', r'C{CP}\1'),
        (rf'(?<![A-Za-z0-9_])BI_{complex_letter}([a-zA-Z])', r'BI_{CP}\1'),
        (rf'(?<![A-Za-z0-9_])BI_{real_letter}([a-zA-Z])', r'BI_{RP}\1'),
    ]


REAL_CLONE_SUBS = _real_clone_subs('double', 'd')
COMPLEX_CLONE_SUBS = _complex_clone_subs('double', 'DCOMPLEX', 'SCOMPLEX', 'z', 'd')
# S/C convergence-only mirror tables: re-derive the Q/X target from the
# S/C half (sources use ``float``/``MPI_FLOAT``/``Cs*``/``Cc*``) so the
# in-memory result can be compared against the canonical produced from
# the D/Z half on disk.
SINGLE_CLONE_SUBS = _real_clone_subs('float', 's')
CSINGLE_CLONE_SUBS = _complex_clone_subs('float', 'SCOMPLEX', 'DCOMPLEX', 'c', 's')


def classify_blacs_stem(stem: str) -> tuple[str, bool, list[tuple[str, str]]] | None:
    """Identify the precision variant of a BLACS-style source stem.

    Returns ``(src_prefix, is_complex, subs)`` where ``src_prefix`` is
    the source precision letter (``d``/``z``/``s``/``c``), or ``None``
    for precision-independent files.  Order matters: check
    BI_-prefixed names before bare single-letter prefixes.
    """
    if stem.startswith('BI_d'):
        return 'd', False, REAL_CLONE_SUBS
    if stem.startswith('BI_z'):
        return 'z', True, COMPLEX_CLONE_SUBS
    if stem.startswith('BI_s'):
        return 's', False, SINGLE_CLONE_SUBS
    if stem.startswith('BI_c'):
        return 'c', True, CSINGLE_CLONE_SUBS
    if stem.startswith('BI_'):
        return None  # precision-independent BI_* helper
    if stem.startswith('d'):
        return 'd', False, REAL_CLONE_SUBS
    if stem.startswith('z'):
        return 'z', True, COMPLEX_CLONE_SUBS
    if stem.startswith('s'):
        return 's', False, SINGLE_CLONE_SUBS
    if stem.startswith('c'):
        return 'c', True, CSINGLE_CLONE_SUBS
    return None


def apply_clone_transform(text: str,
                          subs: list[tuple[str, str]],
                          template_vars: dict[str, str],
                          routine_renames: list[tuple[str, str]] | None = None,
                          ) -> str:
    """Apply clone substitutions plus routine renames to C source text.

    routine_renames is a list of (old_name, new_name) pairs for literal
    routine name replacements (both lowercase and uppercase are applied).

    The renames use a left-side word-boundary negative lookbehind so we
    don't double-rename inside identifiers already produced by the regex
    prefix substitutions above. For example with multifloats, the regex
    pass turns ``Cdgesd2d`` into ``Cddgesd2d``; a plain
    text.replace('dgesd2d', 'ddgesd2d') would then find the substring
    inside ``Cddgesd2d`` and produce ``Cdddgesd2d``. The lookbehind
    prevents that because the 'd' at offset 1 is preceded by another
    word character ('C' / 'd'). Single-char-prefix KIND targets dodge
    this because their regex produces names like ``Cqgesd2d`` that no
    longer contain ``dgesd2d`` as a substring.
    """
    for pattern, replacement in subs:
        expanded = expand_template(replacement, template_vars)
        text = re.sub(pattern, expanded, text, flags=re.MULTILINE)

    if routine_renames:
        for old_name, new_name in routine_renames:
            text = re.sub(
                r'(?<![A-Za-z0-9_])' + re.escape(old_name),
                new_name, text)
            text = re.sub(
                r'(?<![A-Za-z0-9_])' + re.escape(old_name.upper()),
                new_name.upper(), text)
    return text


def clone_c_file(src_path: Path, dst_path: Path,
                 subs: list[tuple[str, str]],
                 template_vars: dict[str, str],
                 routine_renames: list[tuple[str, str]] | None = None,
                 ) -> None:
    """Clone a C file with mechanical text substitutions."""
    text = src_path.read_text(errors='replace')
    text = apply_clone_transform(text, subs, template_vars, routine_renames)
    dst_path.write_text(text)


def derive_routine_renames(old_stem: str, new_stem: str) -> list[tuple[str, str]]:
    """Derive routine name renames from source/target file stems.

    For user-facing files like 'dgesd2d_' → 'qgesd2d_', strips the
    trailing underscore to get the routine base name and returns rename
    pairs.  For BI_-prefixed files the function names are already handled
    by the regex rules, so this returns an empty list.
    """
    if old_stem.startswith('BI_'):
        return []
    # Strip trailing underscore if present (Fortran naming convention)
    old_routine = old_stem.rstrip('_')
    new_routine = new_stem.rstrip('_')
    if old_routine == new_routine:
        return []
    return [(old_routine, new_routine)]


def rename_c_file(name: str, old_prefix: str, new_prefix: str) -> str:
    """Rename a C file by replacing its prefix character."""
    stem = Path(name).stem
    ext = Path(name).suffix

    if stem.startswith(f'BI_{old_prefix}'):
        new_stem = f'BI_{new_prefix}' + stem[len(f'BI_{old_prefix}'):]
    elif stem.startswith(old_prefix):
        new_stem = new_prefix + stem[len(old_prefix):]
    else:
        return name

    return new_stem + ext
