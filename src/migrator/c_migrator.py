"""General-purpose C type migration engine.

Handles C source files that use the template-clone pattern (like BLACS):
type-specific files are near-identical clones differing only in C type
names and MPI datatype constants.

Migration is done by cloning files with mechanical text substitution.
No clang parser needed — C types are unambiguous single tokens.
"""

import re
import shutil
from pathlib import Path

from .prefix_classifier import target_prefix, SymbolClassification
from .target_mode import TargetMode
from .pbcharshim import (
    is_typeset_stem,
    transform_typeset,
    pbcharshim_header_text,
)


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


def _build_sub_vars(target_mode: TargetMode) -> dict[str, str]:
    """Build template substitution variables for a given target mode.

    All values are read from the target's c_interop YAML section
    (populated in TargetMode by ``load_target()``).
    """
    assert target_mode.c_real_type is not None, (
        "target_mode missing c_real_type; ensure the target YAML has a "
        "c_interop section with real_type defined."
    )
    rp = target_prefix(target_mode, is_complex=False).lower()
    cp = target_prefix(target_mode, is_complex=True).lower()
    return {
        'REAL_TYPE': target_mode.c_real_type,
        'COMPLEX_TYPE': target_mode.c_complex_type,
        'C_REAL_TYPE': target_mode.c_c_real_type,
        'MPI_REAL': target_mode.c_mpi_real,
        'MPI_COMPLEX': target_mode.c_mpi_complex,
        'MPI_SUM_REAL': target_mode.c_mpi_sum_real,
        'MPI_SUM_COMPLEX': target_mode.c_mpi_sum_complex,
        'RP': rp,
        'CP': cp,
        'RPU': rp.upper(),
        'CPU': cp.upper(),
    }


# ------------------------------------------------------------------ #
# Source-stem decoration helpers                                     #
# ------------------------------------------------------------------ #

# XBLAS f2c-bridge file naming: ``BLAS_dgemv_x-f2c.c`` carries the
# Fortran-callable wrapper for the C routine ``BLAS_dgemv_x``. The
# bridge files share the routine's stem with a ``-f2c`` decoration,
# so they need the same precision-prefix rewrite as their parent
# routine. The migrator keys on routine stems though, so we strip
# this decoration before the rename-map lookup and re-append it on
# the way out.
_C_FILE_DECORATIONS: tuple[str, ...] = ('-f2c',)


def _strip_decoration(stem: str) -> tuple[str, str]:
    """Split ``BLAS_dgemv_x-f2c`` into (``BLAS_dgemv_x``, ``-f2c``).

    Returns (base_stem, decoration). Returns (stem, '') if no known
    decoration is present.
    """
    for deco in _C_FILE_DECORATIONS:
        if stem.endswith(deco):
            return stem[: -len(deco)], deco
    return stem, ''


def _redist_clone_stem(routine: str,
                       target_mode: TargetMode | None) -> str | None:
    """Clone-stem fallback for REDIST/SRC-style C files whose stems are
    ``p<sdcz><root>`` but whose Fortran-callable exports live behind
    ``#define`` macros (so the symbol scanner never registers them).

    Returns the precision-substituted lowercase stem (``pdgemr2`` →
    ``pqgemr2`` for kind16, ``pcgemr`` → ``pxgemr``), or ``None`` if
    the routine doesn't match the expected ``p<precision-letter>...``
    shape — which signals to skip the file.
    """
    if target_mode is None or len(routine) < 2:
        return None
    upper = routine.upper()
    if upper[0] != 'P' or upper[1] not in ('S', 'D', 'C', 'Z'):
        return None
    from .prefix_classifier import CHAR_TYPE
    family = CHAR_TYPE[upper[1]]
    new_char = target_mode.prefix_map.get(family)
    if new_char is None:
        return None
    return 'p' + new_char.lower() + routine[2:].lower()




def clone_c_file(src_path: Path, dst_path: Path,
                 subs: list[tuple[str, str]],
                 template_vars: dict[str, str],
                 routine_renames: list[tuple[str, str]] | None = None,
                 ) -> None:
    """Clone a C file with mechanical text substitutions.

    routine_renames is a list of (old_name, new_name) pairs for literal
    routine name replacements (both lowercase and uppercase are applied).
    """
    text = src_path.read_text(errors='replace')

    for pattern, replacement in subs:
        # Expand template variables in replacement
        expanded = replacement
        for key, val in template_vars.items():
            expanded = expanded.replace(f'{{{key}}}', val)
        text = re.sub(pattern, expanded, text, flags=re.MULTILINE)

    # Apply routine name renames. Use a left-side word-boundary
    # negative lookbehind so we don't double-rename inside identifiers
    # already produced by the regex prefix substitutions above. For
    # example with multifloats, the regex pass turns ``Cdgesd2d`` into
    # ``Cddgesd2d``; a plain text.replace('dgesd2d', 'ddgesd2d') would
    # then find the substring inside ``Cddgesd2d`` and produce
    # ``Cdddgesd2d``. The lookbehind prevents that because the 'd' at
    # offset 1 is preceded by another word character ('C' / 'd').
    # Single-char-prefix KIND targets dodge this because their regex
    # produces names like ``Cqgesd2d`` that no longer contain
    # ``dgesd2d`` as a substring.
    if routine_renames:
        for old_name, new_name in routine_renames:
            text = re.sub(
                r'(?<![A-Za-z0-9_])' + re.escape(old_name),
                new_name, text)
            text = re.sub(
                r'(?<![A-Za-z0-9_])' + re.escape(old_name.upper()),
                new_name.upper(), text)

    dst_path.write_text(text)


def _routine_renames(old_stem: str, new_stem: str) -> list[tuple[str, str]]:
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


def migrate_c_directory(src_dir: Path, output_dir: Path,
                        target_mode: TargetMode,
                        classification: SymbolClassification | None = None,
                        rename_map: dict[str, str] | None = None,
                        c_type_aliases: list[dict] | None = None,
                        c_pointer_cast_aliases: list[dict] | None = None,
                        header_patches: list[dict] | None = None,
                        overrides: list[tuple[Path, str]] | None = None,
                        extra_c_dirs: list[Path] | None = None,
                        skip_files: set[str] | None = None,
                        copy_files: set[str] | None = None) -> dict:
    """Migrate a C source directory by cloning real/complex variants.

    Two modes:

    - **Generic** (when `classification` and `rename_map` are supplied):
      clones each D/Z member of every precision family found in
      `classification`, applying `rename_map` for all in-file routine
      name substitutions plus C-type upgrades (``double`` →
      ``{REAL_TYPE}``). Used for ScaLAPACK-style libraries like PBLAS
      whose entry points are Fortran-callable names such as
      ``pdgemm_``, ``pzhemm_``, ``pdznrm2_``.

    - **BLACS** (when neither is supplied): hardcoded ``d → q`` / ``z →
      x`` file renames with BLACS-specific ``Cd*``/``BI_d*`` routine
      patterns, MPI type substitutions, ``Bdef.h`` patching, and an
      MPI_REAL16 check module for KIND=16. Preserved for backward
      compatibility.

    ``overrides`` is a list of (src_path, dst_name) pairs that are
    copied verbatim on top of clones after the main migration step
    and after header patches have run. Used for hand-written
    replacement kernels that cannot be produced by regex substitution.

    Returns a summary dict.
    """
    if classification is not None and rename_map is not None:
        result = _migrate_generic_c_directory(
            src_dir, output_dir, target_mode,
            classification, rename_map,
            c_type_aliases=c_type_aliases,
            c_pointer_cast_aliases=c_pointer_cast_aliases,
            header_patches=header_patches,
            extra_c_dirs=extra_c_dirs,
            skip_files=skip_files,
            copy_files=copy_files,
        )
    else:
        result = _migrate_blacs_c_directory(
            src_dir, output_dir, target_mode,
            skip_files=skip_files,
        )
    if overrides:
        applied = _apply_overrides(output_dir, overrides)
        result['overrides'] = applied
    # Multifloats: ``.c`` sources are compiled as C++ for operator
    # overloading on float64x2_t. Wrap each file body in ``extern "C"``
    # *after* the last ``#include`` so the Fortran-callable /
    # C-callable entry points (blacs_*_, Cblacs_*, BI_*, pddgemm_,
    # pzzherk_, …) keep C linkage. Bridge headers'
    # C++ templates and mpicxx.h declarations stay in C++ scope above
    # the wrap. Bdef.h gets the same after-includes wrap so its
    # forward declarations agree with the wrapped definitions. Run
    # this *after* ``_apply_overrides`` so hand-written replacements
    # are wrapped too.
    split_headers = result.get('split_headers') if isinstance(result, dict) \
        else None
    if target_mode is not None and not target_mode.is_kind_based:
        _wrap_extern_c_after_last_include(output_dir, split_headers)
    # Point every migrated translation unit at the precision-prefixed
    # sibling headers produced by the split step (see
    # ``_migrate_generic_c_directory``). Runs last so it also rewrites
    # hand-written override sources and the prefixed siblings' own
    # cross-includes; pristine restored originals are left untouched.
    if split_headers:
        targets = list(output_dir.glob('*.c'))
        targets += [output_dir / s for s in split_headers.values()]
        for t in targets:
            _rewrite_split_includes(t, split_headers)
    return result


def _wrap_extern_c_after_last_include(
        output_dir: Path,
        split_headers: dict[str, str] | None = None) -> None:
    """Give every Fortran-callable function definition C linkage when
    the source is compiled as C++.

    Two passes per file:

    1. **Wrap the body after the last #include** with
       ``#ifdef __cplusplus extern "C" { #endif`` … closing brace at
       EOF. This covers function definitions while keeping the
       includes themselves outside the wrap — necessary because some
       headers (Bdef.h, multifloats_bridge.h, redist.h via the bridge)
       transitively pull in C++ stdlib templates that cannot live
       inside ``extern "C"``.

    2. **Wrap each contiguous block of forward ``extern <type> Foo(…);``
       declarations** in its own ``extern "C"`` block. These appear
       between #include groups in scalapack_c REDIST sources
       (pcgemr.c, pctrmr.c, …) and would otherwise mismatch the C
       linkage of the matching definitions below the cut.

    Idempotent.
    """
    # Guard every migrated *translation unit* uniformly — no hardcoded
    # per-file special-casing. The target set is driven by the migration
    # mechanism:
    #
    #   * every ``.c`` body (Fortran-callable / C-callable entry points
    #     defined there);
    #   * every precision-prefixed *sibling* header produced by the split
    #     step (``split_headers.values()``) — these carry the transformed,
    #     non-pristine content, including template-instantiated definitions
    #     such as ``mwlamov.h``'s ``LAMOV``/``LACPY``;
    #   * ``Bdef.h`` when present — the BLACS path patches it in place
    #     (it is not installed and not split), so its forward declarations
    #     must agree with the wrapped definitions.
    #
    # Deliberately NOT wrapped: restored-pristine split originals and any
    # unchanged Netlib header (e.g. ``tools.h``, which is pure macros and
    # typedefs with no C-linkage symbol). Their bytes must stay identical
    # to the upstream reference, and they need no guard. Because the set is
    # exactly {.c} ∪ {siblings} ∪ {Bdef.h}, a byte-pristine public header is
    # never a target — the pristineness constraint holds by construction,
    # not by a short-circuit. (Headers that already carry their own
    # ``extern "C"`` are still skipped by the check below.)
    targets = list(output_dir.glob('*.c'))
    targets += [output_dir / s for s in (split_headers or {}).values()]
    bdef = output_dir / 'Bdef.h'
    if bdef.exists():
        targets.append(bdef)
    open_block = '#ifdef __cplusplus\nextern "C" {\n#endif\n'
    close_block = '#ifdef __cplusplus\n} /* extern "C" */\n#endif\n'
    for f in sorted(targets):
        text = f.read_text(errors='replace')
        if 'extern "C"' in text:
            continue
        includes = list(re.finditer(r'(?m)^\s*#\s*include\b[^\n]*\n', text))
        if not includes:
            continue
        cut = includes[-1].end()
        body = text[cut:]
        prefix = text[:cut]
        # Pass 2: wrap contiguous forward extern declaration blocks in
        # the prefix region so they share C linkage with definitions.
        # Match either ``extern <type> Name(...);`` or a bare
        # ``<type> Name(...);`` single-line forward declaration. Group
        # consecutive such lines into a single wrap. ``[^;{}\n]*``
        # forbids the parameter list from spanning lines so the regex
        # cannot accidentally swallow Fortran-style ``SUBROUTINE Foo(`
        # text living inside a multi-line C comment.
        decl_one = (r'^[ \t]*(?:extern[ \t]+)?'
                    r'[A-Za-z_][\w\s\*]*?'
                    r'[A-Za-z_]\w*\s*\([^;{}\n]*\)\s*;[ \t]*\n')
        ext_re = re.compile(r'(?m)(' + decl_one + r'(?:' + decl_one + r')*)')
        prefix = ext_re.sub(lambda m: open_block + m.group(1) + close_block,
                            prefix)
        text = prefix + '\n' + open_block + body + '\n' + close_block
        f.write_text(text)


def _apply_overrides(output_dir: Path,
                     overrides: list[tuple[Path, str]]) -> list[str]:
    """Copy hand-written replacement files on top of migrated clones.

    Each entry is ``(src_path, dst_name)``. Any pre-existing
    ``<output_dir>/<dst_name>`` is overwritten. Missing sources raise
    FileNotFoundError so recipe typos are caught early.
    """
    applied: list[str] = []
    for src_path, dst_name in overrides:
        if not src_path.is_file():
            raise FileNotFoundError(
                f"override file not found: {src_path}"
            )
        dst = output_dir / dst_name
        shutil.copy2(src_path, dst)
        applied.append(dst_name)
    return applied


def _rewrite_split_includes(path: Path,
                            split_headers: dict[str, str]) -> None:
    """Rewrite ``#include "name.h"`` → ``#include "<prefix>name.h"``.

    ``split_headers`` maps an original Netlib header name to the
    precision-prefixed sibling that carries this precision's transformed
    content. Migrated sources (and prefixed siblings that cross-include
    another split header) must reference the sibling, never the restored-
    pristine original. Both quote and angle-bracket include forms are
    rewritten (the delimiter is preserved). Idempotent: already-rewritten
    includes no longer match the original name.
    """
    text = path.read_text(errors='replace')
    new = text
    for orig, sibling in split_headers.items():
        esc = re.escape(orig)
        new = re.sub(r'(#\s*include\s*")' + esc + r'(")',
                     r'\1' + sibling + r'\2', new)
        new = re.sub(r'(#\s*include\s*<)' + esc + r'(>)',
                     r'\1' + sibling + r'\2', new)
    if new != text:
        path.write_text(new)


# C identifier tokenizer used by the rename substituter. Cheap DFA
# that the regex engine evaluates in one linear pass — replaces the
# previous 1500+ alternation regex whose backtracking dominated C
# migration time (90%+ of xblas runtime).
_C_IDENT_RE = re.compile(r'[A-Za-z_]\w*')


def _build_rename_regex(rename_map: dict[str, str]) -> tuple[re.Pattern, dict[str, str]]:
    """Build a tokenizer + lookup map for renaming routine names in C text.

    Strategy: scan the text for identifier-shaped tokens and dict-look
    each one up. Returns ``(_C_IDENT_RE, combined)`` where ``combined``
    maps lowercase token → lowercase replacement (with and without the
    Fortran-style trailing underscore that the BLACS C bridge appends
    to call sites).

    The pre-cleanup form built one giant alternation regex with a
    ``(?<![.>])\\b...\\b`` boundary; the substituter (see
    :func:`_make_rename_substituter`) now enforces the same struct-
    member-access guard by inspecting the character preceding the
    matched token via ``Match.string`` / ``Match.start``, preserving
    PBLAS PBTYP_T compatibility where field names like
    ``TypeStruct.Cgesd2d`` must not be renamed.

    Matching is case-insensitive: the substituter normalizes the
    matched token to lowercase before looking it up, then transfers
    each source character's case onto the replacement.
    """
    combined: dict[str, str] = {}
    for old, new in rename_map.items():
        combined[old.lower()] = new.lower()
        combined[old.lower() + '_'] = new.lower() + '_'
    return _C_IDENT_RE, combined


# PBTYP_T C-callable BLACS binding pointers (pblas.h). Their Titlecase
# spellings collide *case-insensitively* with the complex-single BLACS
# Fortran routines ``cgesd2d_`` / ``CGESD2D`` (the leading ``C`` reads as
# precision letter ``c``), so the case-insensitive rename would rewrite
# ``Cgesd2d`` → ``Wgesd2d`` and — via the header-duplication pass — inject
# a dead ``W*`` twin field into the PBTYP_T struct *definition*. That
# grows the struct by 8 bytes per field, shifting every pointer after it
# and desyncing the typed pblas/src/pblas.h layout from the reference
# _pblas_src/pblas.h (and from MKL's stock ScaLAPACK PBTYP_T), so the
# type-generic PTOOLS (PB_CpaxpbyNN, …) read ``TYPE->Cgerv2d`` as NULL and
# SIGSEGV. The struct-member *access* form (``TYPE->Cgesd2d``) is already
# protected by the ``.``/``>`` lookbehind; this guards the bare member
# *declaration*. These exact Titlecase spellings never appear as a BLACS
# call site (those are all-lower ``cgesd2d_`` or all-upper ``CGESD2D``),
# so protecting them case-sensitively leaves the genuine routine renames
# untouched.
_PBTYP_C_BINDING_FIELDS = frozenset({
    'Cgesd2d', 'Cgerv2d', 'Cgebs2d', 'Cgebr2d', 'Cgsum2d',
})


def _make_rename_substituter(pattern: re.Pattern, combined: dict[str, str]):
    """Return a callback for ``pattern.sub`` that renames case-preservingly.

    Renames performed by the migrator change the precision letter (e.g.
    ``D``→``Q`` for KIND, ``D``→``DD`` for multifloats); all other
    characters of the identifier are identical modulo case. So we look
    up the new name by the matched text's lowercase form and then
    transfer each source character's case onto the replacement.

    For single-char swaps (D→Q) the source and target are the same
    length and a positional case copy works directly. For multi-char
    expansions (D→DD) the target is longer; we find the first
    differing position, take the case of the original precision letter
    there, and apply it to all the inserted characters. Suffix after
    the old prefix is copied positionally from the source.

    Examples:
        ``dgemm_`` → ``qgemm_``      (KIND, equal length)
        ``DGEMM_`` → ``QGEMM_``
        ``PB_Cctypeset`` → ``PB_Cxtypeset``
        ``PB_Cdtypeset`` → ``PB_Cddtypeset``  (multifloats, expansion)
    """
    def _sub(m: re.Match) -> str:
        src = m.group(0)
        # Honor the original ``(?<![.>])`` lookbehind: skip identifier
        # tokens preceded by ``.`` or ``->`` (struct field access).
        start = m.start()
        if start > 0 and m.string[start - 1] in '.>':
            return src
        # PBTYP_T BLACS binding fields: protect the bare member
        # declaration (see _PBTYP_C_BINDING_FIELDS above). Matched
        # case-sensitively so genuine ``cgesd2d_`` / ``CGESD2D`` renames
        # still fire.
        if src in _PBTYP_C_BINDING_FIELDS:
            return src
        new_lower = combined.get(src.lower())
        if new_lower is None:
            return src

        if len(src) == len(new_lower):
            # Equal-length: positional case transfer.
            return ''.join(
                c.upper() if s.isupper() else c
                for s, c in zip(src, new_lower)
            )

        # Multi-char prefix expansion: locate the precision letter, copy
        # its case onto the inserted target chars, and stitch the
        # unchanged head and tail back in with positional case transfer.
        src_lower = src.lower()
        i = 0
        end = min(len(src_lower), len(new_lower))
        while i < end and src_lower[i] == new_lower[i]:
            i += 1
        # i is the first differing position. Assume a single-char source
        # prefix is being replaced with N new chars (the only pattern
        # the migrator currently emits): N = len(new) - len(src) + 1.
        n_old = 1
        n_new = len(new_lower) - len(src_lower) + n_old

        head = ''.join(
            c.upper() if s.isupper() else c
            for s, c in zip(src[:i], new_lower[:i])
        )
        # The inserted prefix takes its case from the original precision
        # letter at src[i]. For the typical lowercase identifier this is
        # lowercase; for an UPPERCASE call site (xerbla strings, name
        # mangling defines) it becomes uppercase consistently.
        new_prefix = new_lower[i:i + n_new]
        if i < len(src) and src[i].isupper():
            new_prefix = new_prefix.upper()

        tail_src = src[i + n_old:]
        tail_new = new_lower[i + n_new:]
        tail = ''.join(
            c.upper() if s.isupper() else c
            for s, c in zip(tail_src, tail_new)
        )

        return head + new_prefix + tail
    return _sub


def _expand_template(s: str, template_vars: dict[str, str]) -> str:
    """Expand ``{KEY}`` placeholders in ``s`` using template_vars."""
    for key, val in template_vars.items():
        s = s.replace('{' + key + '}', val)
    return s


# Cost-estimate local variable names used by PBLAS Level-3 entry points
# (pdgemm/pdsymm/pdsyrk/...) for algorithm selection. These are pure
# heuristic doubles that must NOT be promoted to the multifloats struct
# type, otherwise (double) casts and `*=` arithmetic in the cost model
# stop compiling. Survey of pblas/SRC/p[dz]*.c shows just two declaration
# lines containing these names; recognising them by name is sufficient.
#
# Note: ``tmp\d+`` deliberately excluded — PBLAS pairs them with
# ABest/ACest/BCest on the same line (so the line still matches via
# those anchors), while XBLAS uses bare ``double tmp1;`` accumulators
# that DO need to be promoted to float64x2 for working-precision math.
_PBLAS_COST_LOCAL = (
    r'ABest|ACest|BCest|ABestL|ABestR|Best'
)


def _apply_aliases_to_original(text: str, template_vars: dict[str, str],
                                c_type_aliases: list[dict] | None,
                                c_pointer_cast_aliases: list[dict] | None) -> str:
    """Apply recipe-declared aliases to a copy-original C source.

    Limited to (a) type-name aliases (e.g. ``cmplx16`` → ``cmplxQ``) and
    (b) pointer-cast aliases (e.g. ``(double*)`` → ``(quad*)``). Does
    NOT apply the broad ``double``/``float`` → ``REAL_TYPE`` substitution
    that :func:`_apply_c_type_subs` performs because copy-originals
    frequently contain precision-dispatch logic — e.g. ``PB_Cconjg``
    switches on ``TYPE->type`` and uses ``(double*)`` casts in the
    DCPLX/DREAL branches and ``(float*)`` casts in the SCPLX/SREAL
    branches. The bare ``double`` / ``float`` keywords inside those
    branches must stay so the dispatch stays well-formed; only the
    cast-stride needs upgrading so the kind16 (cmplxQ, 32-byte) target
    receives a 16-byte stride per real component instead of 8.
    """
    for rule in c_type_aliases or []:
        target = _expand_template(rule['to'], template_vars)
        for src in rule['from']:
            text = re.sub(r'\b' + re.escape(src) + r'\b', target, text)
    for rule in c_pointer_cast_aliases or []:
        target = _expand_template(rule['to'], template_vars)
        for src in rule['from']:
            text = text.replace(src, target)
    return text


def _apply_c_type_subs(text: str, template_vars: dict[str, str],
                      aliases: list[dict] | None = None,
                      target_mode: TargetMode | None = None) -> str:
    """Upgrade C type names used by precision-specific source files.

    ``aliases`` is a list of recipe-level rename rules of the form
    ``{'from': [names...], 'to': '<target>'}``. The target may contain
    ``{KEY}`` placeholders that expand from ``template_vars``. Applied
    after the built-in double/float/SCOMPLEX/DCOMPLEX substitutions.

    For multifloats targets we additionally:

    - Protect ``(double)`` cast expressions and the standard PBLAS
      cost-estimate local declarations (``double ABest, ACest, ...``)
      from the broad ``double`` -> ``float64x2_t`` substitution. These
      are heuristic algorithm-selection scalars that must stay as
      primitive doubles so the surrounding ``*=``, comparison and
      ``MAX(...)`` arithmetic continues to compile.
    - Rewrite the PBLAS scalar quick-return idioms (``ALPHA[REAL_PART]
      == ZERO`` etc.) into ``MF_IS_ZERO`` / ``MF_IS_ONE`` macro calls,
      because C ``==`` is undefined on the float64x2_t struct.
    """
    is_multifloats = target_mode is not None and not target_mode.is_kind_based

    # Protect cost-estimate idioms before the broad sub.
    cast_marker = '\x00MF_DOUBLE_CAST\x00'
    decl_marker = '\x00MF_DOUBLE_KW\x00'
    if is_multifloats:
        # (double) cast expressions. Negative lookbehind excludes
        # sizeof(double) / alignof(double), where the parens form a
        # sizeof argument rather than a cast — those must promote to
        # sizeof(float64x2_t) so heap allocations get the right size.
        text = re.sub(r'(?<!sizeof)(?<!alignof)\(\s*double\s*\)',
                      cast_marker, text)
        # Local declaration lines for cost-estimate locals. Match the
        # leading 'double' keyword on a line that mentions one of the
        # known cost-estimate names somewhere on the same line.
        text = re.sub(
            r'(?m)^(\s*)double(\s+[^;\n]*\b(?:'
            + _PBLAS_COST_LOCAL
            + r')\b[^;\n]*;)',
            lambda m: m.group(1) + decl_marker + m.group(2),
            text)

    text = re.sub(r'\bdouble\b', template_vars['REAL_TYPE'], text)
    text = re.sub(r'\bfloat\b', template_vars['REAL_TYPE'], text)
    text = re.sub(r'\bDCOMPLEX\b', template_vars['COMPLEX_TYPE'], text)
    text = re.sub(r'\bSCOMPLEX\b', template_vars['COMPLEX_TYPE'], text)
    for rule in aliases or []:
        target = _expand_template(rule['to'], template_vars)
        for src in rule['from']:
            text = re.sub(r'\b' + re.escape(src) + r'\b', target, text)

    if is_multifloats:
        text = text.replace(cast_marker, '(double)')
        text = text.replace(decl_marker, 'double')
        text = _apply_multifloats_pblas_subs(text)
    return text


_MF_PBLAS_PART = r'(\w+)\s*\[\s*(REAL_PART|IMAG_PART)\s*\]'


def _apply_multifloats_pblas_subs(text: str) -> str:
    """Rewrite PBLAS scalar quick-return checks for the float64x2_t
    struct type. C operators ``==`` / ``!=`` are not defined on
    structs, so we replace ``ALPHA[REAL_PART] == ZERO`` with the
    inline macro ``MF_IS_ZERO(ALPHA[REAL_PART])`` (and likewise for
    ``ONE`` / ``IMAG_PART`` / negation).
    """
    # ZERO comparisons
    text = re.sub(
        _MF_PBLAS_PART + r'\s*==\s*ZERO\b',
        r'MF_IS_ZERO(\1[\2])', text)
    text = re.sub(
        _MF_PBLAS_PART + r'\s*!=\s*ZERO\b',
        r'(!MF_IS_ZERO(\1[\2]))', text)
    # ONE comparisons
    text = re.sub(
        _MF_PBLAS_PART + r'\s*==\s*ONE\b',
        r'MF_IS_ONE(\1[\2])', text)
    text = re.sub(
        _MF_PBLAS_PART + r'\s*!=\s*ONE\b',
        r'(!MF_IS_ONE(\1[\2]))', text)
    # Integer truncation. PBLAS packs integer indices into work
    # buffers: ``(Int)(work[1])`` or ``(Int)(work[1][REAL_PART])``
    # becomes invalid on struct types. Rewrite to ``mf_to_int(...)``
    # (defined in both multifloats_c.h and the C++ bridge header).
    text = re.sub(
        r'\(Int\)\(\s*(\w+(?:\[\s*\w+\s*\])+)\s*\)',
        r'mf_to_int(\1)', text)
    return text


def migrate_c_file_to_string(
    src_path: Path,
    target_mode: TargetMode,
    rename_map: dict[str, str] | None = None,
    classification: SymbolClassification | None = None,
    c_type_aliases: list[dict] | None = None,
) -> tuple[str, str] | None:
    """Migrate one C source file in memory — no disk I/O.

    Mirrors :func:`migrate_file_to_string` for Fortran: returns
    ``(target_filename, migrated_text)`` or ``None`` when the file is
    precision-independent / not part of any family. Used by the
    convergence report to re-derive the target from an S/C sibling and
    compare with the D/Z-derived canonical already on disk.

    Two modes, mirroring :func:`migrate_c_directory`:

    - **Generic/scalapack** (both ``rename_map`` and ``classification``
      supplied): file is cloned iff its routine is a family member;
      in-text renames from ``rename_map`` plus C type upgrades applied.
    - **BLACS/direct** (both omitted): file's leading precision prefix
      (``d``/``z``/``s``/``c``, possibly preceded by ``BI_``) selects a
      substitution rule set; routine name is rewritten in-place.
    """
    if src_path.suffix.lower() != '.c':
        return None

    template_vars = _build_sub_vars(target_mode)

    if rename_map is not None and classification is not None:
        # Scalapack mode (PBLAS).
        stem = src_path.stem
        has_underscore = stem.endswith('_')
        routine = stem[:-1] if has_underscore else stem
        upper_routine = routine.upper()
        if upper_routine not in rename_map:
            return None
        if classification.get_family(upper_routine) is None:
            return None

        target_upper = rename_map[upper_routine]
        target_lower = target_upper.lower()
        new_stem = target_lower + ('_' if has_underscore else '')
        new_name = new_stem + src_path.suffix

        pattern, combined = _build_rename_regex(rename_map)
        sub = _make_rename_substituter(pattern, combined)
        text = src_path.read_text(errors='replace')
        text = pattern.sub(sub, text)
        text = _apply_c_type_subs(text, template_vars, c_type_aliases,
                                  target_mode=target_mode)
        return new_name, text

    # Direct/BLACS mode.
    stem = src_path.stem
    rp = template_vars['RP']
    cp = template_vars['CP']

    # Identify precision variant from stem prefix. Order matters:
    # check BI_-prefixed names before bare single-letter prefixes.
    if stem.startswith('BI_d'):
        src_prefix, new_prefix, subs = 'd', rp, REAL_CLONE_SUBS
    elif stem.startswith('BI_z'):
        src_prefix, new_prefix, subs = 'z', cp, COMPLEX_CLONE_SUBS
    elif stem.startswith('BI_s'):
        src_prefix, new_prefix, subs = 's', rp, SINGLE_CLONE_SUBS
    elif stem.startswith('BI_c'):
        src_prefix, new_prefix, subs = 'c', cp, CSINGLE_CLONE_SUBS
    elif stem.startswith('BI_'):
        return None  # precision-independent BI_* helper
    elif stem.startswith('d'):
        src_prefix, new_prefix, subs = 'd', rp, REAL_CLONE_SUBS
    elif stem.startswith('z'):
        src_prefix, new_prefix, subs = 'z', cp, COMPLEX_CLONE_SUBS
    elif stem.startswith('s'):
        src_prefix, new_prefix, subs = 's', rp, SINGLE_CLONE_SUBS
    elif stem.startswith('c'):
        src_prefix, new_prefix, subs = 'c', cp, CSINGLE_CLONE_SUBS
    else:
        return None

    new_name = rename_c_file(src_path.name, src_prefix, new_prefix)
    if new_name == src_path.name:
        return None

    renames = _routine_renames(stem, Path(new_name).stem)
    text = src_path.read_text(errors='replace')
    for pat, repl in subs:
        expanded = _expand_template(repl, template_vars)
        text = re.sub(pat, expanded, text, flags=re.MULTILINE)
    # See clone_c_file for the multi-char-prefix lookbehind rationale.
    for old_name, new_routine in renames:
        text = re.sub(
            r'(?<![A-Za-z0-9_])' + re.escape(old_name),
            new_routine, text)
        text = re.sub(
            r'(?<![A-Za-z0-9_])' + re.escape(old_name.upper()),
            new_routine.upper(), text)
    return new_name, text



def _duplicate_header_lines(text: str,
                            rename_pattern: re.Pattern,
                            rename_sub,
                            type_transform=None) -> str:
    """Duplicate header lines that mention precision-family identifiers.

    For each line (or multi-line declaration block) containing an
    identifier in the rename map, keep the original verbatim AND
    append a copy with rename + ``type_transform`` applied. The type
    transform is applied to the duplicated lines only so that the
    original ``void pdcopy_(double *)`` declaration coexists with a
    new ``void pddcopy_(float64x2_t *)`` declaration without
    introducing a conflicting redeclaration of either name.

    Multi-line declarations are detected by an open paren count: a
    line that opens an unbalanced ``(`` continues until the matching
    ``)`` is reached.
    """
    out_lines: list[str] = []
    lines = text.splitlines(keepends=True)
    i = 0
    while i < len(lines):
        line = lines[i]
        new_line = rename_pattern.sub(rename_sub, line)
        if new_line == line:
            out_lines.append(line)
            i += 1
            continue
        block: list[str] = [line]
        new_block: list[str] = [new_line]
        depth = line.count('(') - line.count(')')
        j = i + 1
        while depth > 0 and j < len(lines):
            block.append(lines[j])
            new_block.append(rename_pattern.sub(rename_sub, lines[j]))
            depth += lines[j].count('(') - lines[j].count(')')
            j += 1
        out_lines.extend(block)
        if type_transform is not None:
            new_block = [type_transform(b) for b in new_block]
        out_lines.extend(new_block)
        i = j
    return ''.join(out_lines)


def _apply_header_patches(output_dir: Path,
                          patches: list[dict],
                          template_vars: dict[str, str],
                          target_mode: TargetMode | None = None) -> None:
    """Insert recipe-declared lines into migrated headers.

    Each patch has ``file`` (relative name under output_dir), ``after``
    (literal anchor line that must be present exactly once) and
    ``insert`` (text to insert on the line after the anchor). The
    ``insert`` text is template-expanded.

    Patches can be gated by an optional ``when`` field whose value is
    matched against the active target_mode:

      - ``when: kind``        - applied for any KIND target (10 / 16)
      - ``when: multifloats`` - applied for multifloats target only
      - ``when: <name>``      - applied if target_mode.name == name
      - (absent)              - always applied (legacy behavior)
    """
    for patch in patches:
        when = patch.get('when')
        if when is not None and target_mode is not None:
            if when == 'kind':
                if not target_mode.is_kind_based:
                    continue
            elif when == 'multifloats':
                if target_mode.is_kind_based:
                    continue
            elif when != target_mode.name:
                continue
        target = output_dir / patch['file']
        if not target.exists():
            continue
        insert = _expand_template(patch['insert'], template_vars).rstrip('\n')
        text = target.read_text(errors='replace')
        if insert in text:
            continue  # already patched (idempotent)
        # Insertion modes:
        #  - ``at_bof: true`` — prepend at file start
        #  - ``at_eof: true`` — append at file end
        #  - ``before: <anchor>`` — insert on the line before the anchor
        #  - ``after: <anchor>`` — insert on the line after the anchor
        if patch.get('at_bof'):
            text = insert + '\n' + text
        elif patch.get('at_eof'):
            sep = '' if text.endswith('\n') else '\n'
            text = text + sep + insert + '\n'
        else:
            before = patch.get('before')
            anchor = before or patch.get('after')
            if anchor is None or anchor not in text:
                continue
            if before is not None:
                text = text.replace(anchor, insert + '\n' + anchor, 1)
            else:
                text = text.replace(anchor, anchor + '\n' + insert, 1)
        target.write_text(text)


def _migrate_generic_c_directory(src_dir: Path, output_dir: Path,
                                 target_mode: TargetMode,
                                 classification: SymbolClassification,
                                 rename_map: dict[str, str],
                                 c_type_aliases: list[dict] | None = None,
                                 c_pointer_cast_aliases: list[dict] | None = None,
                                 header_patches: list[dict] | None = None,
                                 extra_c_dirs: list[Path] | None = None,
                                 skip_files: set[str] | None = None,
                                 copy_files: set[str] | None = None) -> dict:
    """Rename-map-driven C migration for ScaLAPACK-style libraries.

    A file ``foo_.c`` is cloned iff its routine ``FOO`` is a D- or
    Z-precision member of some family in ``classification``. The clone
    is written to ``<target>.c`` where ``<target>`` is the lowercase of
    ``rename_map[FOO]``, with all in-file routine names rewritten via
    ``rename_map`` and C types upgraded to the target precision.

    ``extra_c_dirs`` is a list of additional directories whose .c
    sources are flat-copied (and cloned-as-applicable) into
    ``output_dir`` alongside ``src_dir``. Used by PBLAS to migrate the
    PTOOLS/ helpers.
    """
    template_vars = _build_sub_vars(target_mode)
    output_dir.mkdir(parents=True, exist_ok=True)

    all_src_dirs: list[Path] = [src_dir]
    if extra_c_dirs:
        all_src_dirs.extend(extra_c_dirs)

    # Stage the pass-through files first (headers, copy_files entries,
    # precision-independent dispatchers).
    # When extra_c_dirs sources contain `#include "../foo.h"` paths
    # (PTOOLS uses ../pblas.h etc.), strip the `..` part since we're
    # flattening into a single output dir.
    # K&R-style function definitions are converted to ANSI so originals
    # compile as C++ (mpicxx).
    _skip = skip_files or set()
    _copy = copy_files or set()
    # Only headers, copy_files entries, and precision-independent
    # dispatcher .c files are staged. The std archive (built directly
    # from upstream sources by the CMake side) carries the S/D/C/Z
    # entry points; the migrated archive carries the Q/X/E/Y/M/W
    # clones plus the dispatchers (PB_Cconjg, PB_CpswapNN, BI_BlacsErr,
    # ...) — files whose stems are NOT in rename_map. Dispatchers must
    # ride in the migrated archive because the migrator's
    # _apply_aliases_to_original pass widens their ``(double*)`` /
    # ``(float*)`` pointer-casts to ``(QREAL*)`` / ``(EREAL*)`` /
    # ``(float64x2_t*)`` for KIND targets, so the byte-stride
    # arithmetic matches the wider type carried by callers. The std
    # archive's untouched ``(double*)`` copies would corrupt strides
    # when the migrated entry points dispatch through them.
    for d in all_src_dirs:
        for f in (p for p in sorted(d.iterdir())):
            stem_upper = (f.stem[:-1] if f.stem.endswith('_')
                          else f.stem).upper()
            # Also consider the decoration-stripped stem for skip
            # matching. ``BLAS_cdot_c_s-f2c.c`` is a bridge for the
            # already-skipped ``BLAS_cdot_c_s.c``; both should drop.
            base_stem_for_skip, _deco_for_skip = _strip_decoration(
                f.stem[:-1] if f.stem.endswith('_') else f.stem
            )
            base_stem_upper = base_stem_for_skip.upper()
            is_c_or_h = f.suffix.lower() in ('.c', '.h')
            is_header = f.suffix.lower() == '.h'
            is_copy = stem_upper in _copy
            # A .c file is a "dispatcher" iff it is NOT precision-
            # prefixed: stem not in rename_map AND not matching the
            # ``p<sdcz><root>`` pattern that ``_redist_clone_stem``
            # uses to clone files whose Fortran-callable exports live
            # behind ``#define`` macros (REDIST/SRC entry points like
            # pdgemr.c, pztrmr2.c, ScaLAPACK orphan helpers like
            # pdlaiect.c). Both classes are owned by the std archive;
            # only true dispatchers (PB_Cconjg, BI_BlacsErr, ...) ride
            # in the migrated archive.
            is_dispatcher = False
            if (f.suffix.lower() == '.c'
                    and stem_upper not in (rename_map or {})):
                # XBLAS f2c-bridge files (``BLAS_dgemv_x-f2c.c``) are
                # not in rename_map under their decorated stem, but
                # the cloning loop knows how to rewrite them via
                # ``_strip_decoration``. Don't classify them as
                # dispatchers; the cloning pass will produce the
                # correctly-renamed output and the original copy
                # would just collide with that.
                base_stem, deco = _strip_decoration(
                    f.stem[:-1] if f.stem.endswith('_') else f.stem
                )
                if (deco
                        and base_stem.upper() in (rename_map or {})):
                    pass  # cloned by the rename-map pass
                elif _redist_clone_stem(f.stem, target_mode) is None:
                    is_dispatcher = True
            if not (is_c_or_h or is_copy):
                continue
            if stem_upper in _skip or base_stem_upper in _skip:
                continue
            if not (is_header or is_copy or is_dispatcher):
                continue
            text = f.read_text(errors='replace')
            if is_copy and not is_c_or_h:
                # Fortran copy-files are staged verbatim: no include
                # rewrites, no K&R conversion, no stdc-ifdef folding.
                (output_dir / f.name).write_text(text)
                continue
            if d != src_dir:
                text = re.sub(
                    r'#include\s*"\.\./([^"]+)"',
                    r'#include "\1"', text)
            if f.suffix.lower() == '.c':
                # Aliases (cmplx16 → cmplxQ, (double*) → (REAL_TYPE*))
                # apply to:
                #   (a) precision-independent dispatchers (not in
                #       rename_map, e.g. PB_Cconjg, PB_Ctzher2k) —
                #       on EVERY target, because cloned entry points
                #       call them and need the wider types.
                #   (b) precision-prefixed originals (in rename_map,
                #       e.g. pdgemm_, pcamax_, PB_Cdtypeset) — only on
                #       KIND-based targets, where -freal-8-real-16
                #       promotion means the original is still called
                #       with quad-stride args and its body must
                #       widen too. On multifloats targets the
                #       precision-prefixed originals are dead (callers
                #       route through clones) and their native
                #       double*/float* signatures must stay so the
                #       body still compiles as C++ (struct types
                #       cannot assign to native scalar lvalues).
                is_unprefixed = stem_upper not in rename_map
                if is_unprefixed or target_mode.is_kind_based:
                    text = _apply_aliases_to_original(
                        text, template_vars, c_type_aliases,
                        c_pointer_cast_aliases)
            (output_dir / f.name).write_text(text)

    # Snapshot the pristine (staged, pre-transform) bytes of every
    # copied header. Any header the transforms below modify is emitted
    # under a precision-prefixed sibling name (see the split step before
    # the return) and the original is restored to these bytes, so every
    # installed public Netlib-named header stays byte-identical to the
    # upstream reference across all precisions.
    _pristine_headers = {
        h.name: h.read_bytes()
        for h in output_dir.iterdir()
        if h.suffix.lower() == '.h'
    }

    # Apply recipe-declared header patches to the copied originals so
    # later clones can reference the newly introduced typedefs.
    if header_patches:
        _apply_header_patches(output_dir, header_patches, template_vars,
                              target_mode=target_mode)

    rename_pattern, combined_map = _build_rename_regex(rename_map)
    _rename_sub = _make_rename_substituter(rename_pattern, combined_map)

    def _rename(text: str) -> str:
        return rename_pattern.sub(_rename_sub, text)

    # Header rename pass: each header in output_dir gets its precision-
    # family declarations duplicated. Lines that mention an identifier
    # in the rename map are kept verbatim AND a copy with the rename
    # applied is appended right after. This propagates BLACS / PBLAS
    # function prototypes to both the original (Cdgesd2d) and the
    # cloned (Cddgesd2d) name so cloned source files compile.
    # The duplicated lines also get type subs (double -> float64x2_t,
    # etc.) so the cloned signature uses the right types.
    # Only header declarations are duplicated, not the .c sources --
    # those go through the per-file rename + clone path below.
    def _hdr_type_transform(line: str) -> str:
        return _apply_c_type_subs(line, template_vars, c_type_aliases,
                                  target_mode=target_mode)

    for hdr in sorted(output_dir.iterdir()):
        if hdr.suffix.lower() != '.h':
            continue
        text = hdr.read_text(errors='replace')
        new_text = _duplicate_header_lines(
            text, rename_pattern, _rename_sub,
            type_transform=_hdr_type_transform)
        if new_text != text:
            hdr.write_text(new_text)

    # Insert real-type typedef into pblas.h for KIND targets so that
    # migrated declarations (QREAL, QCOMPLEX etc.) resolve.
    pblas_path = output_dir / 'pblas.h'
    if pblas_path.exists() and target_mode.c_header_mode == 'typedef':
        _patch_pblas_header(pblas_path, template_vars)

    # Process D/Z-sourced files first so they become the canonical
    # output; S/C co-family members are verified against them.
    def _is_double_key(routine_upper: str) -> bool:
        fam = classification.get_family(routine_upper)
        if fam is None:
            return False
        key = next(
            (k for k, v in fam.members.items() if v == routine_upper),
            '',
        )
        return bool(key) and key[0] in ('D', 'Z')

    all_c_files: list[Path] = []
    for d in all_src_dirs:
        all_c_files.extend(f for f in d.iterdir()
                            if f.suffix.lower() == '.c')
    entries = sorted(
        all_c_files,
        key=lambda f: (
            not _is_double_key(
                (f.stem[:-1] if f.stem.endswith('_') else f.stem).upper()
            ),
            f.name,
        ),
    )

    def _normalize_for_compare(s: str) -> str:
        # Strip block + line comments
        s = re.sub(r'/\*.*?\*/', '', s, flags=re.DOTALL)
        s = re.sub(r'//[^\n]*', '', s)
        # Canonicalize prefix-dependent identifiers: leading s/d/c/z
        # becomes '@' so sibling C sources that differ only in the
        # precision prefix of local names collapse together.
        s = re.sub(r'\b[sdczSDCZ]+(?=[A-Za-z])', '@', s)
        # Tokenize each line into words and single non-whitespace
        # characters, then rejoin with a single space. This collapses
        # not just column-aligned padding like ``float          *``
        # vs ``double         *`` but also incidental single-space
        # drift around punctuation: ``(QREAL )(`` vs ``(QREAL)(``,
        # ``if (`` vs ``if(``, ``Mptr(...)))`` vs ``Mptr(... ))``.
        lines = [' '.join(re.findall(r'\w+|\S', ln)) for ln in s.split('\n')]
        return '\n'.join(ln for ln in lines if ln)

    cloned: list[str] = []
    divergences: list[str] = []
    canonical_normalized: dict[str, str] = {}
    canonical_source: dict[str, str] = {}

    for f in entries:
        has_underscore = f.stem.endswith('_')
        routine = f.stem[:-1] if has_underscore else f.stem
        # XBLAS bridge files (``BLAS_dgemv_x-f2c.c``) share the
        # routine stem of their parent ``BLAS_dgemv_x`` plus a
        # decoration; strip it for the rename-map lookup and reapply
        # to the output filename.
        base_routine, decoration = _strip_decoration(routine)
        upper_routine = base_routine.upper()

        if upper_routine in _skip or routine.upper() in _skip:
            continue

        # REDIST-style files like ``pdgemr.c`` / ``pdgemr2.c`` /
        # ``pdtrmr2.c`` expose Fortran-callable entry points only behind
        # ``#define fortran_mr2dnew pdgemr2d_`` macros, so the file
        # stem doesn't match an exported symbol and the symbol scanner
        # never adds it to rename_map. The stem itself is still
        # ``p<precision-letter><root>``, so derive the clone stem
        # directly via the target's precision-family map.
        if (upper_routine in rename_map
                and classification.get_family(upper_routine) is not None):
            new_stem = rename_map[upper_routine].lower() + decoration
        else:
            new_stem = _redist_clone_stem(base_routine, target_mode)
            if new_stem is None:
                continue
            new_stem = new_stem + decoration
        new_name = new_stem + ('_' if has_underscore else '') + f.suffix
        new_path = output_dir / new_name

        text = f.read_text(errors='replace')
        # Files from extra_c_dirs (e.g. PBLAS PTOOLS/) are flattened
        # into output_dir, so #include "../foo.h" must become "foo.h".
        if f.parent != src_dir:
            text = re.sub(
                r'#include\s*"\.\./([^"]+)"',
                r'#include "\1"', text)
        # Apply renames first, then type upgrades — the two domains
        # don't overlap but this order preserves identifier names that
        # happen to coincide with generic type keywords.
        text = _rename(text)
        text = _apply_c_type_subs(text, template_vars, c_type_aliases,
                                  target_mode=target_mode)

        # PBLAS typesets assign the type-generic drivers' BLAS callbacks
        # (PBTYP_T.F<slot>) to migrated gfortran leaves. The generic C
        # drivers call these through prototypes that pass NO hidden Fortran
        # CHARACTER lengths, but the gfortran leaves spill scratch into
        # those (unprovisioned) slots and corrupt the caller frame. Route
        # every char-taking callback through a hidden-length trampoline and
        # drop the shared PBcharshim.h substrate alongside the clone. The
        # MKL-backed s/c/d/z typesets are never cloned, so they never see
        # this pass.
        if is_typeset_stem(new_stem):
            text, _shim_changed = transform_typeset(text)
            if _shim_changed:
                (output_dir / 'PBcharshim.h').write_text(
                    pbcharshim_header_text())

        normalized = _normalize_for_compare(text)
        prior = canonical_normalized.get(new_name)
        if prior is None:
            new_path.write_text(text)
            canonical_normalized[new_name] = normalized
            canonical_source[new_name] = f.name
            cloned.append(f'{f.name} → {new_name}')
        elif prior == normalized:
            pass  # convergence (ignoring comment differences)
        else:
            divergences.append(
                f'{f.name} vs {canonical_source[new_name]} → {new_name}'
            )

    # --- split transformed headers into precision-prefixed siblings ---
    # Every header the migrator changed above (recipe patches, prototype
    # duplication, precision typedefs) is written to a sibling named
    # ``<pair>name.h`` (pair = real+complex prefix, e.g. ``qx``/``ey``/
    # ``mw``); the original ``name.h`` is restored to its pristine
    # upstream bytes. The include-rewrite (in ``migrate_c_directory``,
    # after overrides + the C++ extern-wrap) points migrated sources at
    # the sibling. Result: the installed public ``name.h`` is exactly the
    # Netlib reference, while this precision's build sees the sibling.
    #
    # The split set is closed under the ``#include`` relation: a header
    # that includes a split header must itself be split, otherwise a
    # migrated TU that reaches it (e.g. ``blas_extended.h`` including
    # ``blas_extended_proto.h``) would see the restored-pristine original,
    # which has lost the transforms (missing ``extern "C"`` / typedefs) and
    # fails to compile. Walk the relation to a fixed point.
    pair = template_vars['RP'] + template_vars['CP']
    # Guard idempotency: a prior stage into this tree may have left the
    # precision-prefixed siblings behind. Never treat an already-prefixed
    # sibling as a splittable original, else re-staging doubles the prefix
    # (mwblas_extended.h -> mwmwblas_extended.h).
    split_names = {
        name for name, pristine in _pristine_headers.items()
        if not name.startswith(pair)
        and (output_dir / name).read_bytes() != pristine
    }
    grew = True
    while grew:
        grew = False
        for name in _pristine_headers:
            if name in split_names or name.startswith(pair):
                continue
            text = (output_dir / name).read_text(errors='replace')
            for target in split_names:
                esc = re.escape(target)
                if re.search(r'#\s*include\s*("' + esc + r'"|<' + esc + r'>)',
                             text):
                    split_names.add(name)
                    grew = True
                    break
    split_headers: dict[str, str] = {}
    for name in split_names:
        sibling = pair + name
        # Read the current (transformed for changed headers, pristine for
        # pulled-in intermediates) content before restoring the original.
        (output_dir / sibling).write_bytes((output_dir / name).read_bytes())
        (output_dir / name).write_bytes(_pristine_headers[name])
        split_headers[name] = sibling

    return {
        'cloned': cloned,
        'divergences': divergences,
        'template_vars': template_vars,
        'split_headers': split_headers,
    }


def _migrate_blacs_c_directory(src_dir: Path, output_dir: Path,
                                 target_mode: TargetMode,
                                 skip_files: set[str] | None = None) -> dict:
    """BLACS-specific C migration (original behavior).

    Reads from ``src_dir`` (the staged source tree, already patched).

    ``skip_files`` (uppercased logical stems, no trailing Fortran
    underscore) are dropped from the migrated output — used to remove the
    genuine S/C entry points (``sgesd2d_``, ``BI_cvvsum``, …) that the
    plain ``blacs`` standard archive and MKL already provide, so the
    migrated archive stays strictly extended-precision. Skipped stems
    are also never used as clone sources (they are S/C, not the D/Z the
    clone loops read), so dropping them cannot disturb the e/y clones.
    """
    _skip = skip_files or set()
    template_vars = _build_sub_vars(target_mode)
    rp = template_vars['RP']
    cp = template_vars['CP']

    output_dir.mkdir(parents=True, exist_ok=True)

    sources = sorted(src_dir.iterdir())

    # Copy headers always; drop precision-prefixed .c originals (the
    # std archive owns the S/D/C/Z entry points); copy precision-
    # independent dispatcher .c files (BI_BlacsErr, BI_GetBuff,
    # blacs_setup_, ...). The
    # dispatchers ride in the migrated archive because callers reach
    # them through the BLACS-style symbol map; the std archive picks
    # them up too via the upstream source dir, but the linker only
    # ever pulls one definition (whichever resolves an undefined
    # symbol first), so duplication is benign and the migrated copies
    # match the precision the migrated bodies expect.
    def _is_blacs_precision_prefix(stem: str) -> bool:
        # Mirrors the d/z clone discovery in the loops below.
        if stem.startswith(('d', 'z')) and not stem.startswith(('BI_',)):
            return True
        if stem.startswith(('BI_d', 'BI_z')):
            return True
        return False

    for f in sources:
        ext = f.suffix.lower()
        # ``skip_files`` drops genuine S/C sources from the migrated
        # output entirely (headers are never precision-specific, so they
        # are exempt). Match on the underscore-stripped, uppercased stem
        # to align with the recipe's logical-name convention.
        stem_upper = (f.stem[:-1] if f.stem.endswith('_') else f.stem).upper()
        if ext == '.c' and stem_upper in _skip:
            continue
        if ext == '.h':
            shutil.copy2(f, output_dir / f.name)
        elif ext == '.c':
            if not _is_blacs_precision_prefix(f.stem):
                shutil.copy2(f, output_dir / f.name)

    cloned = []

    # Clone d-variant → real-extended
    for f in sources:
        if f.suffix.lower() != '.c':
            continue
        stem = f.stem
        if stem.startswith('d') or stem.startswith('BI_d'):
            # Skip if it's not a type-specific file (check for 'd' prefix)
            if stem.startswith('BI_d') and stem[3] == 'd':
                pass  # BI_d* files
            elif stem.startswith('d') and not stem.startswith('BI_'):
                pass  # d* files
            else:
                continue

            new_name = rename_c_file(f.name, 'd', rp)
            if new_name == f.name:
                continue
            renames = _routine_renames(stem, Path(new_name).stem)
            clone_c_file(f, output_dir / new_name,
                         REAL_CLONE_SUBS, template_vars, renames)
            cloned.append(f'{f.name} → {new_name}')

    # Clone z-variant → complex-extended
    for f in sources:
        if f.suffix.lower() != '.c':
            continue
        stem = f.stem
        if stem.startswith('z') or stem.startswith('BI_z'):
            new_name = rename_c_file(f.name, 'z', cp)
            if new_name == f.name:
                continue
            renames = _routine_renames(stem, Path(new_name).stem)
            clone_c_file(f, output_dir / new_name,
                         COMPLEX_CLONE_SUBS, template_vars, renames)
            cloned.append(f'{f.name} → {new_name}')

    # Patch Bdef.h with extended-precision type definitions and macros
    bdef_path = output_dir / 'Bdef.h'
    if bdef_path.exists():
        _patch_bdef_header(bdef_path, target_mode, template_vars)


    # Generate MPI datatype availability check when the target requires
    # stock MPI datatypes that may not be universally available (e.g.
    # MPI_REAL16 / MPI_COMPLEX32 for KIND=16).
    if target_mode.c_needs_mpi_check:
        _generate_mpi_real16_check(output_dir)

    return {
        'cloned': cloned,
        'template_vars': template_vars,
    }


def _patch_pblas_header(pblas_path: Path,
                        template_vars: dict[str, str]) -> None:
    """Insert extended-precision typedefs into pblas.h for KIND targets.

    The migrator replaces ``double``/``complex16`` with the target type
    names (e.g. QREAL, QCOMPLEX) in .c files, but pblas.h keeps the
    original types.  We add typedefs so both old and new names compile.
    """
    real_type = template_vars['REAL_TYPE']
    complex_type = template_vars['COMPLEX_TYPE']
    c_type = template_vars['C_REAL_TYPE']

    block = (f'typedef {c_type} {real_type};\n'
             f'typedef struct {{ {real_type} re, im; }} {complex_type};\n')

    text = pblas_path.read_text(errors='replace')
    if real_type in text and 'typedef' in text.split(real_type)[0]:
        return  # already patched
    # Insert just before the first "typedef struct" line
    marker = 'typedef struct'
    idx = text.find(marker)
    if idx >= 0:
        text = text[:idx] + block + '\n' + text[idx:]
        pblas_path.write_text(text)


# The 11 BLACS user-facing routine suffixes (same for each type prefix).
_BLACS_ROUTINE_SUFFIXES = [
    'gesd2d', 'gerv2d', 'gebs2d', 'gebr2d',
    'trsd2d', 'trrv2d', 'trbs2d', 'trbr2d',
    'gsum2d', 'gamx2d', 'gamn2d',
]


def _patch_bdef_header(bdef_path: Path, target_mode: TargetMode,
                       template_vars: dict[str, str]) -> None:
    """Add extended-precision type definitions and macros to Bdef.h."""
    rp = template_vars['RP']
    cp = template_vars['CP']
    real_type = template_vars['REAL_TYPE']      # e.g. QREAL / float64x2_t
    complex_type = template_vars['COMPLEX_TYPE']  # e.g. XCOMPLEX / complex128x2_t

    c_type = template_vars['C_REAL_TYPE']

    if target_mode.c_header_mode == 'include':
        # Module-based target: types come from an external header
        # (e.g. multifloats_bridge.h) which provides full operator
        # overloading via C++ templates.
        type_block = f"""
/*
 *  Companion types for migrated {rp}/{cp} routines.
 *  Provided by {target_mode.c_header}.
 */
#include "{target_mode.c_header}"

void BI_{rp}mvcopy(Int m, Int n, {real_type} *A, Int lda, {real_type} *buff);
void BI_{rp}vmcopy(Int m, Int n, {real_type} *A, Int lda, {real_type} *buff);
"""
    else:
        # typedef mode (KIND targets): emit primitive typedefs.
        type_block = f"""
/*
 *  Extended-precision types for migrated {{prefix}} routines.
 */
typedef {c_type} {real_type};
typedef struct {{{real_type} r, i;}} {complex_type};

void BI_{rp}mvcopy(Int m, Int n, {real_type} *A, Int lda, {real_type} *buff);
void BI_{rp}vmcopy(Int m, Int n, {real_type} *A, Int lda, {real_type} *buff);
""".replace('{prefix}', f'{rp}/{cp}')

    # --- Complex copy macros ---
    macro_block = f"""#define BI_{cp}mvcopy(m, n, A, lda, buff) \\
        BI_{rp}mvcopy(2*(m), (n), ({real_type} *) (A), 2*(lda), ({real_type} *) (buff))
#define BI_{cp}vmcopy(m, n, A, lda, buff) \\
        BI_{rp}vmcopy(2*(m), (n), ({real_type} *) (A), 2*(lda), ({real_type} *) (buff))
"""

    # --- Fortran name mangling defines ---
    def _mangling_block(prefix: str, transform) -> str:
        lines = []
        for suf in _BLACS_ROUTINE_SUFFIXES:
            src = f'{prefix}{suf}_'
            dst = transform(f'{prefix}{suf}')
            lines.append(f'#define {src:19s}{dst}')
        return '\n'.join(lines) + '\n'

    nochange_block = (_mangling_block(rp, lambda s: s) +
                      _mangling_block(cp, lambda s: s))
    upcase_block = (_mangling_block(rp, str.upper) +
                    _mangling_block(cp, str.upper))

    text = bdef_path.read_text(errors='replace')

    # Insert type definitions after the DCOMPLEX/SCOMPLEX typedefs
    text = text.replace(
        'typedef struct {float r, i;} SCOMPLEX;\n',
        'typedef struct {float r, i;} SCOMPLEX;\n' + type_block
    )

    # Insert copy macros after the existing BI_zvmcopy macro
    zvmcopy_line = '#define BI_zvmcopy(m, n, A, lda, buff) \\\n        BI_dvmcopy(2*(m), (n), (double *) (A), 2*(lda), (double *) (buff))'
    text = text.replace(zvmcopy_line, zvmcopy_line + '\n' + macro_block)

    # Insert name mangling defines in NOCHANGE section (after zgamn2d)
    nochange_marker = '#define zgamn2d_   zgamn2d\n'
    text = text.replace(nochange_marker, nochange_marker + nochange_block)

    # Insert name mangling defines in UPCASE section (after ZGAMN2D)
    upcase_marker = '#define zgamn2d_   ZGAMN2D\n'
    text = text.replace(upcase_marker, upcase_marker + upcase_block)

    bdef_path.write_text(text)


def _generate_mpi_real16_check(output_dir: Path) -> None:
    """Generate a CMake module that verifies MPI_REAL16 support."""
    cmake = """\
# CheckMpiReal16.cmake — Verify that the MPI implementation provides
# MPI_REAL16 and MPI_COMPLEX32, which are required by the migrated
# quad-precision BLACS routines.
#
# Usage:
#   include(CheckMpiReal16)
#   check_mpi_real16()          # FATAL_ERROR if unsupported
#
# After a successful check the cache variable HAVE_MPI_REAL16 is set.

include(CheckCSourceCompiles)

function(check_mpi_real16)
    find_package(MPI REQUIRED COMPONENTS C)

    set(CMAKE_REQUIRED_INCLUDES  ${MPI_C_INCLUDE_DIRS})
    set(CMAKE_REQUIRED_LIBRARIES ${MPI_C_LIBRARIES})
    set(CMAKE_REQUIRED_FLAGS     ${MPI_C_COMPILE_FLAGS})

    check_c_source_compiles("
        #include <mpi.h>
        int main(void) {
            MPI_Datatype dt_real = MPI_REAL16;
            MPI_Datatype dt_cplx = MPI_COMPLEX32;
            (void)dt_real; (void)dt_cplx;
            return 0;
        }
    " HAVE_MPI_REAL16)

    if(NOT HAVE_MPI_REAL16)
        message(FATAL_ERROR
            "The MPI implementation does not provide MPI_REAL16 / MPI_COMPLEX32. "
            "Quad-precision BLACS requires an MPI library built with 128-bit real "
            "support (e.g. OpenMPI or MPICH compiled with __float128 enabled).")
    endif()

    set(HAVE_MPI_REAL16 ${HAVE_MPI_REAL16} PARENT_SCOPE)
endfunction()
"""
    (output_dir / 'CheckMpiReal16.cmake').write_text(cmake)
