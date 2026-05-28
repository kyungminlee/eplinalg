"""General-purpose Fortran type migration engine.

Handles fixed-form (.f) and free-form (.f90) Fortran source files.
Works for any library following the S/D/C/Z precision prefix convention
(BLAS, LAPACK, ScaLAPACK, MUMPS, etc.).

Uses a hybrid strategy (per DEVELOPER.md):
  1. Parse with Flang (via subprocess) to get a structured parse tree
  2. Extract facts: type declarations, routine names, call sites,
     literals, intrinsics, XERBLA strings
  3. Apply transformations as source-level replacements, preserving
     all original formatting, comments, and preprocessor directives

Falls back to regex-only scanning when Flang is not available.

Transformations:
  1. Type declarations (DOUBLE PRECISION → REAL(KIND=k), etc.)
  2. Standalone REAL/COMPLEX in declaration context
  3. Floating-point literals (D-exponent → E_kind suffix)
  4. Type-specific intrinsic function calls
  5. INTRINSIC declaration statements
  6. Routine name prefixes (using rename map)
  7. XERBLA string arguments
"""

import functools
import re
from pathlib import Path

from .intrinsics import INTRINSIC_MAP, INTRINSIC_DECL_MAP
from .target_mode import TargetMode


# ---------------------------------------------------------------------------
# Type declaration replacement
# ---------------------------------------------------------------------------

# Pre-compiled patterns for replace_type_decls. The 14 substitutions
# below were string-pattern ``re.sub`` calls; each went through
# Python's regex cache lookup, which profiling showed accumulated to
# ~5 s across 828k LAPACK invocations.
_TD_DECL_TAIL = r'(?=\s+[A-Za-z_]|\s*::|\s*,)'
_TD_DBL_PREC = re.compile(r'DOUBLE\s+PRECISION', re.IGNORECASE)
_TD_DBL_CMPLX = re.compile(r'DOUBLE\s+COMPLEX', re.IGNORECASE)
_TD_CMPLX_STAR16 = re.compile(r'COMPLEX\*16', re.IGNORECASE)
_TD_REAL_STAR8 = re.compile(r'REAL\*8', re.IGNORECASE)
_TD_CMPLX_STAR8 = re.compile(r'COMPLEX\*8', re.IGNORECASE)
_TD_REAL_STAR4 = re.compile(r'REAL\*4', re.IGNORECASE)
_TD_REAL_KIND_E0 = re.compile(
    r'REAL\s*\(\s*kind\s*\(\s*0\.[Ee]0\s*\)\s*\)', re.IGNORECASE)
_TD_CMPLX_KIND_E0 = re.compile(
    r'COMPLEX\s*\(\s*kind\s*\(\s*0\.[Ee]0\s*\)\s*\)', re.IGNORECASE)
_TD_REAL_KIND_D0 = re.compile(
    r'REAL\s*\(\s*kind\s*\(\s*0\.[Dd]0\s*\)\s*\)', re.IGNORECASE)
_TD_CMPLX_KIND_D0 = re.compile(
    r'COMPLEX\s*\(\s*kind\s*\(\s*0\.[Dd]0\s*\)\s*\)', re.IGNORECASE)
_TD_REAL_KIND_WP = re.compile(
    r'REAL\s*\(\s*(?:KIND\s*=\s*)?WP\s*\)', re.IGNORECASE)
_TD_CMPLX_KIND_WP = re.compile(
    r'COMPLEX\s*\(\s*(?:KIND\s*=\s*)?WP\s*\)', re.IGNORECASE)
_TD_REAL_KIND_4 = re.compile(
    r'REAL\s*\(\s*(?:KIND\s*=\s*)?4\s*\)' + _TD_DECL_TAIL, re.IGNORECASE)
_TD_CMPLX_KIND_4 = re.compile(
    r'COMPLEX\s*\(\s*(?:KIND\s*=\s*)?4\s*\)' + _TD_DECL_TAIL, re.IGNORECASE)
_TD_REAL_KIND_8 = re.compile(
    r'REAL\s*\(\s*(?:KIND\s*=\s*)?8\s*\)' + _TD_DECL_TAIL, re.IGNORECASE)
_TD_CMPLX_KIND_8 = re.compile(
    r'COMPLEX\s*\(\s*(?:KIND\s*=\s*)?8\s*\)' + _TD_DECL_TAIL, re.IGNORECASE)
# Cheap early-out gate: lines containing none of these tokens can't
# match any of the 14 type-decl patterns above, so the whole pass is
# skippable. Most LAPACK source lines (executable statements,
# comments, blank lines) hit this gate and return immediately.
_TD_GATE_RE = re.compile(r'\b(?:REAL|COMPLEX|DOUBLE|WP)\b', re.IGNORECASE)


def replace_type_decls(
    line: str,
    target_mode: TargetMode,
    complex_names: set[str] | None = None,
    source_kind: int | None = None,
) -> str:
    """Replace precision type keywords with target form.

    ``source_kind`` (4 / 8 / None) gates which patterns get promoted
    so a kind4 source half (S/C) only promotes kind4 type tokens and
    a kind8 source half (D/Z) only promotes kind8 type tokens. ``None``
    promotes every pattern (precision-independent files / cases where
    half is unknown). The split mirrors the LAPACK / ScaLAPACK / MUMPS
    convention: each precision half declares its working-precision
    variables in the half's native kind; cross-kind references inside
    a half (notably MUMPS S-half ``DOUBLE PRECISION`` for timing
    statistics that should remain kind8 even at kind16 retarget) are
    intentionally NOT working precision and must survive migration.

    In multifloats mode, also filters out variable names that are
    supplied as named constants by the multifloats module (e.g. ZERO,
    ONE, TWO). Those names become module imports and must not appear
    as locally-declared variables — otherwise gfortran complains:
    "Symbol 'mf_one' conflicts with symbol from module 'multifloats'".

    ``complex_names`` is the file-scope set of names declared as
    COMPLEX in any procedure. When provided, names in that set are
    NOT stripped from real declarations — ZERO, ONE etc. that occur
    in BOTH real and complex contexts in the same file are kept as
    locals everywhere, since the global rename to a real
    multifloats constant would mistype the COMPLEX scope.
    """
    # Most LAPACK lines have no type-decl token; bail before any sub.
    if not _TD_GATE_RE.search(line):
        if not target_mode.is_kind_based:
            line = _filter_known_constants_from_decl(
                line, target_mode, complex_names=complex_names,
            )
        return line

    real_target = target_mode.real_type
    complex_target = target_mode.complex_type

    promote_k8 = source_kind in (None, 8)
    promote_k4 = source_kind in (None, 4)

    # Longest patterns first to avoid partial matches.
    if promote_k8:
        line = _TD_DBL_PREC.sub(real_target, line)
        line = _TD_DBL_CMPLX.sub(complex_target, line)
        line = _TD_CMPLX_STAR16.sub(complex_target, line)
        line = _TD_REAL_STAR8.sub(real_target, line)
    if promote_k4:
        line = _TD_CMPLX_STAR8.sub(complex_target, line)
        line = _TD_REAL_STAR4.sub(real_target, line)
    # ``REAL(kind(0.E0))`` (single, kind4) / ``REAL(kind(0.D0))`` (double, kind8)
    # — an older idiom that predates F90's REAL*N form. MUMPS's
    # dmumps_struc.h uses ``REAL(kind(0.E0))`` for single-precision
    # fields. Gate by source kind so a kind4 source half preserves
    # ``REAL(kind(0.D0))`` references and vice versa.
    if promote_k4:
        line = _TD_REAL_KIND_E0.sub(real_target, line)
        line = _TD_CMPLX_KIND_E0.sub(complex_target, line)
    if promote_k8:
        line = _TD_REAL_KIND_D0.sub(real_target, line)
        line = _TD_CMPLX_KIND_D0.sub(complex_target, line)
    # ``REAL(KIND=WP)`` / ``REAL(WP)`` / ``COMPLEX(KIND=WP)`` etc.
    # appear in newer-style LAPACK files (e.g. DGEDMD, DGEDMDQ). The
    # ``wp`` parameter declaration is independently stripped earlier
    # in the pipeline, so the bare ``KIND=WP`` reference would dangle.
    # WP is by convention the half's working precision, so always
    # promote regardless of source_kind.
    line = _TD_REAL_KIND_WP.sub(real_target, line)
    line = _TD_CMPLX_KIND_WP.sub(complex_target, line)
    # Explicit numeric kinds on working-precision-like types. MUMPS's
    # z-half uses ``COMPLEX(kind=8) A(LA)`` where its s/c/d siblings
    # say ``REAL``/``COMPLEX``/``DOUBLE PRECISION``; rewrite all four
    # spellings so they land on the same target type.
    #
    # Restricted to declaration context via the trailing lookahead:
    # the match must be followed by whitespace+identifier (``REAL(4) X``
    # / ``COMPLEX(8) A(LA)``), ``::`` (modern decl), or ``,`` (attribute
    # list). This excludes expression-context ``real(4)`` / ``real(8)``
    # intrinsic calls (e.g. ``real(4)*real(KMAX)``) which have the same
    # token shape but take the integer as an argument, not a kind.
    if promote_k4:
        line = _TD_REAL_KIND_4.sub(real_target, line)
        line = _TD_CMPLX_KIND_4.sub(complex_target, line)
    if promote_k8:
        line = _TD_REAL_KIND_8.sub(real_target, line)
        line = _TD_CMPLX_KIND_8.sub(complex_target, line)

    if not target_mode.is_kind_based:
        line = _filter_known_constants_from_decl(
            line, target_mode, complex_names=complex_names,
        )

    return line


_DECL_START_RE = re.compile(
    r'^(\s+)('
    r'DOUBLE\s+PRECISION|DOUBLE\s+COMPLEX|'
    r'COMPLEX\s*\*\s*16|COMPLEX\s*\*\s*8|'
    r'REAL\s*\*\s*8|REAL\s*\*\s*4|'
    r'TYPE\s*\([^)]+\)|'
    # Bare REAL / COMPLEX without ``*N`` / ``(KIND=...)`` decoration.
    # The negative lookahead rules out ``REAL*N`` / ``REAL(KIND=...)``.
    # ``REAL FUNCTION FOO(...)`` is rejected later by a FUNCTION-bail
    # check on the variable list.
    r'REAL(?!\s*[*(])|COMPLEX(?!\s*[*(])'
    r')(\s*(?:,\s*[A-Za-z][\w]*\s*(?:\([^)]*\))?\s*)*::\s*|\s+)(.+)$',
    re.IGNORECASE,
)


def fix_misdeclared_statement_functions(source: str,
                                          source_kind: int | None = None) -> str:
    """Correct the declared type of statement functions whose body is
    a real-valued expression.

    LAPACK's ``CABS1`` is the textbook example: it is declared
    ``COMPLEX*16`` in ``zla_lin_berr.f`` but its body is
    ``ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) )`` — a real
    expression. F77 tolerated the mismatch via implicit complex←real
    promotion at assignment time, but the post-migration derived-type
    version refuses to assign a ``float64x2`` RHS to a
    ``complex128x2`` LHS inside a statement-function definition.

    This pass scans for any statement function whose RHS is
    syntactically a chain of ``ABS(...)`` terms joined by ``+``/``-``
    — which always evaluates to real — and, if the corresponding
    ``COMPLEX*16`` / ``DOUBLE COMPLEX`` declaration of the function
    name exists in the same file, rewrites it to ``DOUBLE PRECISION``.
    """
    lines = source.splitlines(keepends=True)
    # First pass: find statement function names whose body is real.
    stmt_fn_re = re.compile(
        r'^\s+([A-Za-z_]\w*)\s*\(\s*[A-Za-z_]\w*\s*\)\s*=\s*'
        r'(?:[-+]?\s*ABS\s*\([^()]*(?:\([^()]*\)[^()]*)*\)\s*[-+]?\s*)+$',
        re.IGNORECASE,
    )
    real_names: set[str] = set()
    for raw in lines:
        if raw and raw[0] in ('C', 'c', '*', '!'):
            continue
        body = raw.rstrip()
        m = stmt_fn_re.match(body)
        if m:
            real_names.add(m.group(1).upper())
    if not real_names:
        return source

    # Second pass: demote any single-variable ``COMPLEX*16`` /
    # ``DOUBLE COMPLEX`` declaration of one of those names to
    # ``DOUBLE PRECISION``. Only rewrite declarations whose variable
    # list is exactly a single unadorned identifier — a compound
    # declaration like ``COMPLEX*16  CABS1, CDUM`` would require
    # splitting the line and is left alone.
    cplx_decl_re = re.compile(
        r'^(\s+)(DOUBLE\s+COMPLEX|COMPLEX\s*\*\s*16|COMPLEX\s*\*\s*8|COMPLEX)'
        r'(\s+)([A-Za-z_]\w*)\s*$',
        re.IGNORECASE,
    )
    # Pick a real kind that matches the source half so the rewrite
    # doesn't introduce a kind8 reference into a kind4 source (or
    # vice versa) — see rule (a) discussion in
    # ``replace_type_decls``'s docstring.
    new_decl = 'REAL            ' if source_kind == 4 else 'DOUBLE PRECISION'
    out: list[str] = []
    for raw in lines:
        if raw and raw[0] not in ('C', 'c', '*', '!'):
            m = cplx_decl_re.match(raw.rstrip())
            if m and m.group(4).upper() in real_names:
                nl = '\n' if raw.endswith('\n') else ''
                out.append(f'{m.group(1)}{new_decl}{m.group(3)}{m.group(4)}{nl}')
                continue
        out.append(raw)
    return ''.join(out)


def strip_known_constants_from_decls(
    source: str, target_mode: TargetMode,
) -> tuple[str, dict[str, str]]:
    """Whole-source pre-pass: drop known-constant names from multi-line type decls.

    A continuation line like ``$  DU,GAM,GAMSQ,ONE,RGAMSQ,TWO,ZERO`` is
    invisible to per-line filtering. This pass joins continuation lines
    of a type declaration into one logical statement, removes any plain
    identifier whose uppercase form is in ``target_mode.known_constants``,
    and re-emits the declaration. Items containing parens (array specs)
    or ``=`` (initializers) are preserved verbatim. If every item is
    removed, the entire declaration is dropped.

    Returns ``(new_source, removed_renames)`` where ``removed_renames``
    maps each filtered name (uppercase) to its multifloats replacement
    (e.g. ``'ZERO' -> 'MF_ZERO'``). Callers feed this map to
    :func:`replace_known_constants` so that **only** the names that
    were actually filtered get rewritten in the body — names imported
    via ``USE LA_CONSTANTS_MF, ONLY: zero=>dzero`` are left intact
    because they were never in a local declaration.
    """
    if target_mode.is_kind_based:
        return source, {}
    known = {k.upper(): v for k, v in target_mode.known_constants.items()}
    if not known:
        return source, {}

    # A name is "globally ambiguous" if some other procedure in the
    # same file declares it as COMPLEX. Stripping it here would push a
    # global ZERO -> DD_ZERO rename, but ZERO in the complex procedure
    # is locally COMPLEX — the rename would mistype it. Skip stripping
    # for these names; they stay declared in every scope, with a
    # runtime assignment per scope handling the value. (The two scopes'
    # respective declarations are typed correctly by the per-line type
    # rewrite later.)
    file_complex_names = _scan_complex_var_names(source)

    removed: dict[str, str] = {}
    lines = source.splitlines(keepends=True)
    out: list[str] = []
    i = 0
    while i < len(lines):
        raw = lines[i]
        # Don't touch comment lines
        if raw and raw[0] in ('C', 'c', '*', '!'):
            out.append(raw); i += 1; continue
        rstripped = raw.rstrip('\n')
        m = _DECL_START_RE.match(rstripped)
        if not m:
            out.append(raw); i += 1; continue

        indent, type_text, sep, vars_part = m.groups()

        # Bail on function-return decls: 'REAL FUNCTION FOO(...)'.
        if re.match(r'^\s*FUNCTION\b', vars_part, re.IGNORECASE):
            out.append(raw); i += 1; continue

        # ZERO/ONE/etc. that are declared as a COMPLEX type carry
        # complex semantics. Replacing them globally with MF_ONE
        # (a real float64x2 constant) breaks call sites such as
        # CALL UGEMV(..., ONE, ...) where the dummy is complex.
        # Skip stripping for complex declarations: the local var
        # remains, and convert_parameter_stmts emits an assignment
        # ``ONE = complex128x2(MF_ONE, MF_ZERO)`` later.
        if re.search(r'COMPLEX', type_text, re.IGNORECASE):
            out.append(raw); i += 1; continue

        # Collect continuation lines (fixed-form col-6 marker, OR
        # previous logical line ends with '&').
        stmt_lines = [raw]
        j = i + 1
        prev_amp = vars_part.rstrip().endswith('&')
        while j < len(lines):
            nxt = lines[j]
            if is_continuation_line(nxt):
                stmt_lines.append(nxt); j += 1; continue
            if prev_amp:
                stmt_lines.append(nxt)
                prev_amp = nxt.rstrip('\n').rstrip().endswith('&')
                j += 1; continue
            break

        # Build the joined var-list text. Strip continuation markers
        # ('&' free-form, col-6 char fixed-form) and inline comments.
        def _strip_amp(s: str) -> str:
            s = s.rstrip()
            return s[:-1].rstrip() if s.endswith('&') else s

        full = _strip_amp(vars_part)
        for cl in stmt_lines[1:]:
            body = cl.rstrip('\n')
            if is_continuation_line(body):
                body = body[6:]
            body = body.lstrip()
            if body.startswith('&'):
                body = body[1:]
            full = full + ' ' + _strip_amp(body)

        comment = ''
        bang = full.find('!')
        if bang >= 0:
            comment = full[bang:]
            full = full[:bang]

        # Top-level comma split, respecting parentheses.
        items: list[str] = []
        cur, depth = '', 0
        for ch in full:
            if ch == '(':
                depth += 1; cur += ch
            elif ch == ')':
                depth -= 1; cur += ch
            elif ch == ',' and depth == 0:
                if cur.strip(): items.append(cur.strip())
                cur = ''
            else:
                cur += ch
        if cur.strip():
            items.append(cur.strip())

        # If any item has an '=' initializer, bail out — preserving
        # original is safer than rewriting initializer expressions.
        if any('=' in it for it in items):
            for sl in stmt_lines: out.append(sl)
            i = j; continue

        kept: list[str] = []
        for it in items:
            nm = re.match(r'^([A-Za-z_]\w*)', it)
            if (nm and nm.group(1).upper() in known and it == nm.group(1)
                    and nm.group(1).upper() not in file_complex_names):
                removed[nm.group(1).upper()] = known[nm.group(1).upper()]
                continue  # drop bare known-constant name
            kept.append(it)

        if len(kept) == len(items):
            for sl in stmt_lines: out.append(sl)
            i = j; continue

        if not kept:
            # Whole declaration removed
            i = j; continue

        body = ', '.join(kept)
        rebuilt = f'{indent}{type_text}{sep}{body}'
        if comment:
            rebuilt = rebuilt + ' ' + comment
        nl = '\n' if stmt_lines[0].endswith('\n') else ''
        out.append(rebuilt + nl)
        i = j
    return ''.join(out), removed


# Hoisted patterns reused across hot-path callers (replace_literals,
# replace_generic_conversions, _rewrite_int_of_complex / _rewrite_int_kind_on_real64x2,
# _filter_known_constants_from_decl, _unwrap_redundant_constructors). All
# depend only on constants known at module load time or values that are
# stable for the duration of a run; rebuilding them per line was the
# dominant cost in microbenchmarks of the migrator hot loops.
#
# ``.EQ.`` / ``.AND.`` etc. operator detector used in ``replace_literals``
# to mask operator dots out before the bare-literal pass scans for
# floating-point suffixes.
_FORTRAN_OP_RE = re.compile(
    r'\.\s*(EQ|NE|LT|GT|LE|GE|AND|OR|NOT|TRUE|FALSE|EQV|NEQV)\s*\.',
    re.IGNORECASE,
)
# ``REAL(`` / ``CMPLX(`` paren-call detector for ``replace_generic_conversions``.
_GENERIC_CONV_RE: dict[str, re.Pattern] = {
    _name: re.compile(rf'(?<=[=+\-*/,(.\0])\s*\b({_name})\s*\(', re.IGNORECASE)
    for _name in ('REAL', 'CMPLX')
}
# ``INT(`` / ``NINT(`` detectors for the two _rewrite_int_* helpers.
_INT_CALL_RE = re.compile(r'\bINT\s*\(', re.IGNORECASE)
_NINT_CALL_RE = re.compile(r'\bNINT\s*\(', re.IGNORECASE)


@functools.cache
def _filter_known_decl_re(real_target: str, complex_target: str) -> re.Pattern:
    """Cache the ``_filter_known_constants_from_decl`` matcher per
    (real_target, complex_target) pair. Both come from the run's
    TargetMode and are stable across all per-line invocations."""
    type_alt = f'(?:{re.escape(real_target)}|{re.escape(complex_target)})'
    return re.compile(
        rf'^(\s*)({type_alt})(\s*(?:::)?\s*)(.+?)(\s*(?:!.*)?)$',
        re.IGNORECASE,
    )


@functools.cache
def _unwrap_ctor_re(ctor: str) -> re.Pattern:
    """Cache the ``_unwrap_redundant_constructors`` matcher per
    constructor name (``target_mode.real_constructor`` is run-stable)."""
    return re.compile(rf'\b{re.escape(ctor)}\s*\(\s*([A-Za-z_]\w*)\s*\)')


@functools.cache
def _known_constants_pattern(keys: frozenset) -> re.Pattern:
    """Cache the comma-name alternation in ``replace_known_constants``.

    The ``renames`` dict is invariant for the duration of a file
    (built once before the per-line loop), so the compiled pattern
    can be cached. The function is called per source line; rebuilding
    the alternation on every call was the dominant cost.
    """
    names_alt = '|'.join(re.escape(n) for n in sorted(keys, key=len, reverse=True))
    return re.compile(rf'(?<![A-Za-z0-9_])({names_alt})(?![A-Za-z0-9_])', re.IGNORECASE)


def _filter_known_constants_from_decl(
    line: str,
    target_mode: TargetMode,
    complex_names: set[str] | None = None,
) -> str:
    """Drop known-constant names from a TYPE() declaration's variable list.

    Matches `TYPE(...)` (followed by optional ``::``) followed by a
    comma-separated list of *plain* variable names — not array specs
    like ``A(LDA,*)`` or initializers, which would require a deeper
    parse. Names matching ``target_mode.known_constants`` are removed;
    if every name is removed, the entire declaration is dropped.

    When ``complex_names`` is provided, names in that set are NOT
    removed even from a real-typed declaration — they're declared as
    COMPLEX somewhere else in the file and the multifloats real
    constant would mistype the complex scope.
    """
    real_target = target_mode.real_type
    complex_target = target_mode.complex_type
    known = {k.upper() for k in target_mode.known_constants}
    if complex_names:
        known = known - {n.upper() for n in complex_names}

    m = _filter_known_decl_re(real_target, complex_target).match(line)
    if not m:
        return line
    indent, type_text, sep, vars_part, trailer = m.groups()
    # Skip stripping when the declared type is COMPLEX. ZERO/ONE/etc.
    # declared as complex carry complex semantics; replacing them
    # globally with the real ``MF_ONE``/``MF_ZERO`` would break call
    # sites that pass them as a complex argument (e.g.
    # ``CALL UGEMV(..., ONE, ...)``).
    if complex_target in type_text or 'COMPLEX' in type_text.upper():
        return line
    # Only handle simple comma-separated identifier lists. If any
    # entry contains parens (array spec) or '=' (initializer), bail
    # out — those need a smarter parser and the wholesale rewrite
    # would be unsafe.
    items = [v.strip() for v in vars_part.split(',')]
    if not items or any(('(' in it or '=' in it) for it in items):
        return line
    if not all(re.fullmatch(r'[A-Za-z_]\w*', it) for it in items):
        return line
    kept = [it for it in items if it.upper() not in known]
    if not kept:
        return ''  # entire declaration removed
    if len(kept) == len(items):
        return line
    return f'{indent}{type_text}{sep}{",".join(kept)}{trailer}'


def replace_standalone_real_complex(line: str, target_mode: TargetMode,
                                     source_kind: int | None = None) -> str:
    """Replace standalone REAL/COMPLEX keywords in declaration context.

    Matches both classic F77 syntax (``REAL X, Y``) and modern F90
    attribute-list syntax (``REAL :: X``, ``REAL, POINTER :: X``,
    ``REAL, DIMENSION(:) :: X``). The trailing lookahead requires
    one of:
      ``\\s+[A-Za-z]``  — space + letter (F77 ``REAL X``)
      ``\\s*::``        — F90 ``REAL :: X``
      ``\\s*,``         — F90 ``REAL, attr :: X``

    The leading negative lookahead rejects:
      ``\\s*\\(KIND``   — explicit kind spec like ``REAL(KIND=8)``
      (for COMPLEX, also ``\\*`` and ``\\(`` — ``COMPLEX*16`` and
      ``COMPLEX(KIND=...)`` / ``COMPLEX(x)`` function call).

    Bare REAL / COMPLEX are kind4 by Fortran default. Skip the rewrite
    when ``source_kind == 8`` so a kind8 source half preserves its
    rare bare-REAL declarations as kind4 (rule a).
    """
    if source_kind == 8:
        return line
    real_target = target_mode.real_type
    complex_target = target_mode.complex_type

    # The ``\s*,`` alternative for F90 attribute-list syntax
    # (``REAL, POINTER :: X``) requires ``::`` to appear later on the
    # line to distinguish a type declaration from an intrinsic list
    # like ``INTRINSIC REAL, AIMAG`` (no ``::`` on those).
    _decl_tail = r'(?=\s+[A-Za-z]|\s*::|\s*,[^\n]*::)'
    line = re.sub(
        r'\bREAL\b(?!\s*\(KIND)' + _decl_tail,
        real_target, line, flags=re.IGNORECASE
    )
    line = re.sub(
        r'\bCOMPLEX\b(?!\s*[\*(])' + _decl_tail,
        complex_target, line, flags=re.IGNORECASE
    )
    return line


# ---------------------------------------------------------------------------
# Literal constant replacement
# ---------------------------------------------------------------------------

def replace_literals(line: str, target_mode: TargetMode,
                     source_kind: int | None = None) -> str:
    """Replace floating-point literals with target form.

    KIND mode: 1.0D+0 → 1.0E+0_k, 0.0E+0 → 0.0E+0_k, 0.0 → 0.0E0_k.
    Constructor mode: 1.0D+0 → float64x2('1.0D+0') or float64x2(1.0D0).

    Bare unsuffixed literals (no D/E exponent, no ``_kind``) are also
    promoted so that complex constants like ``(0.0,0.0)`` in C sources
    converge with ``(0.0d0,0.0d0)`` in Z sources after migration.

    Rule (a) gating by ``source_kind`` (mirrors :func:`replace_type_decls`):
    a kind4 source half (S/C) only promotes kind4-shaped literals
    (``1.0E0``, bare ``1.0``); a kind8 source half (D/Z) only promotes
    kind8-shaped literals (``1.0D0``). Cross-kind literals inside a
    half are intentional explicit-kind references and must survive
    migration. ``None`` falls back to the legacy promote-everything
    behavior.
    """
    # Keep-kind preserves the LHS type (DOUBLE PRECISION) or call-site
    # wrapper (dble/dcmplx) on this line. In constructor mode, wrapping
    # an RHS literal as ``real64x2(limbs=[...])`` produces a derived-
    # type value the keep-kind LHS / wrapper cannot absorb (no implicit
    # narrowing exists for derived types). Skip the literal pass on
    # such lines and leave the original DP literal in place — it
    # matches the preserved LHS / wrapper type. The kind_suffix path
    # is unaffected because REAL(KIND=N) → REAL(8) narrowing IS legal.
    if (target_mode.literal_mode == 'constructor'
            and (_KK_SENTINEL in line
                 or _KK_DBLE_SENTINEL in line
                 or _KK_DCMPLX_SENTINEL in line)):
        return line

    # ``INTEGER, PARAMETER :: name = <fp_literal>`` lines exist in
    # upstream MUMPS as a degenerate idiom (e.g. dsol_fwd_aux.F:1093
    # ``INTEGER, PARAMETER :: ZERO = 0.0D0``) — gfortran tolerates the
    # FP literal via integer coercion. In multifloats mode, wrapping
    # the literal to ``real64x2(limbs=[...])`` produces a derived-type
    # value that the INTEGER PARAMETER cannot accept. Leave the
    # original FP literal alone on those lines.
    if (target_mode.literal_mode == 'constructor'
            and re.match(r'^\s*INTEGER\b.*\bPARAMETER\b.*::',
                         line, re.IGNORECASE)):
        return line

    def literal_sub(m):
        mantissa = m.group(1)
        suffix = m.group(2)
        exp_rest = m.group(3)
        if exp_rest.startswith('+'):
            exp_rest = exp_rest[1:]

        # Rule (a): the suffix encodes the literal's source kind
        # (E = kind4, D = kind8). Skip the rewrite if the source half
        # doesn't own that kind.
        suffix_kind = 8 if suffix.upper() == 'D' else 4
        if source_kind is not None and source_kind != suffix_kind:
            return m.group(0)

        if target_mode.literal_mode == 'kind_suffix':
            return f'{mantissa}E{exp_rest}_{target_mode.kind_suffix}'
        else:
            # Named-component (``limbs=...``) structure-constructor
            # form. This is a constant expression and is therefore
            # legal in PARAMETER initializers — the alternative
            # ``float64x2('1.0D0')`` is a function call (the
            # ``mf_from_char`` generic) and would be rejected by the
            # compiler in PARAMETER context.
            return (
                f"{target_mode.real_constructor}"
                f"(limbs=[{mantissa}D{exp_rest}, 0.0_8])"
            )

    parts = re.split(r"('(?:[^']|'')*'|\"(?:[^\"]|\"\")*\")", line)
    for idx in range(0, len(parts), 2):
        seg = parts[idx]
        masked = _FORTRAN_OP_RE.sub(
            lambda m: '\x00' + m.group(1) + '\x00', seg,
        )

        # ``1.23_wp`` style literals are normalized first into a
        # placeholder so that the subsequent passes (D/E exponent and
        # bare-literal substitution) do not see the substituted form
        # and re-wrap its inner ``1.23D0`` text. We restore the
        # placeholders at the very end of the loop.
        wp_placeholders: list[str] = []

        def _wp_sub(m):
            val = m.group(1)
            if 'D' not in val.upper() and 'E' not in val.upper():
                if '.' in val:
                    val += 'D0'
                else:
                    val += '.0D0'
            tok = f"\x01WP{len(wp_placeholders)}\x01"
            wp_placeholders.append(
                f"{target_mode.real_constructor}(limbs=[{val}, 0.0_8])"
            )
            return tok

        if target_mode.literal_mode == 'constructor':
            masked = re.sub(r'(\d+\.\d*|\d*\.\d+)_wp', _wp_sub, masked, flags=re.IGNORECASE)

        # The negative lookbehind also rejects ``[`` (so we don't
        # match the inner ``.0D0`` of a previously-wrapped
        # ``float64x2(limbs=[1.0D0, 0.0_8])`` form) and rejects a
        # leading digit (so we don't match the ``.0D0`` substring of
        # a complete ``1.0D0`` literal that's already been consumed).
        masked = re.sub(
            r'(?<![\[\d])(\d+\.\d*|\d*\.\d+)([DEde])([+-]?\d+)',
            literal_sub, masked,
        )
        # Bare unsuffixed literals (e.g. ``1.0``, ``0.5``) are kind4
        # by Fortran default. Skip when ``source_kind == 8`` so a
        # kind8 source half preserves them as kind4 (rule a).
        if source_kind != 8:
            if target_mode.literal_mode == 'kind_suffix':
                masked = re.sub(
                    r'(?<![.\w])(\d+\.\d*|\d*\.\d+)(?![DdEe\w]|_\d)',
                    rf'\1E0_{target_mode.kind_suffix}', masked,
                )
            else:
                def bare_sub(m):
                    val = m.group(1)
                    return (
                        f"{target_mode.real_constructor}"
                        f"(limbs=[{val}D0, 0.0_8])"
                    )
                # The negative lookbehind on ``[`` skips literals that
                # already live inside a previously-wrapped
                # ``float64x2(limbs=[...])`` form.
                masked = re.sub(
                    r'(?<![.\w\[])(\d+\.\d*|\d*\.\d+)(?![DdEe\w]|_\d)',
                    bare_sub, masked,
                )
            
        if target_mode.literal_mode == 'constructor':
            # Wrap complex constants (float64x2(...), float64x2(...)) in complex128x2(...).
            # An optional unary +/- is allowed before each component to
            # cover cases like ``(-1.0D0, 0.0D0)`` from LAPACK.
            # Bare-name complex literals like ``(MF_ZERO, MF_ZERO)`` are
            # wrapped later in the per-line pipeline by
            # _wrap_bare_complex_literals (after replace_known_constants
            # has produced the MF_* names).
            masked = re.sub(
                r'\(\s*([-+]?\s*' + target_mode.real_constructor + r'\([^)]+\))\s*,\s*([-+]?\s*' + target_mode.real_constructor + r'\([^)]+\))\s*\)',
                rf"{target_mode.complex_constructor}(\1,\2)",
                masked,
                flags=re.IGNORECASE
            )

        # Restore the ``_wp`` literal placeholders.
        if wp_placeholders:
            def _wp_restore(m):
                return wp_placeholders[int(m.group(1))]
            masked = re.sub(r'\x01WP(\d+)\x01', _wp_restore, masked)

        parts[idx] = re.sub(
            r'\x00([A-Za-z]+)\x00', r'.\1.', masked,
        )
    return ''.join(parts)


# ---------------------------------------------------------------------------
# Intrinsic function replacement
# ---------------------------------------------------------------------------

# Precompiled once per process — INTRINSIC_MAP keys are static so the
# per-name patterns never vary. ``replace_intrinsic_calls`` is on the
# hot path (called per statement from migrate_*_form), so rebuilding
# 74 patterns per call was wasteful.
_INTRINSIC_CALL_RE: dict[str, re.Pattern] = {
    _name: re.compile(rf'\b{_name}\s*\(', re.IGNORECASE)
    for _name in INTRINSIC_MAP
}
_INTRINSIC_CALL_RE_REPL: dict[str, re.Pattern] = {
    _name: re.compile(rf'\b({_name})(\s*\()', re.IGNORECASE)
    for _name in INTRINSIC_MAP
}
# Cheap early-out gate: if no intrinsic name from the map appears
# followed by an opening paren on this line, the 74-iteration body
# below has nothing to do. One alternation search per line is much
# cheaper than 74 single-name searches.
_INTRINSIC_CALL_GATE_RE = re.compile(
    r'\b(?:' + '|'.join(INTRINSIC_MAP) + r')\s*\(',
    re.IGNORECASE,
)


def replace_intrinsic_calls(
    line: str,
    target_mode: TargetMode,
    real_names: set[str] | None = None,
    complex_names: set[str] | None = None,
) -> str:
    """Replace type-specific intrinsic function calls."""
    if not _INTRINSIC_CALL_GATE_RE.search(line):
        return line
    for old_name, (new_name, needs_kind) in INTRINSIC_MAP.items():
        pattern = _INTRINSIC_CALL_RE[old_name]
        if needs_kind:
            search_start = 0
            while True:
                m = pattern.search(line, search_start)
                if not m:
                    break
                start = m.start()
                paren_start = line.index('(', m.start())
                depth, pos = 1, paren_start + 1
                while pos < len(line) and depth > 0:
                    if line[pos] == '(': depth += 1
                    elif line[pos] == ')': depth -= 1
                    pos += 1
                if depth == 0:
                    close_pos = pos - 1
                    inner = line[paren_start + 1:close_pos]
                    inner_stripped = inner.strip().upper()
                    old_upper = old_name.upper()
                    is_type_spec_name = old_upper in ('REAL', 'CMPLX', 'COMPLEX')
                    # Bare integer literals are skipped only for the
                    # ambiguous names — ``REAL(3)`` might be the type
                    # spec ``REAL(KIND=3)``, so we leave it alone.
                    # ``DBLE(3)``/``DCMPLX(3)``/etc. are unambiguously
                    # conversion-function calls and must be rewritten
                    # (otherwise an `s*` sibling that wrote ``real(3)``
                    # diverges post-migration: `real(3)` strips to `3`
                    # in the light-diff normalizer; `dble(3)` does not).
                    if (re.match(r'KIND\s*=', inner_stripped)
                            or (is_type_spec_name and inner_stripped.isdigit())
                            or (is_type_spec_name and re.match(r'^[A-Z_]\w*$', inner_stripped))):
                        search_start = pos
                        continue
                    
                    if target_mode.intrinsic_mode == 'add_kind':
                        depth_k = 0
                        has_top_kind = False
                        ii = 0
                        while ii < len(inner):
                            ch = inner[ii]
                            if ch == '(': depth_k += 1
                            elif ch == ')': depth_k -= 1
                            elif (depth_k == 0 and inner[ii:ii + 4].upper() == 'KIND'
                                    and re.match(r'KIND\s*=', inner[ii:], re.IGNORECASE)):
                                has_top_kind = True
                                break
                            ii += 1
                        if has_top_kind:
                            search_start = pos
                            continue
                        replacement = f'{new_name}({inner}, KIND={target_mode.kind_suffix})'
                    else:
                        # wrap_constructor mode (multifloats).
                        # ``DBLE`` / ``DREAL`` / ``REAL`` map to the
                        # ``float64x2(...)`` generic constructor, which
                        # handles every input type uniformly: integer,
                        # real(sp/dp), float64x2 (identity), and complex
                        # (extract real part). Older code routed these
                        # through ``MF_REAL`` but that interface only
                        # accepts float64x2 / complex128x2, breaking
                        # ``DBLE(integer)`` from LAPACK.
                        # ``CMPLX`` / ``DCMPLX`` map to ``complex128x2``.
                        # multifloats's complex128x2 interface only has
                        # overloads for (float64x2[, float64x2]) and
                        # (real(dp)[, real(dp)]) — there's no integer
                        # overload, so a 1-arg call with an integer
                        # argument falls back to the structure
                        # constructor and fails. Pre-wrap any single-arg
                        # cmplx with float64x2(...) so it always picks
                        # the cx_from_mf_1 procedure.
                        # ``AIMAG`` / ``CONJG`` are module generics —
                        # leave the name alone.
                        if old_name.upper() in ('REAL', 'DREAL', 'DBLE'):
                            # Drop the optional kind-spec second argument
                            # when present: `REAL(x, wp)` / `REAL(x, KIND=wp)`
                            # / `DBLE(x, wp)` becomes `real64x2(x)` for
                            # multifloats. The constructor takes a single
                            # value argument; passing the kind-symbol as a
                            # second positional component overflows the
                            # 1-component derived type ("Too many components
                            # in structure constructor"). Split at the
                            # top-level comma respecting nested parens.
                            depth_k = 0
                            split_at = -1
                            for ii, ch in enumerate(inner):
                                if ch == '(': depth_k += 1
                                elif ch == ')': depth_k -= 1
                                elif ch == ',' and depth_k == 0:
                                    split_at = ii
                                    break
                            inner_first = inner[:split_at] if split_at >= 0 else inner
                            replacement = f'{target_mode.real_constructor}({inner_first})'
                        elif old_name.upper() in ('CMPLX', 'DCMPLX'):
                            wrapped = _wrap_complex_args(
                                inner, target_mode, real_names,
                                complex_names=complex_names,
                            )
                            replacement = f'{target_mode.complex_constructor}({wrapped})'
                        else:
                            replacement = f'{new_name}({inner})'

                    line = line[:start] + replacement + line[close_pos + 1:]
                    search_start = start + len(replacement)
                else:
                    if old_name.upper() != new_name.upper():
                        matched = line[start:m.end() - 1]
                        repl = new_name.upper() if matched.isupper() else new_name.lower()
                        line = line[:start] + repl + line[start + len(matched):]
                        search_start = start + len(repl)
                    else:
                        search_start = m.end()
                    continue
        else:
            def _call_replace(m, _new=new_name):
                matched_name = m.group(1)
                rest = m.group(2)
                return (_new.upper() if matched_name.isupper() else _new.lower()) + rest

            line = _INTRINSIC_CALL_RE_REPL[old_name].sub(_call_replace, line)
    return line


def replace_generic_conversions(
    line: str,
    target_mode: TargetMode,
    complex_names: set[str] | None = None,
) -> str:
    """Add KIND (or wrap in constructor) to generic REAL() and CMPLX() calls in expression context."""
    is_fixed_cont = (len(line) >= 6 and line[:5] == '     ' and line[5] not in (' ', '0', '!', 'C', 'c', '*', '\t'))
    cont_marker = line[5] if is_fixed_cont else ''
    if is_fixed_cont:
        line = line[:5] + '\0' + line[6:]

    for name in ('REAL', 'CMPLX'):
        pattern = _GENERIC_CONV_RE[name]
        search_start = 0
        while True:
            m = pattern.search(line, search_start)
            if not m: break
            name_start = m.start(1)
            paren_start = line.index('(', name_start)
            depth, pos = 1, paren_start + 1
            while pos < len(line) and depth > 0:
                if line[pos] == '(': depth += 1
                elif line[pos] == ')': depth -= 1
                pos += 1
            if depth != 0: break
            close_pos = pos - 1
            inner = line[paren_start + 1:close_pos]
            
            if target_mode.intrinsic_mode == 'add_kind':
                if re.search(r'\bKIND\s*=', inner, re.IGNORECASE):
                    search_start = pos
                    continue
                top_commas = 0
                d = 0
                for ch in inner:
                    if ch == '(': d += 1
                    elif ch == ')': d -= 1
                    elif ch == ',' and d == 0: top_commas += 1
                max_args = 1 if name == 'REAL' else 2
                if top_commas >= max_args:
                    search_start = pos
                    continue
                replacement = f'{name}({inner}, KIND={target_mode.kind_suffix})'
            else:
                # multifloats: REAL(x) and CMPLX(...) become the
                # universal generic constructors. float64x2 / complex128x2
                # have overloads for every input type (int, real(sp/dp),
                # float64x2, complex128x2 — extracts real part for the
                # first), which DBLE/REAL/CMPLX in LAPACK rely on.
                # Single-arg CMPLX(int) needs the inner pre-wrapped with
                # float64x2 so it picks the (float64x2)→complex128x2
                # interface procedure instead of the structure
                # constructor (which has no integer overload).
                if name.upper() == 'REAL':
                    # Drop the optional kind-spec second argument
                    # (`REAL(x, wp)` / `REAL(x, KIND=wp)`) — see the
                    # parallel comment in the wrap_constructor branch
                    # of `replace_intrinsic_calls` for the rationale.
                    depth_k = 0
                    split_at = -1
                    for ii, ch in enumerate(inner):
                        if ch == '(': depth_k += 1
                        elif ch == ')': depth_k -= 1
                        elif ch == ',' and depth_k == 0:
                            split_at = ii
                            break
                    inner_first = inner[:split_at] if split_at >= 0 else inner
                    replacement = f'{target_mode.real_constructor}({inner_first})'
                else:
                    wrapped = _wrap_complex_args(
                        inner, target_mode, None,
                        complex_names=complex_names,
                    )
                    replacement = f'{target_mode.complex_constructor}({wrapped})'

            line = line[:name_start] + replacement + line[close_pos + 1:]
            search_start = name_start + len(replacement)
    if is_fixed_cont:
        line = line[:5] + cont_marker + line[6:]
    return line


def replace_intrinsic_decls(line: str, target_mode: TargetMode | None = None) -> str:
    """Replace intrinsic names in INTRINSIC declarations.

    In multifloats mode, additionally drop names that are overloaded by
    the multifloats module — gfortran refuses to accept a USE-associated
    generic name when an INTRINSIC declaration in the same scope binds
    the name to the standard intrinsic of incompatible signature.
    """
    if not re.match(r'\s+INTRINSIC\b', line, re.IGNORECASE):
        return line
    for old_name, new_name in INTRINSIC_DECL_MAP.items():
        line = re.sub(rf'\b{old_name}\b', new_name, line, flags=re.IGNORECASE)

    m = re.match(r'(\s+INTRINSIC\s+)(.*)', line, re.IGNORECASE)
    if m:
        prefix, name_list = m.group(1), m.group(2)
        newline = '\n' if line.endswith('\n') else ''
        stripped = name_list.rstrip().rstrip('\n')
        trail = ''
        if stripped.endswith('&'):
            trail = ' &'
            stripped = stripped[:-1]
        elif stripped.endswith(','):
            trail = ','
            stripped = stripped[:-1]
        # Strip an optional ``::`` separator after the keyword.
        sep_m = re.match(r'\s*::\s*', stripped)
        sep = ''
        if sep_m:
            sep = stripped[:sep_m.end()]
            stripped = stripped[sep_m.end():]
        names = [n.strip() for n in stripped.split(',') if n.strip()]
        seen: set[str] = set()
        deduped: list[str] = []
        for n in names:
            key = n.upper()
            if key not in seen:
                seen.add(key)
                deduped.append(n)
        if target_mode is not None and target_mode.intrinsic_mode == 'wrap_constructor':
            overloaded = frozenset(n.upper() for n in target_mode.module_generic_names)
            deduped = [n for n in deduped if n.upper() not in overloaded]
        if not deduped:
            # Whole declaration empty — drop the line so we don't emit a
            # bare ``INTRINSIC`` keyword. The trailing comment, if any,
            # is preserved as a regular comment line.
            return ''
        line = prefix + sep + ', '.join(deduped) + trail + newline
    return line



def _dedup_intrinsic_stmts(text: str, target_mode: TargetMode | None = None) -> str:
    """Remove duplicate names from multi-line INTRINSIC statements.

    In multifloats mode, also drop names that are overloaded by the
    multifloats module — they conflict with the use-associated generics.
    """
    lines = text.split('\n')
    result: list[str] = []
    i = 0
    while i < len(lines):
        line = lines[i]
        if not re.match(r'\s+INTRINSIC\b', line, re.IGNORECASE):
            result.append(line)
            i += 1
            continue

        stmt_lines = [line]
        j = i + 1
        while j < len(lines):
            next_line = lines[j]
            is_fixed_cont = (len(next_line) > 6 and next_line[:5] == '     ' and next_line[5] not in (' ', '0') and next_line[0] not in ('C', 'c', '*', '!'))
            if is_fixed_cont:
                stmt_lines.append(next_line)
                j += 1
            elif stmt_lines[-1].rstrip().endswith('&'):
                stmt_lines.append(next_line)
                j += 1
            else:
                break

        m = re.match(r'(\s+INTRINSIC\s+)', stmt_lines[0], re.IGNORECASE)
        if not m:
            result.extend(stmt_lines)
            i = j
            continue

        prefix = m.group(1)
        name_part = stmt_lines[0][len(prefix):]
        for sl in stmt_lines[1:]:
            if len(sl) > 5 and sl[0] == ' ' and sl[5] not in (' ', '0', ''):
                name_part += ' ' + sl[6:]
            else:
                stripped = sl.lstrip()
                if stripped.startswith('&'): stripped = stripped[1:]
                name_part += ' ' + stripped

        name_part = name_part.replace('&', ' ')
        # Drop an optional ``::`` separator following the keyword.
        sep_m = re.match(r'\s*::\s*', name_part)
        if sep_m:
            name_part = name_part[sep_m.end():]
        names = [n.strip() for n in name_part.split(',') if n.strip()]
        seen: set[str] = set()
        deduped: list[str] = []
        for n in names:
            key = n.upper()
            if key not in seen:
                seen.add(key)
                deduped.append(n)
        if target_mode is not None and target_mode.intrinsic_mode == 'wrap_constructor':
            overloaded = frozenset(n.upper() for n in target_mode.module_generic_names)
            deduped = [n for n in deduped if n.upper() not in overloaded]

        if not deduped:
            # Whole INTRINSIC statement empty: drop the entire (possibly
            # multi-line) declaration.
            i = j
            continue

        body = ', '.join(deduped)
        full = prefix + body
        if len(full) <= 72:
            result.append(full)
        else:
            cont_prefix = '     $                   '
            first_cap, cont_cap = 72 - len(prefix) - 1, 72 - len(cont_prefix) - 1
            chunks: list[str] = []
            cur = ''
            for name in deduped:
                addition = (', ' + name) if cur else name
                cap = first_cap if not chunks else cont_cap
                if not cur: cur = name
                elif len(cur) + len(addition) <= cap: cur += addition
                else:
                    chunks.append(cur)
                    cur = name
            if cur: chunks.append(cur)
            for ci, chunk in enumerate(chunks):
                trail = ',' if ci < len(chunks) - 1 else ''
                result.append((prefix if ci == 0 else cont_prefix) + chunk + trail)
        i = j
    return '\n'.join(result)


_ROUTINE_NAME_TOK_RE = re.compile(r'[A-Za-z_]\w*')


def replace_routine_names(line: str, rename_map: dict[str, str]) -> str:
    """Replace routine names using the rename map (case-preserving).

    Tokenize-then-lookup: scan ``line`` for identifier-shaped tokens
    and dict-look each one up. This sidesteps the multi-thousand
    alternation regex that Python's ``re`` engine evaluates by
    backtracking through every alternative at every position, which
    profiling showed dominated migration runtime (~60% of total time
    for LAPACK kind16). One small DFA-friendly tokenizer plus O(1)
    dict lookups is dramatically faster (70× on a representative
    LAPACK rename_map) and semantically equivalent: identifier tokens
    captured by ``[A-Za-z_]\\w*`` match exactly the ``\\b...\\b``-
    bounded forms the alternation regex used to find.
    """
    _, upper_map = _get_rename_pattern(rename_map)
    if not upper_map:
        return line

    def case_replace(m):
        tok = m.group(0)
        new = upper_map.get(tok.upper())
        if new is None:
            return tok
        if tok.isupper():
            return new.upper()
        if tok.islower():
            return new.lower()
        return new.capitalize()

    return _ROUTINE_NAME_TOK_RE.sub(case_replace, line)


_INCLUDE_RE = re.compile(
    r'''^(?P<lead>\s*(?:INCLUDE|include|Include))\s+(?P<q>['"])(?P<name>[^'"]+)(?P=q)(?P<tail>.*)$''',
)


def replace_include_filenames(line: str, rename_map: dict[str, str]) -> str:
    """Rewrite ``INCLUDE 'xxx.h'`` when ``xxx`` is in the rename_map.

    Needed for precision-prefixed headers (e.g. MUMPS's ``dmumps_struc.h``,
    which the ``mumps_struc`` recipe migrates to ``ddmumps_struc.h`` in
    multifloats mode). The including ``.F`` file's literal filename
    string isn't touched by :func:`replace_routine_names` because it's
    inside a quoted Fortran string, so it needs a dedicated rewrite.
    """
    # Fast reject: ``INCLUDE`` statements are rare (≈1 per file) so we
    # avoid the per-call rename_map scan until we know the line matches.
    m = _INCLUDE_RE.match(line)
    if not m:
        return line
    name = m.group('name')
    stem = Path(name).stem
    ext = Path(name).suffix
    upper_stem = stem.upper()
    # rename_map keys are already uppercase (built by build_rename_map);
    # fall back to a case-insensitive lookup if not.
    new_stem = rename_map.get(upper_stem) or rename_map.get(stem)
    if new_stem is None:
        return line
    new_name = (new_stem.upper() if stem.isupper() else
                (new_stem.lower() if stem.islower() else new_stem)) + ext
    return f'{m.group("lead")} {m.group("q")}{new_name}{m.group("q")}{m.group("tail")}'

_RENAME_PATTERN_CACHE: dict[tuple[int, int], tuple[re.Pattern, dict[str, str]]] = {}

def _get_rename_pattern(rename_map: dict[str, str]) -> tuple[re.Pattern, dict[str, str]]:
    # Cache key is (id(rename_map), len(rename_map)). The migrator
    # builds rename_map once per run and never mutates it, so id is
    # stable for the run's duration. The length tiebreaker guards
    # against the (vanishingly unlikely) case where a GC'd dict's id
    # is reused for a same-content map; in practice the caller holds
    # rename_map live, so this is belt-and-suspenders.
    #
    # The previous implementation re-built ``{k.upper() for k in
    # rename_map.keys()}`` on every cache hit to verify the cache —
    # that O(N) scan dominated for MUMPS (~6k renames × thousands of
    # per-line calls per file).
    key = (id(rename_map), len(rename_map))
    cached = _RENAME_PATTERN_CACHE.get(key)
    if cached is not None:
        return cached
    upper_map = {k.upper(): v for k, v in rename_map.items()}
    names = sorted(upper_map.keys(), key=len, reverse=True)
    pattern = re.compile(r'\b(' + '|'.join(re.escape(n) for n in names) + r')\b', re.IGNORECASE) if names else re.compile(r'(?!x)x')
    _RENAME_PATTERN_CACHE[key] = (pattern, upper_map)
    return pattern, upper_map


def _all_known_constant_renames(target_mode: TargetMode) -> dict[str, str]:
    """Combined uppercase rename map for known + la_constants names."""
    out: dict[str, str] = {}
    for k, v in target_mode.known_constants.items():
        out[k.upper()] = v
    for k, v in target_mode.la_constants_map.items():
        out[k.upper()] = v
    return out


_DECL_KEYWORD_RE = re.compile(
    r'^\s*(?:TYPE\s*\(|INTEGER|REAL|COMPLEX|LOGICAL|CHARACTER|'
    r'DOUBLE\s+PRECISION|DOUBLE\s+COMPLEX|PARAMETER\b|DATA\b|'
    r'IMPLICIT\b|DIMENSION\b|EXTERNAL\b|INTRINSIC\b|SAVE\b|'
    r'COMMON\b|EQUIVALENCE\b|USE\b|'
    # Keep-kind sentinels — the pre-pass masks ``DOUBLE PRECISION`` /
    # ``dble`` / ``dcmplx`` on keep-kind lines, and these decl-keyword
    # detection points must treat the masked tokens the same as the
    # un-masked originals so that subsequent passes (e.g.
    # ``replace_known_constants``) don't try to rewrite identifiers in
    # what is logically still a declaration line.
    r'__KEEPKIND_DP__\b)',
    re.IGNORECASE,
)


def replace_known_constants(
    line: str,
    target_mode: TargetMode,
    renames: dict[str, str] | None = None,
) -> str:
    """Substitute names in ``renames`` in the code portion of ``line`` only.

    ``renames`` is the per-file set of names that were removed from a
    local declaration by :func:`strip_known_constants_from_decls`. They
    are exactly the names that need module-import substitution. Names
    bound via a ``USE LA_CONSTANTS_MF, ONLY: zero=>dzero`` alias are
    NOT in this set and therefore are not touched here.

    When ``renames`` is ``None`` the function falls back to the union
    of ``target_mode.known_constants`` and ``target_mode.la_constants_map``;
    that legacy mode is preserved for callers that haven't yet been
    threaded through the new per-file API.

    Skips: comment lines (fixed-form column-1 ``C/c/*/!`` and bare
    ``!``-prefixed free-form), inline comment text after ``!``, and
    declaration / PARAMETER / DATA / USE statement lines (those need
    structural handling, not regex). String literals are also masked
    out so an English ``ZERO`` inside an XERBLA-style message stays
    untouched.
    """
    if target_mode.is_kind_based:
        return line
    if not line or is_comment_line(line):
        return line
    if _DECL_KEYWORD_RE.match(line):
        return line

    if renames is None:
        renames = _all_known_constant_renames(target_mode)
    if not renames:
        return line
    renames = {k.upper(): v for k, v in renames.items()}

    # An assignment of the form ``ZERO = <expr>`` typically came from
    # ``convert_parameter_stmts`` keeping a name as a runtime local
    # because it carries complex semantics (e.g. ``COMPLEX(kind=8) ZERO;
    # PARAMETER(ZERO = (0.0D0, 0.0D0))``). Renaming the LHS to
    # ``DD_ZERO`` (a multifloats public constant) would turn the
    # assignment into a write to an imported PARAMETER, which gfortran
    # rejects as "Named constant in variable definition context". Skip
    # the rename for the LHS name on such lines (use sites on the
    # right-hand side are still a problem at the procedure-scope level
    # but are out of scope for this per-line fix).
    assn_m = re.match(r'^\s*([A-Za-z_]\w*)\s*=(?!=)', line)
    if assn_m and assn_m.group(1).upper() in renames:
        renames = {
            k: v for k, v in renames.items()
            if k != assn_m.group(1).upper()
        }
        if not renames:
            return line

    # Split off any inline comment so we don't substitute inside it.
    # We must respect string literals when locating the '!'.
    code_end = len(line)
    in_str, qch = False, ''
    i = 0
    while i < len(line):
        ch = line[i]
        if in_str:
            if ch == qch:
                if i + 1 < len(line) and line[i + 1] == qch:
                    i += 2  # doubled-quote escape stays inside the string
                    continue
                in_str = False
        elif ch in ("'", '"'):
            in_str, qch = True, ch
        elif ch == '!':
            code_end = i
            break
        i += 1
    code, tail = line[:code_end], line[code_end:]

    # Mask out string literal interiors in code segment.
    parts = re.split(r"('(?:[^']|'')*'|\"(?:[^\"]|\"\")*\")", code)
    pattern = _known_constants_pattern(frozenset(renames.keys()))

    def _sub(m):
        return renames[m.group(1).upper()]

    for idx in range(0, len(parts), 2):
        parts[idx] = pattern.sub(_sub, parts[idx])
    return ''.join(parts) + tail


_COMPLEX_DECL_RE = re.compile(
    r'^\s+(?:DOUBLE\s+COMPLEX|COMPLEX\s*\*\s*(?:8|16)'
    r'|COMPLEX\s*\(\s*(?:KIND\s*=\s*)?\w+\s*\)'
    r'|COMPLEX(?!\s*[*(])'
    r'|TYPE\s*\(\s*cmplx64x2\s*\))',
    re.IGNORECASE,
)
_REAL_DECL_RE = re.compile(
    r'^\s+(?:DOUBLE\s+PRECISION|REAL\s*\*\s*(?:4|8)'
    r'|REAL\s*\(\s*(?:KIND\s*=\s*)?\w+\s*\)'
    r'|REAL(?!\s*[*(])'
    r'|TYPE\s*\(\s*real64x2\s*\))',
    re.IGNORECASE,
)
_INTEGER_DECL_RE = re.compile(
    r'^\s+INTEGER(?!\s*\()',
    re.IGNORECASE,
)


def _scan_typed_var_names(source: str, decl_re: re.Pattern) -> set[str]:
    """Return the set of (uppercase) variable names declared by lines
    matching ``decl_re``. Multi-line decls are joined first.
    """
    lines = source.splitlines()
    out: set[str] = set()
    i = 0
    while i < len(lines):
        line = lines[i]
        if line and line[0] in ('C', 'c', '*', '!'):
            i += 1; continue
        m = decl_re.match(line)
        if not m:
            i += 1; continue
        joined = line
        j = i + 1
        while j < len(lines) and is_continuation_line(lines[j]):
            joined += ' ' + lines[j][6:]
            j += 1
        m = decl_re.match(joined)
        rest = joined[m.end():]
        if '::' in rest:
            rest = rest.split('::', 1)[1]
        if '!' in rest:
            rest = rest.split('!', 1)[0]
        # Top-level comma split, ignoring parens
        items, cur, depth = [], '', 0
        for ch in rest:
            if ch == '(':
                depth += 1; cur += ch
            elif ch == ')':
                depth -= 1; cur += ch
            elif ch == ',' and depth == 0:
                items.append(cur); cur = ''
            else:
                cur += ch
        if cur.strip():
            items.append(cur)
        for it in items:
            nm = re.match(r'\s*([A-Za-z_]\w*)', it)
            if nm:
                out.add(nm.group(1).upper())
        i = j
    return out


def _scan_complex_var_names(source: str) -> set[str]:
    return _scan_typed_var_names(source, _COMPLEX_DECL_RE)


def _scan_real_var_names(source: str) -> set[str]:
    return _scan_typed_var_names(source, _REAL_DECL_RE)


def _scan_integer_var_names(source: str) -> set[str]:
    return _scan_typed_var_names(source, _INTEGER_DECL_RE)


# Note: an experimental ``_force_int_assignment`` pass was prototyped
# here that wraps the RHS of ``INT_VAR = ...`` with ``INT(...)`` when
# the RHS mentions a known float64x2 variable. It was removed because
# the heuristic ("any token in real_names") misclassifies the case
# where a float64x2 variable is *passed* to an integer-returning
# function (e.g. ``JP = J - 1 + IUAMAX(M-J+1, A(J,J), 1)`` where A is
# float64x2 but IUAMAX returns INTEGER). Reliable handling needs
# semantic facts (the migrated function's return type), which is
# Phase 1.5 work.


def _rewrite_int_kind_on_real64x2(
    line: str,
    target_mode: TargetMode,
    real_names: set[str] | None = None,
) -> str:
    """Multifloats only — rewrite ``INT(expr, KIND_INT)`` to
    ``INT(dble(expr), KIND_INT)`` when ``expr`` contains a
    ``real64x2(`` token, or when ``expr`` is a single identifier
    declared as ``TYPE(real64x2)`` in the file (``real_names``).

    Multifloats publishes ``int`` only as ``dd_int(real64x2)`` returning
    default-kind integer — no ``int(real64x2, kind)`` overload exists.
    The migrator's wrap_constructor pass replaces literal real arguments
    with ``real64x2(...)`` calls but leaves the surrounding
    ``INT(..., 8)`` shape intact, leading to "Generic function 'int' is
    not consistent with a specific intrinsic interface" build failures
    in MUMPS files like ``mana_reordertree.F`` / ``mumps_ooc.F``.

    Routing the inner through ``dble`` (multifloats's public generic
    over real64x2 → real(dp), via ``dd_to_dp``) lets the standard
    Fortran INT intrinsic handle the kind selector. Values that flow
    through this pattern (cost counters, matrix size estimates) fit
    comfortably in double precision. ``_build_use_only_clause`` adds
    ``dble`` to the only-list whenever ``INT`` and ``real64x2(``
    co-occur, so the import is in place when this rewrite fires.
    """
    if target_mode.intrinsic_mode != 'wrap_constructor':
        return line
    if not line or line[0] in ('C', 'c', '*', '!'):
        return line
    lower = line.lower()
    # ``int (`` with whitespace before the paren is also a valid call;
    # use the regex below as the authoritative gate, not a substring
    # check that misses the spaced form.
    if not re.search(r'\bint\s*\(', lower):
        return line
    real_names_upper = {n.upper() for n in (real_names or ())}
    if 'real64x2' not in lower and not real_names_upper:
        return line
    pattern = _INT_CALL_RE
    out: list[str] = []
    i = 0
    while i < len(line):
        m = pattern.search(line, i)
        if not m:
            out.append(line[i:])
            break
        out.append(line[i:m.start()])
        paren_open = m.end() - 1
        depth, j = 1, paren_open + 1
        while j < len(line) and depth > 0:
            ch = line[j]
            if ch == '(': depth += 1
            elif ch == ')': depth -= 1
            j += 1
        if depth != 0:
            out.append(line[m.start():])
            break
        inner = line[paren_open + 1:j - 1]
        parts: list[str] = []
        cur, d = '', 0
        for ch in inner:
            if ch == '(': d += 1; cur += ch
            elif ch == ')': d -= 1; cur += ch
            elif ch == ',' and d == 0:
                parts.append(cur); cur = ''
            else:
                cur += ch
        if cur:
            parts.append(cur)
        if (len(parts) == 2
                and re.fullmatch(r'\s*(?:KIND\s*=\s*)?\w+\s*', parts[1])):
            expr = parts[0].strip()
            ident_m = re.fullmatch(r'[A-Za-z_]\w*', expr)
            inner_is_real = (
                'real64x2(' in parts[0].lower()
                or (ident_m is not None
                    and ident_m.group(0).upper() in real_names_upper)
            )
            if inner_is_real:
                head = line[m.start():m.end()]
                kind_arg = parts[1].strip()
                replacement = f'{head}dble({expr}), {kind_arg})'
                out.append(replacement)
            else:
                out.append(line[m.start():j])
        else:
            out.append(line[m.start():j])
        i = j
    return ''.join(out)


def _rewrite_int_of_complex(line: str, complex_names: set[str]) -> str:
    """Wrap the argument of ``INT(...)`` / ``NINT(...)`` with ``MF_REAL``
    when its leading identifier is a known complex variable.

    Multifloats's ``int`` interface only accepts float64x2; calling it
    on complex128x2 fails to dispatch (gfortran does NOT fall back to
    the standard intrinsic for derived-type args, the way it does for
    e.g. integer args). ``MF_REAL`` has overloads for both float64x2
    (identity) and complex128x2 (real-part extraction), so wrapping
    the argument is type-safe.
    """
    if not complex_names or not line:
        return line
    if line[0] in ('C', 'c', '*', '!'):
        return line

    def _process(name_re: re.Pattern, src: str) -> str:
        out = src
        pos = 0
        while True:
            m = name_re.search(out, pos)
            if not m:
                break
            paren_open = m.end() - 1
            depth = 1
            i = paren_open + 1
            while i < len(out) and depth > 0:
                ch = out[i]
                if ch == '(':
                    depth += 1
                elif ch == ')':
                    depth -= 1
                i += 1
            if depth != 0:
                break
            paren_close = i - 1
            inner = out[paren_open + 1:paren_close]
            head = re.match(r'\s*([A-Za-z_]\w*)', inner)
            if head and head.group(1).upper() in complex_names:
                # Use multifloats's overloaded ``real`` generic
                # (cdd_to_dd / dd_to_dd) — it extracts the real
                # component of cmplx64x2 as a real64x2, which the
                # ``int`` / ``nint`` overloads then accept. The earlier
                # ``DD_REAL`` spelling targeted an upstream symbol
                # that no longer exists; ``real`` is in the multifloats
                # generic_names list and gets imported via the USE
                # only-clause whenever a file references it.
                replacement = f'real({inner.strip()})'
                out = out[:paren_open + 1] + replacement + out[paren_close:]
                pos = paren_open + 1 + len(replacement) + 1
            else:
                pos = paren_close + 1
        return out

    line = _process(_INT_CALL_RE, line)
    line = _process(_NINT_CALL_RE, line)
    return line


def _wrap_complex_args(
    inner: str, target_mode: TargetMode, real_names: set[str] | None,
    complex_names: set[str] | None = None,
) -> str:
    """Pre-wrap each top-level argument of a CMPLX-style call so the
    multifloats complex128x2 interface can dispatch.

    multifloats's ``complex128x2`` interface only has overloads for
    (float64x2[, float64x2]) and (real(dp)[, real(dp)]). For LAPACK's
    common patterns ``CMPLX(N)`` and ``CMPLX(N, 0)`` with integer
    arguments, the call falls back to the structure constructor and
    fails. Wrapping each arg with ``float64x2(...)`` redirects to
    ``cx_from_mf_*``.

    Skip the wrap when the arg is already a float64x2 expression
    (bare ``MF_*`` constant, ``float64x2(...)`` call, or a known
    float64x2 local) — wrapping again would fail because float64x2
    has no identity constructor.

    Also skip the wrap when the arg is itself a complex expression
    (any token in ``complex_names``). ``CMPLX(D, KIND=N)`` /
    ``DCMPLX(D)`` on a complex ``D`` is the Fortran identity (preserves
    Re/Im). Wrapping with ``float64x2(...)`` would dispatch to
    ``dd_from_cdd`` which discards the imaginary part — silently
    corrupting the result.
    """
    parts: list[str] = []
    cur, depth = '', 0
    for ch in inner:
        if ch == '(':
            depth += 1; cur += ch
        elif ch == ')':
            depth -= 1; cur += ch
        elif ch == ',' and depth == 0:
            parts.append(cur); cur = ''
        else:
            cur += ch
    if cur:
        parts.append(cur)

    # Drop a trailing ``KIND=...`` argument. The kind selector is only
    # meaningful in the original CMPLX call signature; the multifloats
    # ``complex128x2`` interface has no equivalent.
    parts = [p for p in parts if not re.match(r'\s*KIND\s*=', p, re.IGNORECASE)]

    # Type-conversion intrinsics that replace_generic_conversions() will
    # handle later — don't pre-wrap them or we get double wrapping.
    _CONVERSION_INTRINSICS = frozenset({
        'REAL', 'DBLE', 'SNGL', 'DREAL', 'DFLOAT', 'FLOAT',
        'CMPLX', 'DCMPLX',
    })

    def _wrap_one(arg: str) -> str:
        s = arg.strip()
        if not s:
            return arg
        # Strip a leading unary +/- and surrounding whitespace.
        body = s.lstrip('+-').lstrip()
        head = re.match(r'([A-Za-z_]\w*)', body)
        # Detect any float64x2 token in the expression — if any operand
        # is a known float64x2 / MF_* / float64x2() call, the whole
        # expression is float64x2 (multifloats provides operator
        # overloads for *, /, +, -). Wrapping the whole expression
        # again with float64x2(...) would fail.
        if body.lower().startswith(target_mode.real_constructor.lower() + '('):
            return arg
        # Skip type-conversion intrinsics — they will be replaced by
        # replace_generic_conversions() later in the pipeline.
        if head and head.group(1).upper() in _CONVERSION_INTRINSICS:
            after_name = body[head.end():]
            if after_name.lstrip().startswith('('):
                return arg
        if real_names:
            for tok in re.finditer(r'\b([A-Za-z_]\w*)\b', body):
                u = tok.group(1).upper()
                if u in real_names or u.startswith('MF_'):
                    return arg
        if complex_names:
            for tok in re.finditer(r'\b([A-Za-z_]\w*)\b', body):
                if tok.group(1).upper() in complex_names:
                    return arg
        if head and head.group(1).upper().startswith('MF_'):
            return arg
        return f'{target_mode.real_constructor}({s})'

    return ','.join(_wrap_one(p) for p in parts)


def _unwrap_redundant_constructors(
    line: str, target_mode: TargetMode, real_names: set[str] | None = None,
) -> str:
    """Drop ``float64x2(arg)`` wrappers when ``arg`` is already float64x2.

    The migrator generates ``float64x2(arg)`` for ``DBLE(arg)`` /
    ``REAL(arg)`` calls because the float64x2 interface is the
    universal converter. After the rename pass, ``arg`` may be a bare
    ``MF_*`` constant (which is itself float64x2), or a known
    float64x2 local variable. Wrapping those in ``float64x2(...)``
    fails because multifloats has no identity constructor.
    """
    if target_mode.is_kind_based or not target_mode.real_constructor:
        return line
    ctor = target_mode.real_constructor
    pattern = _unwrap_ctor_re(ctor)

    def _sub(m):
        name = m.group(1)
        upper = name.upper()
        if upper.startswith('MF_'):
            return name
        if real_names and upper in real_names:
            return name
        return m.group(0)

    return pattern.sub(_sub, line)


def _wrap_bare_complex_literals(line: str, target_mode: TargetMode) -> str:
    """Wrap ``(MF_NAME, MF_NAME)`` style complex literals.

    After replace_known_constants has rewritten ``(ZERO, ZERO)`` (etc.)
    to ``(MF_ZERO, MF_ZERO)``, the result is still a Fortran complex
    literal whose components are non-numeric — gfortran rejects this
    in PARAMETER and array-constructor contexts. Wrap with
    ``complex128x2(...)`` so the structure constructor takes over.

    The same fix-up applies to literals whose components have been
    rewritten to ``float64x2('...')`` calls (e.g. when the original
    source used ``(1.0_WP, 0.0_WP)``): the parenthesized pair is no
    longer a recognized Fortran complex literal once the components
    contain function calls, so the outer wrap is necessary.
    """
    if target_mode.is_kind_based or not target_mode.complex_constructor:
        return line
    # Use named-component form (``re=..., im=...``) so the result
    # binds to the complex128x2 *structure constructor* and not to the
    # generic interface — the latter is a function call and is illegal
    # inside PARAMETER initializers (e.g. ``ZONE = complex128x2(re=...,
    # im=...)`` in DGEDMD/UGEDMD). The negative lookbehind on the
    # opening paren prevents re-wrapping an existing
    # ``complex128x2(...)`` call, e.g. one that the migrator already
    # produced from ``DCMPLX(ONE, ZERO)``.
    line = re.sub(
        r'(?<![A-Za-z0-9_])\(\s*([-+]?\s*(?:MF|DD)_[A-Z][A-Z0-9_]*)\s*,\s*'
        r'([-+]?\s*(?:MF|DD)_[A-Z][A-Z0-9_]*)\s*\)',
        rf"{target_mode.complex_constructor}(re=\1, im=\2)",
        line,
    )
    real_ctor = re.escape(target_mode.real_constructor or '')
    line = re.sub(
        rf"(?<![A-Za-z0-9_])\(\s*([-+]?\s*{real_ctor}\([^()]*\))\s*,"
        rf"\s*([-+]?\s*{real_ctor}\([^()]*\))\s*\)",
        rf"{target_mode.complex_constructor}(re=\1, im=\2)",
        line,
    )
    return line


_XERBLA_STR_RE = re.compile(r"'([A-Za-z][A-Za-z0-9_]*)( ?)'")


def replace_xerbla_strings(line: str, rename_map: dict[str, str]) -> str:
    """Replace routine names inside XERBLA string arguments.

    XERBLA's first argument is a Fortran string literal naming the
    failing routine — e.g. ``CALL XERBLA('DGEMM ', 1)``. A routine
    rename must rewrite that literal too. Previously this function
    looped over the whole rename_map and called ``str.replace`` twice
    per entry per line, which is O(N * lines) — catastrophic once the
    map hits MUMPS scale (6k+ entries). The regex version below is
    O(lines) with a dict lookup per quoted identifier found.
    """
    def sub(m):
        name = m.group(1)
        upper = name.upper()
        new = rename_map.get(upper)
        if new is None:
            return m.group(0)
        return f"'{new.upper()}{m.group(2)}'"
    return _XERBLA_STR_RE.sub(sub, line)


_FP_VALUE_RE = re.compile(
    r'^[+-]?\s*(?:'
    r'\d+\.\d*|\d*\.\d+|\d+'
    r')(?:[DdEe][+-]?\d+)?\s*$'
)


def _is_fp_value(text: str, known_constants: dict) -> bool:
    """Heuristic: does the trimmed PARAMETER/DATA value look like an
    FP literal, a known-constant reference, or an expression that
    contains one?

    Tighter than a substring scan: rejects identifiers like ``ELEMENT``
    that happen to contain ``E``. The composite-expression case splits
    on Fortran operators and checks each operand: if ANY operand is a
    known FP constant or an FP literal, the whole expression counts.
    """
    s = text.strip()
    if not s:
        return False
    if _FP_VALUE_RE.match(s):
        # Pure numeric token: FP iff it has a decimal point or exponent.
        return ('.' in s) or ('D' in s.upper()) or ('E' in s.upper())
    # Tokenize on identifiers + numeric literals; ignore operators.
    tokens = re.findall(r'[A-Za-z_]\w*|\d+\.\d*[DdEe][+-]?\d+|\d*\.\d+[DdEe]?[+-]?\d*|\d+\.\d*|\d*\.\d+', s)
    for tok in tokens:
        if tok.upper() in known_constants:
            return True
        if _FP_VALUE_RE.match(tok) and (('.' in tok) or ('D' in tok.upper()) or ('E' in tok.upper())):
            return True
    return False


def _join_continued_lines(lines: list[str], start: int) -> tuple[str, int]:
    """Join a fixed-form statement starting at ``lines[start]`` with any
    continuation lines (column-6 marker). Returns ``(joined_text, end)``
    where ``end`` is the index of the first line NOT consumed.
    """
    out = lines[start].rstrip('\n')
    j = start + 1
    while j < len(lines):
        nxt = lines[j]
        if is_continuation_line(nxt):
            out += ' ' + nxt[6:].rstrip('\n')
            j += 1
            continue
        break
    return out, j


_PROC_HEADER_RE_SCOPE = re.compile(
    r'^\s*(?:RECURSIVE\s+|PURE\s+|ELEMENTAL\s+)*'
    r'(?:(?:INTEGER|REAL|COMPLEX|LOGICAL|CHARACTER|TYPE\s*\([^)]+\)|DOUBLE\s+PRECISION|DOUBLE\s+COMPLEX)'
    r'(?:\s*\*\s*\d+)?\s+)?'
    r'(?:PROGRAM|SUBROUTINE|FUNCTION|MODULE|BLOCK\s+DATA)\b',
    re.IGNORECASE,
)


_INTERFACE_BEGIN_RE = re.compile(r'^\s*(?:ABSTRACT\s+)?INTERFACE\b', re.IGNORECASE)
_INTERFACE_END_RE = re.compile(r'^\s*END\s*INTERFACE\b', re.IGNORECASE)


def _scope_indices(lines: list[str]) -> list[int]:
    """Return a per-line list mapping line index → procedure scope index.

    Scope index 0 is the first SUBROUTINE/FUNCTION header encountered
    when scanning forward from the top of the file. Lines before any
    header are scope -1 (module/global level). Used by
    ``convert_parameter_stmts`` / ``convert_data_stmts`` to tag each
    converted assignment with the scope it belongs to, so that
    ``insert_use_multifloats`` can insert only the assignments
    belonging to the current scope.

    INTERFACE-block inner SUBROUTINE/FUNCTION declarations are NOT
    counted as new scopes — they declare prototypes for external
    procedures, not local scopes that take runtime assignments.

    Replaces the prior per-call ``_scope_index_at`` helper that
    rescanned ``lines[0..i]`` on every invocation (O(N²) when called
    from a per-statement outer loop). Compute the full vector once
    in O(N) and index into it.
    """
    scopes: list[int] = []
    scope = -1
    in_interface = 0
    for ln in lines:
        if _INTERFACE_BEGIN_RE.match(ln):
            in_interface += 1
        elif _INTERFACE_END_RE.match(ln):
            if in_interface > 0:
                in_interface -= 1
        elif in_interface == 0 and _PROC_HEADER_RE_SCOPE.match(ln):
            scope += 1
        scopes.append(scope)
    return scopes


def convert_parameter_stmts(
    source: str, target_mode: TargetMode,
) -> tuple[str, list[tuple[int, str]], dict[str, str]]:
    """Convert floating-point PARAMETER statements to executable assignments.

    Returns ``(new_source, fp_assignments, dropped_known)`` where
    ``fp_assignments`` is a list of ``(scope_index, assignment_text)``
    tuples so the caller can insert each assignment into the correct
    procedure scope. ``dropped_known`` maps each known-constant name
    skipped from the PARAMETER list to its multifloats replacement (so
    the caller can add it to the per-file rename set).

    Multi-line PARAMETER statements (fixed-form column-6 continuation)
    are joined into a single logical statement before parsing. The
    original line(s) are replaced as a unit so the line count of the
    output may differ from the input.
    """
    if target_mode.is_kind_based:
        return source, [], {}

    # Pre-scan declarations so we can tell whether a name like ``ONE``
    # was declared COMPLEX. The original LAPACK convention in
    # Z-prefixed routines is ``COMPLEX*16 ONE; PARAMETER(ONE = 1.0D+0)``
    # — the value is a real literal but it carries complex semantics
    # because Fortran promotes the literal to the declared type. The
    # multifloats migrator must NOT fold such names into the real
    # ``MF_ONE`` constant.
    complex_names = _scan_complex_var_names(source)

    lines = source.splitlines(keepends=True)
    scope_vec = _scope_indices(lines)
    result, fp_assignments = [], []
    dropped_known: dict[str, str] = {}
    param_re = re.compile(r'^(\s{6,}|^\s*)PARAMETER\s*\((.*)\)\s*(!.*)?$', re.IGNORECASE)
    # Combined declaration-attribute form: ``TYPE, PARAMETER :: name = val``.
    # The type-spec captured in group 2 is whatever comes before the
    # ``, PARAMETER ::`` token. After matching, the type-spec is
    # preserved for the new declaration line and each ``name = val``
    # entry is split off as a runtime assignment — same logic as the
    # standalone-PARAMETER branch.
    combined_param_re = re.compile(
        r'^(\s+)('
        r'(?:DOUBLE\s+PRECISION|DOUBLE\s+COMPLEX'
        r'|REAL\s*\*\s*\d+|COMPLEX\s*\*\s*\d+'
        r'|REAL\s*\(\s*[^)]*\)|COMPLEX\s*\(\s*[^)]*\)'
        r'|TYPE\s*\(\s*[A-Za-z_]\w*\s*\)'
        r'|INTEGER\s*\(\s*[^)]*\)|INTEGER'
        r'|REAL|COMPLEX)'
        r')\s*,\s*PARAMETER\s*::\s*(.+?)\s*(!.*)?$',
        re.IGNORECASE,
    )

    # CPP guard stack: each frame is the effective ``#if ...`` directive
    # line (without trailing newline) that should wrap any runtime
    # assignment emitted under the current source position, plus the
    # history of conditions seen in this if-block so ``#elif``/``#else``
    # branches can negate them.
    #
    # When a PARAMETER is converted to a runtime assignment, the
    # assignment is wrapped in the same #if/#endif so it stays in scope
    # only when the original declaration is in scope — otherwise gfortran
    # sees an assignment to an undeclared symbol. Tracking
    # ``#else``/``#elif`` is required so an assignment in the ``#else``
    # branch isn't wrapped under the (unrelated) ``#if`` true branch.
    cpp_stack: list[dict] = []

    def _cpp_cond(directive: str) -> str:
        s = directive.strip()
        if s.startswith('#ifdef'):
            return f'defined({s[len("#ifdef"):].strip()})'
        if s.startswith('#ifndef'):
            return f'!defined({s[len("#ifndef"):].strip()})'
        if s.startswith('#elif'):
            return s[len('#elif'):].strip()
        if s.startswith('#if'):
            return s[len('#if'):].strip()
        return ''

    i = 0
    while i < len(lines):
        line = lines[i]
        stripped_for_cpp = line.lstrip()
        if stripped_for_cpp.startswith('#if'):
            cpp_stack.append({
                'wrap': line.rstrip('\n'),
                'history': [_cpp_cond(line)],
            })
        elif stripped_for_cpp.startswith('#elif'):
            if cpp_stack:
                frame = cpp_stack[-1]
                negated = ' || '.join(f'({c})' for c in frame['history'])
                new_cond = _cpp_cond(line)
                frame['history'].append(new_cond)
                frame['wrap'] = f'#if !({negated}) && ({new_cond})'
        elif stripped_for_cpp.startswith('#else'):
            if cpp_stack:
                frame = cpp_stack[-1]
                negated = ' || '.join(f'({c})' for c in frame['history'])
                frame['wrap'] = f'#if !({negated})'
        elif stripped_for_cpp.startswith('#endif'):
            if cpp_stack: cpp_stack.pop()

        # Try matching a single-line PARAMETER first; if not, try
        # joining continuation lines and matching the joined form.
        joined, next_i = line.rstrip('\n'), i + 1
        if param_re.match(joined) is None and re.match(r'^\s{6,}PARAMETER\b', joined, re.IGNORECASE):
            joined, next_i = _join_continued_lines(lines, i)

        # Combined-form ``TYPE, PARAMETER :: name = val [, ...]``.
        # Split into (a) a plain declaration line listing only the
        # variable names and (b) a runtime assignment for each entry.
        # For known-constant names with non-complex values, drop the
        # entry entirely (the multifloats module supplies the constant).
        cm = combined_param_re.match(joined)
        if cm:
            indent  = cm.group(1)
            type_sp = cm.group(2)
            items   = cm.group(3)
            comment = cm.group(4) or ''

            entries: list[str] = []
            cur, depth = '', 0
            for ch in items:
                if ch == '(': depth += 1; cur += ch
                elif ch == ')': depth -= 1; cur += ch
                elif ch == ',' and depth == 0:
                    if cur.strip(): entries.append(cur.strip())
                    cur = ''
                else:
                    cur += ch
            if cur.strip():
                entries.append(cur.strip())

            kept_names: list[str] = []
            line_assignments: list[str] = []
            line_dropped_known: dict[str, str] = {}
            type_is_complex = bool(re.search(r'COMPLEX|cmplx64x2', type_sp, re.IGNORECASE))
            type_is_integer = bool(re.match(r'\s*INTEGER', type_sp, re.IGNORECASE))

            # Standard INTEGER PARAMETER lines (``INTEGER, PARAMETER ::
            # N = 5`` / ``wp = kind(1.d0)`` etc.) are constant
            # expressions and must stay as PARAMETERs — converting them
            # to runtime assignments would break things like the
            # ``kind(1.d0)`` idiom (the literal gets wrapped to a
            # derived-type value the ``kind()`` intrinsic rejects).
            #
            # The exception is the upstream-bug case
            # ``INTEGER, PARAMETER :: ZERO = (0.0D0, 0.0D0)``
            # (zsol_fwd_aux.F:1095): gfortran tolerates these via
            # numeric coercion, but multifloats mode rewrites the
            # literal to ``cmplx64x2(...)`` which is not a constant
            # expression for INTEGER. When we detect a complex-literal
            # value on an INTEGER LHS, override the declared type to
            # the multifloats complex type so the runtime assignment
            # is type-correct (the use sites are downstream complex
            # contexts in every observed case).
            if type_is_integer:
                has_cx_lit = any(
                    ('(' in v.split('=', 1)[1] and ',' in v.split('=', 1)[1])
                    for v in entries if '=' in v
                )
                if has_cx_lit:
                    type_sp = target_mode.complex_type
                    type_is_complex = True
                    type_is_integer = False
                else:
                    # Bail out — the standalone-PARAMETER branch below
                    # won't match either, so the line falls through to
                    # the unchanged-emit path.
                    result.append(line)
                    i += 1
                    continue
            for entry in entries:
                if '=' in entry:
                    name, val = entry.split('=', 1)
                    name, val = name.strip(), val.strip()
                else:
                    # No initializer — a PARAMETER without a value is
                    # not legal Fortran, so this branch should not
                    # trigger in practice. Bail out by emitting the
                    # original line unchanged.
                    kept_names = []
                    break

                if not _is_fp_value(val, target_mode.known_constants):
                    # INTEGER PARAMETER, character literal, etc. don't
                    # need conversion — leave the combined-form line
                    # alone (we'll fall through to the unchanged emit
                    # path below).
                    kept_names = []
                    break

                cx_ctor = (target_mode.complex_constructor or '').lower()
                is_cx_value = (
                    type_is_complex
                    or ('(' in val and ',' in val)
                    or (cx_ctor and cx_ctor in val.lower())
                    or 'cmplx' in val.lower()
                    or 'dcmplx' in val.lower()
                    or name.upper() in complex_names
                )
                if (name.upper() in target_mode.known_constants
                        and not is_cx_value):
                    line_dropped_known[name.upper()] = (
                        target_mode.known_constants[name.upper()]
                    )
                    continue

                kept_names.append(name)
                assn = f"{indent}{name} = {val}{comment}\n"
                if cpp_stack:
                    pre = ''.join(f["wrap"] + '\n' for f in cpp_stack)
                    post = '#endif\n' * len(cpp_stack)
                    assn = pre + assn + post
                line_assignments.append(assn)
            else:
                # Loop completed without `break`: all entries were
                # FP-valued and either dropped or converted. Emit the
                # decl line + assignments and skip past the consumed
                # source lines.
                scope = scope_vec[i]
                fp_assignments.extend((scope, a) for a in line_assignments)
                dropped_known.update(line_dropped_known)
                if kept_names:
                    decl = f"{indent}{type_sp} :: {', '.join(kept_names)}{comment}\n"
                    result.append(decl)
                else:
                    # Every entry was a known constant — drop the line
                    # entirely.
                    pass
                i = next_i
                continue
            # Bailed out — leave the source as-is for this line.

        m = param_re.match(joined)
        if m:
            indent, params_content, comment = m.group(1), m.group(2), m.group(3) or ''
            parts, current, depth = [], [], 0
            for char in params_content:
                if char == '(': depth += 1
                elif char == ')': depth -= 1
                if char == ',' and depth == 0:
                    parts.append(''.join(current))
                    current = []
                else: current.append(char)
            if current: parts.append(''.join(current))

            kept_parts = []
            line_assignments: list[str] = []
            line_dropped_known: dict[str, str] = {}
            for part in parts:
                if '=' in part:
                    name, val = part.split('=', 1)
                    name, val = name.strip(), val.strip()
                    if _is_fp_value(val, target_mode.known_constants):
                        # A known-constant name carries complex
                        # semantics if either (a) the value is a
                        # complex literal / constructor, or (b) the
                        # local declaration of the name is COMPLEX.
                        # In either case it must NOT be folded into
                        # the multifloats real-constant rename map —
                        # we keep it as a runtime assignment so the
                        # variable retains its complex type.
                        cx_ctor = (target_mode.complex_constructor or '').lower()
                        is_cx_value = (
                            ('(' in val and ',' in val) or
                            (cx_ctor and cx_ctor in val.lower()) or
                            'cmplx' in val.lower() or
                            'dcmplx' in val.lower() or
                            name.upper() in complex_names
                        )
                        if (name.upper() in target_mode.known_constants
                                and not is_cx_value):
                            line_dropped_known[name.upper()] = target_mode.known_constants[name.upper()]
                            continue
                        assn = f"{indent}{name} = {val}{comment}\n"
                        if cpp_stack:
                            pre = ''.join(f["wrap"] + '\n' for f in cpp_stack)
                            post = '#endif\n' * len(cpp_stack)
                            assn = pre + assn + post
                        line_assignments.append(assn)
                    else: kept_parts.append(part)
                else: kept_parts.append(part)

            scope = scope_vec[i]
            fp_assignments.extend((scope, a) for a in line_assignments)
            dropped_known.update(line_dropped_known)
            if kept_parts:
                result.append(f"{indent}PARAMETER ({', '.join(kept_parts)}){comment}\n")
            elif line_assignments:
                # Some FP entries became runtime assignments — leave a
                # short marker comment so reviewers can find the source.
                result.append(f"{indent}! Converted to assignments below: {joined.strip()}\n")
            # else: every entry was a known constant supplied by the
            # multifloats module — drop the line entirely (no comment).
            i = next_i
            continue
        result.append(line)
        i += 1
    return "".join(result), fp_assignments, dropped_known


def convert_data_stmts(
    source: str, target_mode: TargetMode,
) -> tuple[str, list[tuple[int, str]], dict[str, str]]:
    """Convert floating-point DATA statements to executable assignments.

    Returns ``(new_source, fp_assignments, dropped_known)`` — see
    :func:`convert_parameter_stmts` for the meaning of the third tuple
    element.
    """
    if target_mode.is_kind_based:
        return source, [], {}

    lines = source.splitlines(keepends=True)
    scope_vec = _scope_indices(lines)
    result, fp_assignments = [], []
    dropped_known: dict[str, str] = {}
    data_re = re.compile(r'^(\s{6,}|^\s*)DATA\s+([^/]+)/\s*([^/]+)\s*/\s*(!.*)?$', re.IGNORECASE)

    i = 0
    while i < len(lines):
        line = lines[i]
        joined, next_i = line.rstrip('\n'), i + 1
        if data_re.match(joined) is None and re.match(r'^\s{6,}DATA\b', joined, re.IGNORECASE):
            joined, next_i = _join_continued_lines(lines, i)

        m = data_re.match(joined)
        if m:
            indent, vars_part, vals_part, comment = m.group(1), m.group(2).strip(), m.group(3).strip(), m.group(4) or ''
            # Each value is FP iff it independently matches the FP
            # heuristic; this is more discriminating than scanning the
            # whole vals_part for ``D``/``E`` substrings (which would
            # falsely match identifiers).
            tmp_vals: list[str] = []
            cur, depth = '', 0
            for ch in vals_part:
                if ch == '(': depth += 1
                elif ch == ')': depth -= 1
                if ch == ',' and depth == 0:
                    tmp_vals.append(cur.strip()); cur = ''
                else:
                    cur += ch
            if cur.strip(): tmp_vals.append(cur.strip())
            any_fp = any(_is_fp_value(v, target_mode.known_constants) for v in tmp_vals)

            if any_fp:
                vars_list = [v.strip() for v in vars_part.split(',')]
                if len(vars_list) == len(tmp_vals):
                    line_assignments: list[str] = []
                    for v, val in zip(vars_list, tmp_vals):
                        if v.upper() in target_mode.known_constants:
                            dropped_known[v.upper()] = target_mode.known_constants[v.upper()]
                            continue
                        line_assignments.append(f"{indent}{v} = {val}{comment}\n")
                    scope = scope_vec[i]
                    fp_assignments.extend((scope, a) for a in line_assignments)
                    if line_assignments:
                        result.append(f"{indent}! Converted to assignments below: {joined.strip()}\n")
                    # else: every name was a known constant — drop the line
                    i = next_i
                    continue
            result.append(line)
            i += 1
            continue
        result.append(line)
        i += 1
    return "".join(result), fp_assignments, dropped_known


_STMT_FN_RE = re.compile(r'^[A-Za-z_]\w*\s*\(\s*[A-Za-z_]\w*\s*\)\s*=')


def _looks_like_statement_function(stripped: str, lines: list[str], k: int) -> bool:
    """Heuristic: does line k look like a LAPACK statement function
    definition that the declaration-block walker should step over?

    A statement function has the form ``NAME(SCALAR_ARG) = expression``
    and appears in the LAPACK source between a comment marker like
    ``*     .. Statement Function definitions ..`` and the executable
    statements section. We detect by looking back at recent lines for
    that marker (within ~10 lines).
    """
    if not _STMT_FN_RE.match(stripped):
        return False
    look = max(0, k - 12)
    for kk in range(k - 1, look - 1, -1):
        prev = lines[kk]
        if not prev.strip():
            continue
        if 'Statement Function' in prev:
            return True
        if prev and prev[0] in ('C', 'c', '*', '!'):
            continue
        if prev.lstrip().startswith('!'):
            continue
        # Hit a code line that wasn't a statement function — stop.
        if not _STMT_FN_RE.match(prev.lstrip()):
            return False
    return False


# Module public names (type names, constants, generics, operator generics)
# are now loaded from the target YAML via TargetMode fields:
#   target_mode.module_type_names
#   target_mode.module_constant_names
#   target_mode.module_generic_names
#   target_mode.module_public_names  (union of the above three)
#   target_mode.module_operator_generics

_DECL_LINE_RE = re.compile(
    r'^\s+(?:TYPE\s*\([^)]*\)|INTEGER\b|REAL\b|COMPLEX\b|LOGICAL\b|'
    r'CHARACTER\b|DOUBLE\s+PRECISION\b|DOUBLE\s+COMPLEX\b)',
    re.IGNORECASE,
)
_IDENT_RE = re.compile(r"\b([A-Za-z_]\w*)\b")
_STRING_RE = re.compile(r"'(?:[^']|'')*'|\"(?:[^\"]|\"\")*\"")


def _strip_strings_and_comments(line: str) -> str:
    """Drop string literals and trailing inline comments."""
    out = _STRING_RE.sub('', line)
    bang = out.find('!')
    if bang >= 0:
        out = out[:bang]
    return out


def _scan_local_declared_names(proc_lines: list[str]) -> set[str]:
    """Collect local variable names from type-declaration statements
    inside a procedure. Used to suppress matching multifloats public
    names so that the local variable can shadow the use-associated
    generic interface (gfortran refuses if the name is in scope).
    """
    names: set[str] = set()
    joined: list[str] = []
    cur = ''
    for raw in proc_lines:
        # Skip pure comment lines (fixed and free form).
        if raw[:1] in ('C', 'c', '*', '!'):
            continue
        if cur and is_continuation_line(raw):
            # Fixed-form continuation.
            cur += ' ' + raw[6:].rstrip('\n')
            continue
        if cur and cur.rstrip().endswith('&'):
            cur = cur.rstrip().rstrip('&') + ' ' + raw.lstrip().lstrip('&').rstrip('\n')
            continue
        if cur:
            joined.append(cur)
        cur = raw.rstrip('\n')
    if cur:
        joined.append(cur)

    for stmt in joined:
        if not _DECL_LINE_RE.match(stmt):
            continue
        # Strip the type prefix (everything up to and including the
        # first ``::`` if present, otherwise up to the type keyword).
        if '::' in stmt:
            tail = stmt.split('::', 1)[1]
        else:
            m = _DECL_LINE_RE.match(stmt)
            tail = stmt[m.end():]
        tail = _strip_strings_and_comments(tail)
        # Drop array specs and KIND parameters in parentheses so we
        # don't pick up bound expressions (``WNRM(MAX(M,N))``) as
        # local-variable names. Iterate until no more nested parens
        # remain.
        while True:
            new_tail = re.sub(r'\([^()]*\)', '', tail)
            if new_tail == tail:
                break
            tail = new_tail
        # Drop ``= initializer`` clauses (PARAMETER initializers may
        # contain references to module names like ``complex128x2``
        # that should not be treated as locally-declared variables).
        # Split at top-level commas, then drop everything from ``=``
        # onward in each item.
        items = []
        cur = ''
        for ch in tail + ',':
            if ch == ',':
                items.append(cur)
                cur = ''
            else:
                cur += ch
        for item in items:
            lhs = item.split('=', 1)[0]
            for m in _IDENT_RE.finditer(lhs):
                names.add(m.group(1).lower())
    return names


def _scan_referenced_identifiers(proc_lines: list[str]) -> set[str]:
    """Lower-cased identifiers referenced anywhere in the procedure
    body, excluding comments and string literals."""
    names: set[str] = set()
    for raw in proc_lines:
        if raw[:1] in ('C', 'c', '*', '!'):
            continue
        cleaned = _strip_strings_and_comments(raw)
        for m in _IDENT_RE.finditer(cleaned):
            names.add(m.group(1).lower())
    return names


# Procedure header for SUBROUTINE/FUNCTION/PROGRAM/MODULE/BLOCK DATA.
# ``MODULE PROCEDURE`` (inside an INTERFACE block) is NOT a procedure
# header and is excluded — injecting USE between ``MODULE PROCEDURE foo``
# and ``END INTERFACE`` is illegal Fortran.
_PROC_HEADER_RE = re.compile(
    r'^(\s{6,}|^\s*)(?:RECURSIVE\s+|PURE\s+|ELEMENTAL\s+)*'
    r'(?:(?:INTEGER|REAL|COMPLEX|LOGICAL|CHARACTER|TYPE\s*\([^)]+\)|DOUBLE\s+PRECISION|DOUBLE\s+COMPLEX)'
    r'(?:\s*\*\s*\d+)?\s+)?'
    r'(?:PROGRAM|SUBROUTINE|FUNCTION|MODULE(?!\s+PROCEDURE\b)|BLOCK\s+DATA)\b',
    re.IGNORECASE,
)
# ``END SUBROUTINE FOO``, ``END FUNCTION``, plain ``END``. Crucially
# NOT ``END IF`` / ``END DO`` / ``END SELECT`` — the keyword whitelist
# must be required when any word follows ``END``, otherwise inner
# control-flow ENDs would falsely terminate body scans.
_END_PROC_RE = re.compile(
    r'^\s*END(?:\s+(?:PROGRAM|SUBROUTINE|FUNCTION|MODULE|BLOCK\s*DATA)'
    r'(?:\s+\w+)?)?\s*(?:!.*)?$',
    re.IGNORECASE,
)
def specialize_use_module(source: str, target_mode: TargetMode, fixed_form: bool) -> str:
    """Replace bare ``USE <module>`` clauses with explicit ``only:``
    lists tailored to each procedure.

    Operates on the fully-migrated source so that scanned identifiers
    reflect the post-transform names. Each procedure body between its
    header and matching END statement is scanned for referenced module
    names; the only-list excludes any name that the procedure declares
    as a local variable, so that LAPACK locals like ``SUM``/``SCALE``
    etc. do not collide with use-associated generic interfaces.
    """
    if not target_mode.module_name or target_mode.intrinsic_mode != 'wrap_constructor':
        return source

    use_mod_re = re.compile(
        rf'^(\s*)USE\s+{re.escape(target_mode.module_name)}\s*(?:!.*)?$',
        re.IGNORECASE,
    )

    lines = source.splitlines(keepends=True)
    # Find all procedure boundaries. INTERFACE-inner SUBROUTINE/
    # FUNCTION declarations are NOT real procedures — they declare
    # external prototypes. Skip them so the outer procedure's
    # body-text scan covers the real implementation, not just the
    # span up to the first interface-inner END.
    headers: list[int] = []
    ends: list[int] = []
    in_interface = 0
    for idx, ln in enumerate(lines):
        if _INTERFACE_BEGIN_RE.match(ln):
            in_interface += 1
            continue
        if _INTERFACE_END_RE.match(ln):
            if in_interface > 0:
                in_interface -= 1
            continue
        if in_interface > 0:
            continue
        if _PROC_HEADER_RE.match(ln):
            headers.append(idx)
        if _END_PROC_RE.match(ln):
            ends.append(idx)

    # Pair headers with their matching END (the next END at or after
    # the header). This is approximate but works for BLAS/LAPACK where
    # CONTAINS is rare.
    out = list(lines)
    for h_idx, h in enumerate(headers):
        # Determine the end of this procedure.
        next_end = next((e for e in ends if e >= h), len(lines) - 1)
        proc_lines = lines[h:next_end + 1]
        only_clause = _build_use_only_clause(proc_lines, target_mode)
        if not only_clause:
            continue
        # Replace the bare ``USE <module>`` line(s) inside this
        # procedure with the only-form. There should be exactly one.
        for k in range(h, next_end + 1):
            m = use_mod_re.match(out[k])
            if not m:
                continue
            indent = m.group(1) or ('      ' if fixed_form else '    ')
            body = f"{target_mode.module_name}{only_clause}"
            wrapped = _wrap_use_clause(indent, body, fixed_form)
            out[k] = wrapped
    return ''.join(out)


def _wrap_use_clause(indent: str, body: str, fixed_form: bool) -> str:
    """Wrap a long ``USE multifloats, only: a, b, c, ...`` statement so
    each emitted line fits within Fortran source line limits.

    ``body`` is the post-USE text (``multifloats, only: a, b, ...``).
    Items are kept comma-separated; when adding the next item would
    overflow the line, the current line is flushed and a new
    continuation line is started.
    """
    cap = 72 if fixed_form else 132
    cont_prefix = '     +' if fixed_form else (indent + '   ')

    full = indent + 'USE ' + body
    if len(full) <= cap:
        return full + '\n'

    # Split body into ``head`` (``multifloats``) and the only-list.
    head, sep, rest = body.partition(', only:')
    if not sep:
        return full + '\n'

    parts: list[str] = []
    cur = ''
    depth = 0
    for ch in rest.strip():
        if ch == ',' and depth == 0:
            if cur.strip():
                parts.append(cur.strip())
            cur = ''
            continue
        if ch == '(':
            depth += 1
        elif ch == ')':
            depth -= 1
        cur += ch
    if cur.strip():
        parts.append(cur.strip())

    out_lines: list[str] = []
    cur_line = indent + 'USE ' + head + ', only:'
    # Reserve room on each line for the continuation overhead: a
    # trailing ``,`` (fixed form) or ``, &`` (free form).
    cont_overhead = 1 if fixed_form else 3
    for part in parts:
        addition = (' ' if cur_line.rstrip().endswith(':') else ', ') + part
        if len(cur_line) + len(addition) + cont_overhead <= cap:
            cur_line += addition
            continue
        # Flush — both forms need a trailing comma so that the next
        # only-list item is correctly comma-separated. Free form also
        # appends an explicit ``&`` continuation marker.
        if not cur_line.rstrip().endswith(','):
            cur_line = cur_line.rstrip() + ','
        if fixed_form:
            out_lines.append(cur_line + '\n')
        else:
            out_lines.append(cur_line + ' &\n')
        cur_line = cont_prefix + part
    out_lines.append(cur_line + '\n')
    return ''.join(out_lines)


def _build_use_only_clause(proc_lines: list[str], target_mode: TargetMode) -> str:
    """Compute the ``, only:`` clause for a module-based USE statement.

    Returns the empty string if the target does not use a module. Otherwise
    returns ``", only: name1, name2, ..., operator(+), ..."`` listing
    the module public names referenced by ``proc_lines`` (minus any name
    that the procedure declares as a local variable, so that the local
    declaration is not shadowed by the use-associated generic interface).
    """
    if target_mode.intrinsic_mode != 'wrap_constructor':
        return ''
    referenced = _scan_referenced_identifiers(proc_lines)
    # If the procedure body pulls in content via ``INCLUDE 'xxx.h'``,
    # the walker can't see inside the header — but migrated struct
    # headers (e.g. MUMPS's ddmumps_struc.h) reference the target-mode
    # type names. Unconditionally surface those type names so the
    # host module's USE clause imports them. Type names rarely collide
    # with local variables so the over-inclusion is safe.
    for raw in proc_lines:
        if raw[:1] in ('C', 'c', '*', '!'):
            continue
        if _INCLUDE_RE.match(raw):
            referenced |= {n.lower() for n in target_mode.module_type_names}
            break
    # Restore identifiers masked by the keep-kind sentinel before the
    # only-list is built. specialize_use_module runs while sentinels are
    # still in place (the restore happens after migrate_*_form returns),
    # so a body like ``CALL F(__KEEPKIND_DBLE__(PEAK))`` would otherwise
    # not contribute ``dble`` to the only-list — and the post-restore
    # ``CALL F(dble(PEAK))`` then dispatches to gfortran's intrinsic
    # ``dble`` which doesn't accept the multifloats real64x2 type.
    body_text = ''.join(proc_lines)
    if _KK_DBLE_SENTINEL in body_text:
        referenced.add('dble')
    # Predict the ``dble`` calls that ``_rewrite_int_kind_on_real64x2``
    # will inject in a later post-pass: ``INT(real64x2_expr, K)`` →
    # ``INT(dble(real64x2_expr), K)``. The rewrite runs AFTER the
    # only-clause is built, so the scanner can't see the literal ``dble``
    # call yet. Conservatively add ``dble`` whenever both ``int(`` and
    # ``real64x2(`` appear in the body — false positives just import an
    # unused name, which gfortran tolerates.
    if (target_mode.intrinsic_mode == 'wrap_constructor'
            and re.search(r'\bint\s*\(', body_text, re.IGNORECASE)
            and 'real64x2' in body_text.lower()):
        # Match either the literal-wrapped form ``real64x2(...)`` or the
        # bare type ``TYPE(real64x2)`` declaration — the latter signals
        # that some local is real64x2-typed and may be passed through
        # ``int(var, K)``, which the rewrite then routes via ``dble``.
        referenced.add('dble')
    declared = _scan_local_declared_names(proc_lines)
    # Determine constant name prefix for sorting (e.g. 'mf_' for multifloats)
    const_prefixes = set()
    for cn in target_mode.module_constant_names:
        idx = cn.find('_')
        if idx >= 0:
            const_prefixes.add(cn[:idx + 1])
    def _sort_key(s: str) -> tuple:
        return (any(s.startswith(p) for p in const_prefixes), s)
    selected = sorted(
        (referenced & target_mode.module_public_names) - declared,
        key=_sort_key,
    )
    parts = list(selected) + list(target_mode.module_operator_generics)
    return ', only: ' + ', '.join(parts) if parts else ''


def insert_use_multifloats(source: str, target_mode: TargetMode,
                           extra_lines: list[tuple[int, str]] | list[str] | None = None) -> str:
    """Insert USE multifloats statement and extra assignments after procedure headers.

    The USE statement is emitted with an explicit ``only:`` clause that
    lists exactly the multifloats public names referenced by the
    enclosing procedure (plus the operator/assignment generics, which
    are always included). This avoids importing names like ``sum``,
    ``scale``, ``gamma``, ``tiny``, ``nint``, etc. that LAPACK uses as
    local variable names — gfortran rejects local declarations that
    collide with use-associated generic interfaces.
    """
    if not target_mode.module_name and not extra_lines:
        return source

    lines = source.splitlines(keepends=True)
    result = []
    proc_header_re = _PROC_HEADER_RE
    end_proc_re = _END_PROC_RE

    # Normalise extra_lines: support both scoped (int, str) tuples and
    # legacy flat strings (all go to scope -1 which matches every scope
    # for backward compat with callers that don't use scoping).
    scoped_extras: list[tuple[int, str]] = []
    if extra_lines:
        for item in extra_lines:
            if isinstance(item, tuple):
                scoped_extras.append(item)
            else:
                scoped_extras.append((-1, item))

    scope_counter = -1
    in_interface = 0

    i = 0
    while i < len(lines):
        line = lines[i]
        result.append(line)
        # INTERFACE blocks bracket prototype SUBROUTINE/FUNCTION
        # declarations. The prototypes still need USE statements (so
        # they can see ``real64x2`` / ``cmplx64x2``) because each
        # interface body is its own scope and does NOT inherit the
        # enclosing procedure's USE clauses. But they must NOT receive
        # runtime assignments, and they don't count as procedure scopes
        # for assignment placement.
        if _INTERFACE_BEGIN_RE.match(line):
            in_interface += 1
            i += 1
            continue
        if _INTERFACE_END_RE.match(line):
            if in_interface > 0:
                in_interface -= 1
            i += 1
            continue
        m = proc_header_re.match(line)
        if m:
            inside_interface = in_interface > 0
            if not inside_interface:
                scope_counter += 1
            j = i + 1
            # Walk past continuation lines so the USE clause is
            # inserted AFTER the entire procedure header, not in the
            # middle of a continued SUBROUTINE/FUNCTION declaration.
            # Both fixed-form (column-6 marker) and free-form
            # (trailing ``&``) continuations are recognized.
            #
            # The SUBROUTINE header may also span CPP ``#if``/``#endif``
            # blocks that conditionalize formal arguments (common in
            # MUMPS — e.g. ``mana_aux.F`` toggles ``METIS_OPTIONS`` via
            # ``#if defined(metis)``). Track parenthesis depth across
            # the header lines: as long as the formal-arg-list paren is
            # still open, we treat the next non-CPP line as part of the
            # header even when neither the previous nor current line
            # carries a continuation marker.
            paren_depth = _count_open_parens(line)
            prev_has_amp = result[-1].rstrip().rstrip('\n').endswith('&')
            while j < len(lines):
                next_line = lines[j]
                if next_line.lstrip().startswith('#'):
                    result.append(next_line)
                    j += 1
                    continue
                if (is_continuation_line(next_line) or prev_has_amp
                        or paren_depth > 0):
                    result.append(next_line)
                    prev_has_amp = next_line.rstrip().rstrip('\n').endswith('&')
                    paren_depth += _count_open_parens(next_line)
                    j += 1
                else:
                    break

            if target_mode.module_name:
                indent = m.group(1)
                use_line = (
                    f"{indent if indent.strip() else '      '}"
                    f"USE {target_mode.module_name}\n"
                )
                already_has = any(f"USE {target_mode.module_name}".upper() in lines[kk].upper() for kk in range(j, min(j+20, len(lines))))
                if not already_has: result.append(use_line)

            # Filter extra_lines to those belonging to this scope
            # (scope_counter). Entries tagged with -1 are unscoped
            # (legacy) and go into every scope. Interface-inner
            # prototypes never get assignments — they declare
            # external interfaces, not local scopes.
            if inside_interface:
                scope_lines = []
            else:
                scope_lines = [text for sc, text in scoped_extras
                               if sc == scope_counter or sc == -1]

            if scope_lines:
                # Walk past the declaration block (blank lines, comments,
                # and any line whose first token is a declaration keyword
                # or a fixed-form continuation of one). Insert the
                # extra assignments at the END of the declaration block,
                # i.e. just before the first executable statement —
                # otherwise the assignments would be parsed as
                # implicit-typed executables before declarations.
                k = j
                prev_amp_decl = False
                while k < len(lines):
                    raw = lines[k]
                    stripped = raw.lstrip()
                    if not stripped.strip():
                        k += 1; continue
                    if raw and raw[0] in ('C', 'c', '*', '!'):
                        k += 1; continue
                    # An all-comment line in fixed-form may also start
                    # with whitespace then ``!`` (the inline-comment
                    # marker is legal at any column).
                    if stripped.startswith('!'):
                        k += 1; continue
                    # Preprocessor directives (``#if defined(_OPENMP)`` /
                    # ``#endif`` / ``#include``) appear inside the decl
                    # block in capital-F sources (e.g. dsytrd_sb2st.F has
                    # ``#if defined(_OPENMP) / use omp_lib / #endif``
                    # between the USE clause and IMPLICIT NONE). Without
                    # this, the walker would stop at the ``#if`` line
                    # — treating it as an executable — and insert the
                    # PARAMETER assignments above IMPLICIT NONE, which
                    # gfortran rejects as ``Unexpected IMPLICIT NONE
                    # statement``. The body of ``#if/#endif`` blocks in
                    # this position is itself decl-only (``use omp_lib``)
                    # so skipping the directive markers is enough.
                    if stripped.startswith('#'):
                        k += 1; continue
                    # Free-form continuation: previous code line ended
                    # in ``&``. The current line is part of the same
                    # logical statement (typically the trailing names
                    # on a multi-line ``INTEGER, INTENT(IN) :: ...``
                    # argument declaration). Skip it so the walker
                    # doesn't treat the continuation prefix as the
                    # start of an executable statement.
                    if prev_amp_decl:
                        prev_amp_decl = raw.rstrip().rstrip('\n').endswith('&')
                        k += 1; continue
                    if is_continuation_line(raw):
                        k += 1; continue
                    # INTERFACE / END INTERFACE blocks live in the
                    # declaration section. Walk past the entire block
                    # so the runtime assignments don't get inserted
                    # into a prototype body or between IMPLICIT NONE
                    # and the rest of the declarations. While walking,
                    # inject ``USE <module>`` after each inner
                    # SUBROUTINE/FUNCTION prototype header — interface
                    # body scopes do not inherit the enclosing
                    # procedure's USE clauses, so the migrated type
                    # references inside the prototype need their own
                    # module visibility.
                    if _INTERFACE_BEGIN_RE.match(raw):
                        k += 1
                        depth = 1
                        while k < len(lines) and depth > 0:
                            inner_ln = lines[k]
                            if _INTERFACE_BEGIN_RE.match(inner_ln):
                                depth += 1
                            elif _INTERFACE_END_RE.match(inner_ln):
                                depth -= 1
                            elif (target_mode.module_name
                                    and proc_header_re.match(inner_ln)):
                                # Walk continuation lines of the inner
                                # header so the USE lands AFTER the
                                # full prototype header, not inside
                                # the formal-arg list.
                                k_after = k + 1
                                paren_d = _count_open_parens(inner_ln)
                                prev_amp = inner_ln.rstrip().rstrip('\n').endswith('&')
                                while k_after < len(lines):
                                    nxt = lines[k_after]
                                    if nxt.lstrip().startswith('#'):
                                        k_after += 1
                                        continue
                                    if (is_continuation_line(nxt)
                                            or prev_amp or paren_d > 0):
                                        prev_amp = nxt.rstrip().rstrip('\n').endswith('&')
                                        paren_d += _count_open_parens(nxt)
                                        k_after += 1
                                    else:
                                        break
                                # Insert the bare ``USE`` line
                                # immediately after the prototype
                                # header. specialize_use_module's
                                # outer-procedure pass then rewrites
                                # the bare form to a ``, only:`` clause.
                                pm = proc_header_re.match(inner_ln)
                                hdr_indent = pm.group(1) if pm else '      '
                                use_line_inner = (
                                    f"{hdr_indent if hdr_indent.strip() else '      '}"
                                    f"USE {target_mode.module_name}\n"
                                )
                                # Splice the USE line into ``lines`` so
                                # the outer walker copies it through
                                # to ``result`` along with the rest of
                                # the prototype. Search 20 lines ahead
                                # for an existing USE to skip duplicates.
                                already_inner = any(
                                    f"USE {target_mode.module_name}".upper()
                                    in lines[kk].upper()
                                    for kk in range(k_after, min(k_after + 20, len(lines)))
                                )
                                if not already_inner:
                                    lines.insert(k_after, use_line_inner)
                                k = k_after
                                continue
                            k += 1
                        continue
                    l = stripped.upper()
                    if any(l.startswith(p) for p in (
                        'REAL', 'DOUBLE', 'COMPLEX', 'INTEGER', 'LOGICAL',
                        'CHARACTER', 'TYPE', 'USE', 'IMPLICIT', 'PARAMETER',
                        'DATA', 'INTRINSIC', 'EXTERNAL', 'DIMENSION', 'SAVE',
                        'EQUIVALENCE', 'COMMON', 'INCLUDE',
                        # Keep-kind sentinels appear at this stage in the
                        # pipeline — the migrator masks ``DOUBLE PRECISION`` /
                        # ``dble`` / ``dcmplx`` with sentinels around lines
                        # in recipes/*/keep-kind.manifest so they survive
                        # the migration unchanged. The sentinels get
                        # restored after migration but BEFORE this walker
                        # only when the migrator doesn't run a USE pass
                        # (e.g. the source-emit path), and the walker runs
                        # on the still-sentineled text in the multifloats
                        # path. Treat them as declaration keywords so the
                        # walker doesn't stop at them and insert
                        # PARAMETER-derived assignments in the middle of
                        # the declaration block.
                        '__KEEPKIND_DP__',
                        '__KEEPKIND_DBLE__',
                        '__KEEPKIND_DCMPLX__',
                    )):
                        prev_amp_decl = raw.rstrip().rstrip('\n').endswith('&')
                        k += 1; continue
                    # LAPACK statement-function definitions look like
                    # ``CABS1( ZDUM ) = ABS( ... )`` and live between the
                    # type-decl block and the executable statements. They
                    # are NOT executable, so the walker should also walk
                    # past them. We detect by looking back: if the
                    # previous non-blank line was a comment marked
                    # ``Statement Function`` (LAPACK convention) or this
                    # line is the only thing between two ``*     ..``
                    # separator comments, treat as still in decl section.
                    if _looks_like_statement_function(stripped, lines, k):
                        prev_amp_decl = raw.rstrip().rstrip('\n').endswith('&')
                        k += 1; continue
                    break
                # Copy declaration block as-is, then emit assignments.
                for kk in range(j, k):
                    result.append(lines[kk])
                for al in scope_lines:
                    result.append(al)
                result.append("\n")
                i = k
            else:
                i = j
            continue
        i += 1
    return "".join(result)


def is_comment_line(line: str) -> bool:
    return bool(line) and line[0] in ('C', 'c', '*', '!')

def _iter_outside_strings(text: str):
    """Yield ``(i, ch)`` for each character in ``text`` that is NOT
    inside a Fortran string literal. Recognizes both ``'`` and ``"``
    delimiters and the Fortran doubled-quote escape (``''`` inside a
    ``'``-string, ``""`` inside a ``"``-string).

    Canonical primitive shared by every scanner that needs to find
    inline ``!`` comments, count parens, or build a string-mask while
    respecting quoted-literal boundaries.
    """
    n = len(text)
    in_string = False
    quote = ''
    i = 0
    while i < n:
        ch = text[i]
        if in_string:
            if ch == quote:
                if i + 1 < n and text[i + 1] == quote:
                    # Doubled-quote escape — both chars stay in-string.
                    i += 2
                    continue
                in_string = False
        else:
            if ch in ("'", '"'):
                in_string = True
                quote = ch
            else:
                yield i, ch
        i += 1


def _find_inline_bang(text: str) -> int:
    """Return the index of the first inline ``!`` comment marker in
    ``text``, or ``len(text)`` if none. Quote-aware (doubled-quote
    escape included).

    Fast path: when the line has no string-delimiter characters at
    all (the common case for code-bearing source lines), fall through
    to ``str.find`` which runs in C. Only when a quote appears do we
    need the full quote-aware scan.
    """
    if "'" not in text and '"' not in text:
        idx = text.find('!')
        return idx if idx != -1 else len(text)
    for i, ch in _iter_outside_strings(text):
        if ch == '!':
            return i
    return len(text)


def _count_open_parens(line: str) -> int:
    """Net paren delta for ``line`` ignoring quoted strings and inline
    fixed-form ``!`` / ``C`` / ``*`` comments. Used to track when a
    SUBROUTINE/FUNCTION formal-arg list is still open across CPP
    ``#if/#endif`` blocks."""
    if not line:
        return 0
    if line[0] in ('C', 'c', '*'):
        return 0
    depth = 0
    for _, ch in _iter_outside_strings(line):
        if ch == '!':
            break
        if ch == '(':
            depth += 1
        elif ch == ')':
            depth -= 1
    return depth


def is_continuation_line(line: str) -> bool:
    # Fixed-form continuation: cols 1-5 blank (spaces only — NOT tabs), col 6
    # holds a non-blank, non-zero marker. A tab in col 1 is the gfortran/intel
    # extension that means "rest of line is normal source from col 7"; such a
    # line is a fresh statement, not a continuation, even when col 6 (i.e. the
    # 6th character after the tab) happens to be alphabetic. See dlarre2.f in
    # ScaLAPACK 2.2.3 for the wild form: `\t    CALL DLARRC(...)`.
    if not line or line[0] == '\t':
        return False
    return len(line) > 5 and line[0:5] == '     ' and line[5] not in (' ', '0', '')

def _build_split_mask(body: str) -> list[bool]:
    """Boolean mask over ``body`` — True at positions safe to split a
    fixed-form line at (i.e. not inside a string literal). The opening
    quote of a literal is itself marked unsafe, matching legacy
    behavior so the splitter never lands a continuation marker
    immediately before a string."""
    mask = [False] * len(body)
    for i, _ in _iter_outside_strings(body):
        mask[i] = True
    return mask

def reformat_fixed_line(line: str, cont_char: str = '+') -> str:
    # Preprocessor directives (``#if``, ``#include``, ``#define`` ...) are
    # not bound by fixed-form column 72 and must not be split into
    # continuation lines — doing so produces a truncated directive on
    # the first line (e.g. ``#if A || B ||`` dangling) and a second line
    # the preprocessor doesn't understand. Leave them alone regardless
    # of length.
    if line.lstrip().startswith('#'):
        return line
    if len(line) <= 72 or is_comment_line(line) or (len(line) > 6 and line[6:].lstrip().startswith('!')):
        return line
    # If an inline ``!`` comment sits within the first 72 columns, keep
    # the whole line intact — fixed-form Fortran ignores columns past
    # 72, and we must NOT split across the comment (the text after
    # ``!`` would otherwise land on a continuation line as code).
    # Scan for the first ``!`` outside a string literal. If it sits
    # within col 72 we must keep the line intact — splitting it would
    # land the comment text on a continuation chunk.
    bang = _find_inline_bang(line)
    if bang < len(line) and bang < 72:
        return line
    prefix, body = line[:6] if len(line) >= 6 else line.ljust(6), line[6:]
    safe = _build_split_mask(body)
    chunks = []
    while len(body) > 66:
        split_pos = 66
        for i in range(65, max(35, 65 - 30), -1):
            if body[i] in (',', ' ') and safe[i]:
                split_pos = i + 1
                break
        else:
            for i in range(65, 0, -1):
                if safe[i]:
                    split_pos = i
                    break
        chunks.append(body[:split_pos])
        body, safe = body[split_pos:], safe[split_pos:]
    chunks.append(body)
    result_lines = [prefix + chunks[0]]
    for chunk in chunks[1:]: result_lines.append('     ' + cont_char + chunk)
    return '\n'.join(result_lines)


def _strip_inline_comment(text: str) -> str:
    """Strip an inline ``!`` comment from a Fortran code line, respecting
    string literals. The trailing whitespace before the ``!`` is also
    trimmed.

    Used by :func:`_segment_fixed_form_statements` when joining
    continuation lines: an inline ``!`` mid-statement would otherwise
    swallow every continuation that follows once
    :func:`reformat_fixed_line` re-splits the joined string at column
    66 and the ``!`` lands on a chunk that includes content from the
    next continuation line. The comment is irretrievably lost from the
    joined / reformatted output, which is the price of correctness.
    Single-physical-line statements with ``s == joined`` go through the
    no-transform path and keep their comments verbatim from ``lines``.
    """
    bang = _find_inline_bang(text)
    if bang < len(text):
        return text[:bang].rstrip()
    return text


def _segment_fixed_form_statements(
    physical: list[str],
) -> list[tuple[str, list[str], list[str], str]]:
    """Group physical fixed-form lines into logical statements.

    Each entry is ``(kind, lines, terminators, joined)`` where ``kind`` is
    ``'blank' | 'comment' | 'pp' | 'code'``. ``lines`` and ``terminators``
    are aligned slices of ``physical`` (text without newline / the
    original line terminator). For ``'code'`` statements with continuation
    lines, ``joined`` is the head plus each continuation's column-7+
    content concatenated with single spaces — with any inline ``!``
    comments stripped (see :func:`_strip_inline_comment`). This is what
    the per-line transform passes operate on, so paren-walkers can match
    across the physical line break. For other kinds, ``joined`` is the
    head text (with inline comments preserved — single-line statements
    don't risk the swallow-the-next-continuation failure mode).
    """
    out: list[tuple[str, list[str], list[str], str]] = []
    i = 0
    while i < len(physical):
        raw = physical[i]
        if raw.endswith('\r\n'):
            text, term = raw[:-2], '\r\n'
        elif raw.endswith('\n') or raw.endswith('\r'):
            text, term = raw[:-1], raw[-1]
        else:
            text, term = raw, ''
        if not text.strip():
            out.append(('blank', [text], [term], text))
            i += 1
            continue
        if text[0] in 'Cc*!':
            out.append(('comment', [text], [term], text))
            i += 1
            continue
        if text.lstrip().startswith('#'):
            out.append(('pp', [text], [term], text))
            i += 1
            continue
        # Code head — absorb any immediately-following continuation lines.
        lines = [text]
        terms = [term]
        joined = text
        j = i + 1
        while j < len(physical):
            nxt = physical[j]
            if nxt.endswith('\r\n'):
                ntext, nterm = nxt[:-2], '\r\n'
            elif nxt.endswith('\n') or nxt.endswith('\r'):
                ntext, nterm = nxt[:-1], nxt[-1]
            else:
                ntext, nterm = nxt, ''
            if (len(ntext) > 5 and ntext[:1] != '\t' and ntext[:5] == '     '
                    and ntext[5:6] not in (' ', '0', '\t')):
                lines.append(ntext)
                terms.append(nterm)
                # Strip the previous segment's inline ``!`` comment
                # before appending the next continuation, otherwise the
                # comment swallows every following chunk once
                # ``reformat_fixed_line`` re-splits the joined string.
                # On the first continuation, this also strips any inline
                # comment from the head line.
                joined = _strip_inline_comment(joined) + ' ' + ntext[6:]
                j += 1
            else:
                break
        # Strip any inline comment from the trailing continuation too
        # (so the rebuilt single-line / reformatted output isn't truncated
        # at the comment when reformat_fixed_line walks it).
        if len(lines) > 1:
            joined = _strip_inline_comment(joined)
        out.append(('code', lines, terms, joined))
        i = j
    return out


def migrate_fixed_form(source: str, rename_map: dict[str, str], target_mode: TargetMode,
                        source_kind: int | None = None) -> str:
    complex_names = _scan_complex_var_names(source) if not target_mode.is_kind_based else set()
    real_names = _scan_real_var_names(source) if not target_mode.is_kind_based else set()
    source = fix_misdeclared_statement_functions(source, source_kind=source_kind)
    source, removed_known = strip_known_constants_from_decls(source, target_mode)
    source, param_assignments, dropped_p = convert_parameter_stmts(source, target_mode)
    source, data_assignments, dropped_d = convert_data_stmts(source, target_mode)
    removed_known.update(dropped_p)
    removed_known.update(dropped_d)
    source = _dedup_intrinsic_stmts(source, target_mode)
    source = insert_use_multifloats(source, target_mode, extra_lines=param_assignments + data_assignments)
    physical = source.splitlines(keepends=True)
    statements = _segment_fixed_form_statements(physical)
    result = []
    for kind, lines, terms, joined in statements:
        if kind == 'blank':
            result.append(lines[0] + terms[0])
            continue
        if kind == 'pp':
            result.append(lines[0] + terms[0])
            continue
        if kind == 'comment':
            s = replace_routine_names(lines[0], rename_map)
            s = replace_type_decls(s, target_mode, complex_names=complex_names,
                                    source_kind=source_kind)
            result.append(s + terms[0])
            continue
        # 'code' — apply all per-line transforms to the joined logical
        # line so paren-walking passes (e.g. replace_intrinsic_calls)
        # can match calls that span fixed-form continuations. Single-
        # physical-line statements have ``joined == lines[0]`` and pass
        # through identically.
        s = replace_type_decls(joined, target_mode, complex_names=complex_names,
                                source_kind=source_kind)
        if not s:
            continue
        s = replace_standalone_real_complex(s, target_mode, source_kind=source_kind)
        s = replace_literals(s, target_mode, source_kind=source_kind)
        s = replace_intrinsic_calls(
            s, target_mode, real_names=real_names, complex_names=complex_names,
        )
        s = replace_intrinsic_decls(s, target_mode)
        s = replace_generic_conversions(s, target_mode, complex_names=complex_names)
        s = replace_routine_names(s, rename_map)
        s = replace_include_filenames(s, rename_map)
        s = replace_xerbla_strings(s, rename_map)
        s = replace_known_constants(s, target_mode, renames=removed_known)
        s = _rewrite_int_of_complex(s, complex_names)
        s = _rewrite_int_kind_on_real64x2(s, target_mode, real_names=real_names)
        s = _wrap_bare_complex_literals(s, target_mode)
        s = _unwrap_redundant_constructors(s, target_mode, real_names=real_names)
        if len(lines) > 1 and s == joined:
            # Multi-line statement, no transforms applied — emit the
            # original physical lines verbatim to avoid reformat churn.
            for line, term in zip(lines, terms):
                result.append(line + term)
            continue
        s = reformat_fixed_line(s)
        result.append(s + terms[0])

    source = ''.join(result)
    if not target_mode.is_kind_based:
        source = re.sub(r'! !    integer, parameter :: wp = kind\(1\.d0\)',
                        '!    integer, parameter :: wp = kind(1.d0)', source)

    source = _dedup_intrinsic_stmts(source, target_mode)
    source = specialize_use_module(source, target_mode, fixed_form=True)
    return source


_KIND_PARAM_NAMES = r'(?:wp|sp|dp)'
_KIND_PARAM_RE = re.compile(rf'(integer\s*,\s*parameter\s*::\s*{_KIND_PARAM_NAMES}\s*=\s*)(?:kind\s*\(\s*1\.[de]0\s*\)|real(?:32|64|128))', re.IGNORECASE)

def _replace_kind_parameter(line: str, target_mode: TargetMode) -> str:
    if target_mode.is_kind_based: return _KIND_PARAM_RE.sub(rf'\g<1>{target_mode.kind_suffix}', line)
    return ('! ' + line) if _KIND_PARAM_RE.search(line) else line

_ISO_USE_ONLY_RE = re.compile(r'^(?P<lead>\s*)USE\s*,\s*INTRINSIC\s*::\s*ISO_FORTRAN_ENV\s*,\s*ONLY\s*:\s*(?P<names>[^\n!]*?)\s*(?P<tail>!.*)?$', re.IGNORECASE)

def _strip_iso_fortran_env_realN(line: str) -> str:
    m = _ISO_USE_ONLY_RE.match(line)
    if not m: return line
    names = [n.strip() for n in m.group('names').split(',') if n.strip()]
    kept = [n for n in names if not re.fullmatch(r'real(?:32|64|128)', n, re.IGNORECASE)]
    if not kept: return ''
    tail = (' ' + m.group('tail')) if m.group('tail') else ''
    return f'{m.group("lead")}USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: {", ".join(kept)}{tail}'


def rewrite_la_constants_use(source: str, target_mode: TargetMode) -> str:
    """Rewrite ``USE LA_CONSTANTS`` clauses for the chosen target.

    KIND mode (extended precision): the LAPACK la_constants module is
    cloned to ``LA_CONSTANTS_EP`` by the migrator, so we just rename
    the module reference and rename each constant to its EP-prefixed
    equivalent (E*/Y*/Q*/X*).

    Multifloats mode: there is no ``la_constants_mf`` module — instead,
    we rewrite the import to point at the real ``multifloats`` module
    and rename each LAPACK constant (``dzero``, ``dsafmin``, ...) to its
    multifloats equivalent (``MF_ZERO``, ``MF_SAFMIN``, ...). The
    ``wp=>dp`` rename entry is dropped because ``wp`` is no longer
    meaningful once the type becomes ``TYPE(float64x2)``.
    """
    const_renames = _la_constants_rename_map(target_mode)
    lines, result, in_use_stmt = source.split('\n'), [], False
    suffix = target_mode.la_constants_suffix
    target_module_upper = f'LA_CONSTANTS{suffix}'
    target_module_lower = f'la_constants{suffix.lower()}'
    target_xisnan_upper = f'LA_XISNAN{suffix}'
    target_xisnan_lower = f'la_xisnan{suffix.lower()}'

    for line in lines:
        upper = line.upper().lstrip()
        if re.search(r'\bUSE\s+LA_XISNAN\b', upper) and target_xisnan_upper not in upper:
            line = re.sub(
                r'(?i)\bLA_XISNAN\b',
                lambda m: target_xisnan_lower if m.group().islower() else target_xisnan_upper,
                line,
            )
        if re.search(r'\bUSE\s+LA_CONSTANTS\b', upper) and target_module_upper not in upper:
            in_use_stmt = True
            line = re.sub(
                r'(?i)\bLA_CONSTANTS\b',
                lambda m: target_module_lower if m.group().islower() else target_module_upper,
                line,
            )
        if in_use_stmt:
            line = replace_routine_names(line, const_renames)
            if target_mode.is_kind_based:
                # Rename ``wp=>dp`` / ``wp=>sp`` to ``wp=>qp`` (kind16) or
                # ``wp=>ep`` (kind10).  The target kind parameter is the
                # real prefix lowercased + "p".
                target_kp = target_mode.prefix_map['R'].lower() + 'p'
                for kindname in ('dp', 'sp'):
                    line = re.sub(
                        rf'(?i)\b{kindname}\b',
                        lambda m, kp=target_kp: kp if m.group().islower() else kp.upper(),
                        line,
                    )
            else:
                # Strip ``wp=>dp`` (D-source) and ``wp=>sp`` (S-source)
                # entries — both become meaningless after the migrator
                # collapses both halves to float64x2.
                for kindname in ('dp', 'sp'):
                    line = re.sub(rf',\s*wp\s*=>\s*{kindname}\s*,', ',', line, flags=re.IGNORECASE)
                    line = re.sub(rf',\s*wp\s*=>\s*{kindname}\s*(?=[!&]|$)', '', line, flags=re.IGNORECASE)
                    line = re.sub(rf'(ONLY\s*:\s*)wp\s*=>\s*{kindname}\s*,', r'\1', line, flags=re.IGNORECASE)
                    line = re.sub(rf'(ONLY\s*:\s*)wp\s*=>\s*{kindname}\s*(?=[!&]|$)', r'\1', line, flags=re.IGNORECASE)
            if not line.rstrip().endswith('&'):
                in_use_stmt = False
        result.append(line)
    return '\n'.join(result)

_LA_CONSTANTS_REAL_NAMES = (
    'ZERO', 'HALF', 'ONE', 'TWO', 'THREE', 'FOUR', 'EIGHT', 'TEN',
    'PREFIX', 'ULP', 'EPS',
    'SAFMIN', 'SAFMAX', 'SMLNUM', 'BIGNUM',
    'RTMIN', 'RTMAX',
    'TSML', 'TBIG', 'SSML', 'SBIG',
)
_LA_CONSTANTS_COMPLEX_NAMES = ('ZERO', 'HALF', 'ONE', 'PREFIX')


def _la_constants_rename_map(target_mode: TargetMode) -> dict[str, str]:
    """Build a rename map for the RHS of LA_CONSTANTS USE-clause aliases.

    Maps the LAPACK la_constants names ``DZERO``, ``DSAFMIN``, ``ZZERO``,
    etc. to the equivalent names exported by the target la_constants
    auxiliary module:

      KIND=10  → ``ezero``, ``ysafmin`` (la_constants_ep)
      KIND=16  → ``qzero``, ``xsafmin`` (la_constants_ep)
      multifloats → ``ddzero``, ``zzsafmin`` (la_constants_mf)

    Only the prefixed names are mapped — the LHS aliases ``zero``,
    ``half`` etc. are intentionally left untouched so the body of the
    routine continues to reference them through the local alias.
    """
    # Build S/D → target_real_prefix, C/Z → target_complex_prefix map
    # from the target's prefix_map (which maps R→prefix, C→prefix).
    real_pfx = target_mode.prefix_map.get('R', 'Q')
    cmplx_pfx = target_mode.prefix_map.get('C', 'X')
    pmap = {'S': real_pfx, 'D': real_pfx, 'C': cmplx_pfx, 'Z': cmplx_pfx}

    renames: dict[str, str] = {}
    for p in ('S', 'D'):
        for base in _LA_CONSTANTS_REAL_NAMES:
            renames[p + base] = pmap[p] + base
    for p in ('C', 'Z'):
        for base in _LA_CONSTANTS_COMPLEX_NAMES:
            renames[p + base] = pmap[p] + base
    return renames


def migrate_free_form(source: str, rename_map: dict[str, str], target_mode: TargetMode,
                       source_kind: int | None = None) -> str:
    complex_names = _scan_complex_var_names(source) if not target_mode.is_kind_based else set()
    real_names = _scan_real_var_names(source) if not target_mode.is_kind_based else set()
    source = rewrite_la_constants_use(source, target_mode)
    source = fix_misdeclared_statement_functions(source, source_kind=source_kind)
    source, removed_known = strip_known_constants_from_decls(source, target_mode)
    if not target_mode.is_kind_based:
        lines_tmp = source.splitlines()
        res_tmp = []
        in_comment_block = False
        # Names of free-form file-scope PARAMETER declarations that are
        # supplied as MF_* constants by the multifloats module. RTMIN and
        # RTMAX are intentionally excluded: in several LAPACK routines
        # they are local variables computed at runtime, not PARAMETERs.
        # The mapping mirrors la_constants_map but with explicit
        # multifloats target names for free-form Pattern A files (those
        # use ``USE multifloats`` directly, not ``USE LA_CONSTANTS_MF``).
        nuke_renames = {
            'zero': 'DD_ZERO', 'one': 'DD_ONE', 'czero': 'DD_ZERO',
            'safmin': 'DD_SAFMIN', 'safmax': 'DD_SAFMAX',
            'tsml': 'DD_TSML', 'tbig': 'DD_TBIG',
            'ssml': 'DD_SSML', 'sbig': 'DD_SBIG',
            # dnrm2.f90's local ``maxN = huge(0.0_wp)`` is equivalent to
            # DD_SAFMAX (which is itself defined as huge(0.0_dp) packed
            # into float64x2's high limb).
            'maxn': 'DD_SAFMAX',
        }
        nuke_names = set(nuke_renames.keys())

        for line in lines_tmp:
            stripped = line.strip().lower()
            is_decl_start = re.match(r'^\s*(?:real|complex|integer|type|parameter).*?::', line, re.IGNORECASE) or \
                            re.match(r'^\s*parameter\s*\(', line, re.IGNORECASE)

            contains_nuke = False
            matched_names: list[str] = []
            for n in nuke_names:
                if re.search(rf'\b{n}\b', stripped):
                    contains_nuke = True
                    matched_names.append(n)

            if not in_comment_block and is_decl_start and contains_nuke:
                res_tmp.append('! ' + line)
                if line.rstrip().endswith('&'): in_comment_block = True
                for n in matched_names:
                    removed_known[n.upper()] = nuke_renames[n]
            elif in_comment_block:
                res_tmp.append('! ' + line)
                if not line.rstrip().endswith('&'): in_comment_block = False
                for n in matched_names:
                    removed_known[n.upper()] = nuke_renames[n]
            else: res_tmp.append(line)
        source = '\n'.join(res_tmp)

    source, param_assignments, dropped_p = convert_parameter_stmts(source, target_mode)
    source, data_assignments, dropped_d = convert_data_stmts(source, target_mode)
    removed_known.update(dropped_p)
    removed_known.update(dropped_d)

    source = _dedup_intrinsic_stmts(source, target_mode)
    source = insert_use_multifloats(source, target_mode, extra_lines=param_assignments + data_assignments)
    lines = source.splitlines(keepends=True)
    result = []
    for line in lines:
        stripped = line.rstrip('\n\r')
        nl = '\n' if line.endswith('\n') else ''
        stripped = _replace_kind_parameter(stripped, target_mode)
        if not stripped.lstrip().startswith('!'):
            stripped = _strip_iso_fortran_env_realN(stripped)
            if not stripped: continue
            stripped = replace_intrinsic_calls(
                stripped, target_mode, real_names=real_names,
                complex_names=complex_names,
            )
            stripped = replace_intrinsic_decls(stripped, target_mode)
            stripped = replace_generic_conversions(
                stripped, target_mode, complex_names=complex_names,
            )
        stripped = replace_routine_names(stripped, rename_map)
        stripped = replace_include_filenames(stripped, rename_map)
        if stripped.lstrip().startswith('!'):
            stripped = replace_type_decls(stripped, target_mode, complex_names=complex_names,
                                           source_kind=source_kind)
        else:
            if not target_mode.is_kind_based:
                stripped = re.sub(r'REAL\s*\(\s*(?:KIND\s*=\s*)?' + _KIND_PARAM_NAMES + r'\s*\)', target_mode.real_type, stripped, flags=re.IGNORECASE)
                stripped = re.sub(r'COMPLEX\s*\(\s*(?:KIND\s*=\s*)?' + _KIND_PARAM_NAMES + r'\s*\)', target_mode.complex_type, stripped, flags=re.IGNORECASE)
            # Rewrite floating-point literals (D/E exponents and the
            # ``1.23_wp`` form). Without this, ``( 1.0_WP, 0.0_WP )``
            # complex constants in modern LAPACK files (DGEDMD/DGEDMDQ)
            # would be left untouched and gfortran would reject the
            # KIND parameter once ``wp`` itself has been stripped.
            if not target_mode.is_kind_based:
                stripped = replace_literals(stripped, target_mode, source_kind=source_kind)
            stripped = replace_known_constants(stripped, target_mode, renames=removed_known)
            stripped = _rewrite_int_of_complex(stripped, complex_names)
            stripped = _rewrite_int_kind_on_real64x2(stripped, target_mode, real_names=real_names)
            stripped = _wrap_bare_complex_literals(stripped, target_mode)
            stripped = _unwrap_redundant_constructors(stripped, target_mode, real_names=real_names)
        result.append(stripped + nl)

    source = ''.join(result)
    if not target_mode.is_kind_based:
        source = re.sub(r'(?i)!\s*!\s*integer\s*,\s*parameter\s*::\s*wp\s*=',
                        '!    integer, parameter :: wp =', source)

    source = _dedup_intrinsic_stmts(source, target_mode)
    source = specialize_use_module(source, target_mode, fixed_form=False)
    return source


def target_filename(name: str, rename_map: dict[str, str],
                    target_mode: TargetMode | None = None) -> str:
    stem, ext = Path(name).stem, Path(name).suffix
    if stem.upper() in rename_map:
        new = rename_map[stem.upper()]
        return (new.upper() if stem.isupper() else (new.lower() if stem.islower() else new.upper())) + ext
    # Fallback for libraries whose filenames encode the arithmetic only in
    # the first character (e.g. MUMPS: ``dana_aux.F``, ``zfac_driver.F``).
    # When the stem is not a routine name, translate the leading s/d/c/z
    # into the target's arithmetic letter via the family prefix map.
    if target_mode is not None and stem and stem[0].upper() in ('S', 'D', 'C', 'Z'):
        from .prefix_classifier import CHAR_TYPE
        family = CHAR_TYPE[stem[0].upper()]         # 'R' or 'C'
        new_char = target_mode.prefix_map.get(family)
        if new_char:
            first = new_char if stem[0].isupper() else new_char.lower()
            return first + stem[1:] + ext
    return name


_KK_SENTINEL = '__KEEPKIND_DP__'
_KK_DBLE_SENTINEL = '__KEEPKIND_DBLE__'
_KK_DCMPLX_SENTINEL = '__KEEPKIND_DCMPLX__'


def _apply_keep_kind_sentinel(source: str, keep_kind_lines: frozenset[int]) -> str:
    """Replace DP-defining tokens with non-matching sentinels on each
    1-based line number in ``keep_kind_lines``. Protects the line from
    every migrator regex that rewrites those tokens; restored by
    :func:`_restore_keep_kind_sentinel` on the migrated output.

    Protected tokens: ``DOUBLE PRECISION`` (type declaration),
    ``dble(`` (convert-to-DP intrinsic), ``dcmplx(`` (convert-to-DC
    intrinsic). The last two are protected on call-site lines so that
    callers of verbatim (copy_files) DP-stable routines keep passing
    DP values instead of being rewritten to ``REAL(x, KIND=16)``.
    """
    import re as _re
    lines = source.splitlines(keepends=True)
    _dp = _re.compile(r'DOUBLE\s+PRECISION', _re.IGNORECASE)
    _dble = _re.compile(r'\bdble\s*\(', _re.IGNORECASE)
    _dcmplx = _re.compile(r'\bdcmplx\s*\(', _re.IGNORECASE)
    for ln in keep_kind_lines:
        if 1 <= ln <= len(lines):
            t = _dp.sub(_KK_SENTINEL, lines[ln - 1])
            t = _dble.sub(_KK_DBLE_SENTINEL + '(', t)
            t = _dcmplx.sub(_KK_DCMPLX_SENTINEL + '(', t)
            lines[ln - 1] = t
    return ''.join(lines)


def _restore_keep_kind_sentinel(source: str) -> str:
    source = source.replace(_KK_SENTINEL, 'DOUBLE PRECISION')
    source = source.replace(_KK_DBLE_SENTINEL, 'dble')
    source = source.replace(_KK_DCMPLX_SENTINEL, 'dcmplx')
    return source


def _strip_roundup_lwork(source: str, target_mode) -> str:
    """At kind ≥ 8 the SROUNDUP_LWORK / DROUNDUP_LWORK wrapper is dead
    code: it guards a float→int round-trip bug that only bites when
    LWORK exceeds the float mantissa (2**24 at single precision). Any
    target with mantissa ≥ 32 bits (kind8 and above, plus multifloats)
    represents 32-bit INTEGER exactly, so the IF branch is unreachable
    and the function is the identity.

    Strip the call wrapper and its declarations so:
      - the S-derived re-migration matches the D-derived canonical
        (which already uses bare ``WORK(1)=LWKOPT`` upstream);
      - the dead ``qroundup_lwork.f`` symbol drops out of the link.

    Skipped at kind4 (single precision is exactly where the wrapper
    matters). Operates on the post-rename text, matching any single-
    letter prefix so it works regardless of target prefix (Q at
    kind16, M at multifloats, etc.).
    """
    if target_mode.is_kind_based and (target_mode.kind_suffix or 0) < 8:
        return source
    import re as _re
    name_re = _re.compile(r'\b([A-Z])ROUNDUP_LWORK\b')
    if not name_re.search(source):
        return source

    # 1) Strip call wrappers: <PREFIX>ROUNDUP_LWORK( <expr> )  ->  <expr>
    out = []
    i = 0
    n = len(source)
    while i < n:
        m = name_re.search(source, i)
        if not m:
            out.append(source[i:])
            break
        out.append(source[i:m.start()])
        # Skip the function definition itself (FUNCTION QROUNDUP_LWORK(...))
        prefix = source[max(0, m.start() - 16):m.start()].upper()
        if _re.search(r'\bFUNCTION\s*$', prefix):
            out.append(m.group(0))
            i = m.end()
            continue
        # Find following '(' (possibly with whitespace)
        j = m.end()
        while j < n and source[j] in ' \t':
            j += 1
        if j >= n or source[j] != '(':
            # Not a call (e.g. EXTERNAL listing) — keep token, handled below.
            out.append(m.group(0))
            i = m.end()
            continue
        # Walk balanced parens
        depth = 0
        k = j
        while k < n:
            c = source[k]
            if c == '(':
                depth += 1
            elif c == ')':
                depth -= 1
                if depth == 0:
                    break
            k += 1
        if k >= n:
            # Unbalanced — bail, keep original
            out.append(source[m.start():])
            break
        out.append(source[j + 1:k])  # the inner expression
        i = k + 1
    source = ''.join(out)

    # 2) Drop the now-orphan declarations. Logical-line-aware: joins
    # fixed-form ($/+/&) and free-form (& at EOL) continuations, strips
    # the dead token from any REAL/EXTERNAL list, then re-emits — also
    # dropping lines whose entire list became empty.
    physical = source.splitlines(keepends=True)
    # Identify continuation runs.
    def is_cont_fixed(line: str) -> bool:
        # Fixed-form: column 6 has a non-space, non-zero character.
        # We accept whitespace-leading lines starting with $/+/& at col 6.
        return len(line) > 5 and line[5] not in (' ', '0', '\t', '\n')
    def is_cont_free(prev: str) -> bool:
        # Free-form: previous line ended (before any trailing comment) in '&'.
        s = prev.rstrip('\n').rstrip()
        # strip a trailing '! ...' comment for the heuristic
        if '!' in s:
            s = s.split('!', 1)[0].rstrip()
        return s.endswith('&')

    out_lines: list[str] = []
    i = 0
    while i < len(physical):
        line = physical[i]
        # Comment lines pass through unchanged.
        if line[:1] in ('C', 'c', '*', '!'):
            out_lines.append(line); i += 1; continue
        # Build logical-line group [i, j)
        j = i + 1
        while j < len(physical):
            nxt = physical[j]
            if nxt[:1] in ('C', 'c', '*', '!'):
                # In fixed form, comment lines may interleave continuations;
                # but if comment line, stop the group here.
                break
            if is_cont_fixed(nxt) or is_cont_free(physical[j - 1]):
                j += 1
            else:
                break
        group = physical[i:j]
        # Reconstruct logical content (strip continuation glue for matching).
        logical = group[0]
        for k in range(1, len(group)):
            ln = group[k]
            # strip continuation marker (col 6 in fixed form) — we just
            # blank out positions 0..5 of continuation lines.
            if is_cont_fixed(ln):
                logical = logical.rstrip('\n').rstrip() + ' ' + ln[6:].lstrip()
            else:
                # free-form: drop trailing & on prior, drop leading & on this
                logical = logical.rstrip('\n').rstrip()
                if logical.endswith('&'):
                    logical = logical[:-1].rstrip()
                cur = ln.lstrip()
                if cur.startswith('&'):
                    cur = cur[1:].lstrip()
                logical = logical + ' ' + cur
        # Match on the joined logical line.
        joined = logical.rstrip('\n').rstrip()

        # REAL(KIND=N) <list>   — strip ROUNDUP token from list
        m_real = _re.match(
            r'^(\s*REAL(?:\s*\(\s*KIND\s*=\s*\d+\s*\))?\s+)(.*)$',
            joined, _re.IGNORECASE)
        m_ext = _re.match(r'^(\s*EXTERNAL\s+)(.*)$',
                          joined, _re.IGNORECASE)
        # Modern F90 attribute-list form: ``REAL, EXTERNAL :: NAME``
        # and ``DOUBLE PRECISION, EXTERNAL :: NAME``. Migrator-level
        # strip needs to handle both to drop the orphan ROUNDUP_LWORK
        # declaration after the call sites are stripped.
        m_attr = _re.match(
            r'^(\s*(?:REAL|DOUBLE\s+PRECISION)'
            r'(?:\s*\(\s*KIND\s*=\s*\d+\s*\))?'
            r'\s*,\s*EXTERNAL\s*::\s*)(.*)$',
            joined, _re.IGNORECASE)
        is_attr = bool(m_attr)
        is_real = m_real and not m_ext and not is_attr
        is_ext  = bool(m_ext) and not is_attr
        if is_ext or is_real or is_attr:
            m = m_attr if is_attr else (m_ext if is_ext else m_real)
            head, names = m.group(1), m.group(2)
            tokens = [t.strip() for t in names.split(',')]
            keep = [t for t in tokens
                    if t and not _re.fullmatch(
                        r'[A-Z]ROUNDUP_LWORK', t, _re.IGNORECASE)]
            if len(keep) == len(tokens):
                # Token not present — emit group unchanged.
                out_lines.extend(group)
            elif not keep:
                # Drop entire logical-line group.
                pass
            else:
                # Re-emit as a single line (give up on preserving the
                # original continuation layout — the result is still
                # legal Fortran in either fixed or free form, since
                # EXTERNAL/REAL declaration lines are not column-7
                # constrained).
                # Keep leading whitespace from original line for cosmetics.
                lead = _re.match(r'^(\s*)', group[0]).group(1)
                new_logical = head + ', '.join(keep) + '\n'
                # Preserve original leading indent from first phys line.
                # head may not have leading whitespace; ensure it does:
                if not new_logical.startswith(lead):
                    new_logical = lead + new_logical.lstrip()
                out_lines.append(new_logical)
        else:
            out_lines.extend(group)
        i = j
    return ''.join(out_lines)


# Fortran-side MPI datatype name rewriter. In an `s*`/`c*` source MUMPS
# uses ``MPI_REAL`` and ``MPI_COMPLEX``; in `d*`/`z*` it uses
# ``MPI_DOUBLE_PRECISION`` and ``MPI_DOUBLE_COMPLEX``. After migration
# both halves must refer to the target's wider datatype (e.g.
# ``MPI_REAL16``/``MPI_COMPLEX32`` for kind16). Without this pass the two
# halves' outputs disagree on every MPI call, even though they are
# semantically identical. Word boundaries keep the rewrite idempotent —
# ``MPI_REAL16`` already has no ``\b`` after ``REAL`` so it is not
# rematched.
_MPI_DOUBLE_COMPLEX_RE = re.compile(r'\bMPI_DOUBLE_COMPLEX\b')
_MPI_DOUBLE_PRECISION_RE = re.compile(r'\bMPI_DOUBLE_PRECISION\b')
_MPI_COMPLEX_RE = re.compile(r'\bMPI_COMPLEX\b')
_MPI_REAL_RE = re.compile(r'\bMPI_REAL\b')


# Fortran ``INCLUDE 'foo.h'`` whose filename starts with an arithmetic
# letter. MUMPS uses this for per-arithmetic C-interop headers
# (``INCLUDE 'dmumps_struc.h'`` in dmumps_struc_def.F). After migration
# the file on disk is renamed by target_filename's first-char fallback
# (dmumps_struc.h → qmumps_struc.h); the INCLUDE string must be
# rewritten to match, otherwise the compiler can't find it.
_INCLUDE_PREFIXED_H_RE = re.compile(
    r"(INCLUDE\s*['\"])([SsDdCcZz])([\w]*\.h)(['\"])",
    re.IGNORECASE,
)


def _rewrite_prefixed_includes(source: str, target_mode: TargetMode) -> str:
    from .prefix_classifier import CHAR_TYPE
    pmap = target_mode.prefix_map

    def sub(m: re.Match) -> str:
        head, prefix, rest, tail = m.group(1), m.group(2), m.group(3), m.group(4)
        family = CHAR_TYPE.get(prefix.upper())
        new = pmap.get(family) if family else None
        if not new:
            return m.group(0)
        first = new if prefix.isupper() else new.lower()
        return head + first + rest + tail

    return _INCLUDE_PREFIXED_H_RE.sub(sub, source)


def _rewrite_mpi_datatypes(source: str, target_mode: TargetMode) -> str:
    mpi_real = target_mode.c_mpi_real
    mpi_complex = target_mode.c_mpi_complex
    if not mpi_real and not mpi_complex:
        return source
    if mpi_complex:
        source = _MPI_DOUBLE_COMPLEX_RE.sub(mpi_complex, source)
        source = _MPI_COMPLEX_RE.sub(mpi_complex, source)
    if mpi_real:
        source = _MPI_DOUBLE_PRECISION_RE.sub(mpi_real, source)
        source = _MPI_REAL_RE.sub(mpi_real, source)
    return source


# Token-context MPI reduction-op rewriter for multifloats. The target
# ships user-defined ops (MPI_DD_SUM / MPI_DD_AMX / MPI_DD_AMN, plus
# the ZZ-prefixed complex variants) that operate on the user-defined
# float64x2 / complex64x2 datatypes; the stock MPI_SUM / MPI_MAX /
# MPI_MIN are undefined for those datatypes and the runtime aborts.
# The C migrator handles this by file-arithmetic context (real /
# complex sibling files use separate rule lists) but Fortran sources
# from MUMPS mix integer reductions (MPI_INTEGER + MPI_SUM, which
# MUST stay MPI_SUM) with floating-point reductions in the same
# translation unit, so a blunt file-prefix rewrite would corrupt the
# integer calls. Instead we match each MPI_(I?All)Reduce[_scatter]
# call as a whole, run *after* _rewrite_mpi_datatypes, and only
# rewrite the op inside a call whose argument list contains the
# rewritten floating datatype token (e.g. MPI_FLOAT64X2). Integer-
# typed calls are left alone by construction. For KIND targets the
# c_mpi_(sum|max|min)_* fields expand to 'MPI_(SUM|MAX|MIN)' or are
# unset, and the pass is a no-op.
_MPI_REDUCE_CALL_RE = re.compile(
    r'\bMPI_(?:I?ALL)?REDUCE(?:_SCATTER)?\b\s*'
    r'\((?:[^()]|\([^()]*\))*\)',
    re.IGNORECASE,
)
_MPI_OP_RE = {
    'MPI_SUM': re.compile(r'\bMPI_SUM\b'),
    'MPI_MAX': re.compile(r'\bMPI_MAX\b'),
    'MPI_MIN': re.compile(r'\bMPI_MIN\b'),
}


def _rewrite_mpi_sum(source: str, target_mode: TargetMode) -> str:
    real_tok = target_mode.c_mpi_real
    complex_tok = target_mode.c_mpi_complex
    rules: list[tuple[re.Pattern, str | None, str | None]] = []
    for stock, mf_real, mf_complex in (
        ('MPI_SUM', target_mode.c_mpi_sum_real, target_mode.c_mpi_sum_complex),
        ('MPI_MAX', target_mode.c_mpi_max_real, target_mode.c_mpi_max_complex),
        ('MPI_MIN', target_mode.c_mpi_min_real, target_mode.c_mpi_min_complex),
    ):
        r = mf_real if mf_real and mf_real != stock else None
        c = mf_complex if mf_complex and mf_complex != stock else None
        if r or c:
            rules.append((_MPI_OP_RE[stock], r, c))
    if not rules:
        return source

    def rewrite(m: re.Match) -> str:
        call = m.group(0)
        is_complex = bool(complex_tok) and complex_tok in call
        is_real = bool(real_tok) and real_tok in call
        if not (is_complex or is_real):
            return call
        for pat, r, c in rules:
            if is_complex and c:
                call = pat.sub(c, call)
            elif is_real and r:
                call = pat.sub(r, call)
        return call

    return _MPI_REDUCE_CALL_RE.sub(rewrite, source)


# Multifloats MPI handle tokens emitted by ``_rewrite_mpi_datatypes`` for
# the multifloats target (``MPI_DOUBLE_PRECISION`` → ``MPI_FLOAT64X2``,
# ``MPI_DOUBLE_COMPLEX`` → ``MPI_COMPLEX64X2``) plus the real + complex
# reduction op handles the bridge exposes. Migrated MUMPS Fortran calls
# MPI directly with these names; without a Fortran module declaring
# them, gfortran rejects every site with "has no IMPLICIT type". The
# ``multifloats_mpi_f`` module
# (external/multifloats-mpi/multifloats_mpi_f.f90) bind-to-C declares
# them as default-INTEGER variables populated at MPI init time.
_MF_MPI_TOKENS = (
    'MPI_FLOAT64X2', 'MPI_COMPLEX64X2',
    'MPI_DD_SUM', 'MPI_DD_AMX', 'MPI_DD_AMN',
    'MPI_ZZ_SUM', 'MPI_ZZ_AMX', 'MPI_ZZ_AMN',
)
_MF_MPI_ANY_RE = re.compile(r'\b(?:' + '|'.join(_MF_MPI_TOKENS) + r')\b')


def insert_use_multifloats_mpi_f(source: str, target_mode: TargetMode) -> str:
    """For multifloats targets only: inject ``USE multifloats_mpi_f``
    after each procedure header in sources that reference custom MPI
    handles. Runs after ``_rewrite_mpi_datatypes`` so the rewritten
    tokens are visible. The module exports only ``MPI_*``-prefixed
    names plus ``multifloats_mpi_init``, so no ONLY clause is needed
    — collision risk is nil and the bare USE keeps fixed-form lines
    well under column 72 even when many tokens are involved."""
    if target_mode.is_kind_based:
        return source
    if not _MF_MPI_ANY_RE.search(source):
        return source

    lines = source.splitlines(keepends=True)
    proc_header_re = _PROC_HEADER_RE
    end_proc_re = _END_PROC_RE

    result: list[str] = []
    i = 0
    while i < len(lines):
        line = lines[i]
        result.append(line)
        m = proc_header_re.match(line)
        if not m:
            i += 1
            continue
        # Walk past the procedure header (continuations + CPP blocks),
        # mirroring insert_use_multifloats so paren-balanced argument
        # lists that span ``#if defined(metis)`` blocks are handled.
        j = i + 1
        paren_depth = _count_open_parens(line)
        prev_has_amp = line.rstrip().rstrip('\n').endswith('&')
        while j < len(lines):
            next_line = lines[j]
            if next_line.lstrip().startswith('#'):
                result.append(next_line)
                j += 1
                continue
            if (is_continuation_line(next_line) or prev_has_amp
                    or paren_depth > 0):
                result.append(next_line)
                prev_has_amp = next_line.rstrip().rstrip('\n').endswith('&')
                paren_depth += _count_open_parens(next_line)
                j += 1
            else:
                break
        # Find this procedure's matching END via depth tracking. Many
        # MUMPS files put an INTERFACE block inside the outer SUBROUTINE
        # whose internal SUBROUTINE signatures look syntactically like
        # procedure headers; their END SUBROUTINE would falsely
        # terminate a naive linear scan, missing the MPI_FLOAT64X2
        # sites further down in the outer body. Counting nested opens
        # vs closes lets us land on the actual matching END at depth 0.
        depth = 1
        k = j
        while k < len(lines):
            if proc_header_re.match(lines[k]):
                depth += 1
            elif end_proc_re.match(lines[k]):
                depth -= 1
                if depth == 0:
                    break
            k += 1
        body_text = ''.join(lines[j:k + 1])
        if _MF_MPI_ANY_RE.search(body_text):
            already_has = any('USE multifloats_mpi_f' in lines[kk]
                              for kk in range(j, min(j + 30, len(lines))))
            if not already_has:
                indent = m.group(1)
                ind = indent if indent.strip() else '      '
                result.append(f"{ind}USE multifloats_mpi_f\n")
        i = j
    return ''.join(result)


def _source_kind_from_filename(name: str) -> int | None:
    """Return the source half's native kind (4 or 8) or None.

    Strips an optional ScaLAPACK ``p`` prefix; the next character
    determines the precision letter:
      ``s``/``c`` (single, single-complex) → kind 4
      ``d``/``z`` (double, double-complex) → kind 8
    Other prefixes (precision-independent files like ``lsame.f``,
    ``mumps_*.F``) return None — promote every kind, current default.

    Used by :func:`migrate_file_to_string` to scope ``replace_type_decls``
    so a kind4 source half preserves kind8 references and vice versa
    (rule (a) of the per-half promotion convention).
    """
    stem = Path(name).stem.lower()
    if not stem:
        return None
    # ScaLAPACK 2-letter (P + S/D/C/Z)
    if stem[0] == 'p' and len(stem) > 1 and stem[1] in 'sdcz':
        letter = stem[1]
    elif stem[0] in 'sdcz':
        letter = stem[0]
    else:
        return None
    return 4 if letter in ('s', 'c') else 8


def migrate_file_to_string(src_path: Path, rename_map: dict[str, str], target_mode: TargetMode, parser: str | None = None, parser_cmd: str | None = None, keep_kind_lines: frozenset[int] | None = None) -> tuple[str, str] | None:
    ext, source = src_path.suffix.lower(), src_path.read_text(errors='replace')
    facts = None
    if parser == 'flang':
        from .flang_parser import scan_file as flang_scan, find_flang
        cmd = parser_cmd or find_flang()
        if cmd: facts = flang_scan(src_path, cmd)
    elif parser == 'gfortran':
        from .gfortran_parser import scan_file as gfortran_scan, find_gfortran
        cmd = parser_cmd or find_gfortran()
        if cmd: facts = gfortran_scan(src_path, cmd)

    if keep_kind_lines:
        source = _apply_keep_kind_sentinel(source, keep_kind_lines)

    source_kind = _source_kind_from_filename(src_path.name)

    if facts is not None: migrated = _migrate_with_flang(source, ext, rename_map, target_mode, facts, source_kind=source_kind)
    elif ext in ('.f', '.for', '.h'): migrated = migrate_fixed_form(source, rename_map, target_mode, source_kind=source_kind)
    elif ext in ('.f90', '.f95', '.F90'): migrated = migrate_free_form(source, rename_map, target_mode, source_kind=source_kind)
    else: return None

    if keep_kind_lines:
        migrated = _restore_keep_kind_sentinel(migrated)

    migrated = _rewrite_mpi_datatypes(migrated, target_mode)
    migrated = _rewrite_mpi_sum(migrated, target_mode)
    migrated = insert_use_multifloats_mpi_f(migrated, target_mode)
    migrated = _rewrite_prefixed_includes(migrated, target_mode)
    migrated = _strip_roundup_lwork(migrated, target_mode)

    out_name = target_filename(src_path.name, rename_map, target_mode)

    # Type-conversion intrinsics frequently drift asymmetrically between
    # co-family halves: D/Z-half upstream tends to declare ``REAL`` /
    # ``CMPLX`` / ``DBLE`` / ``DIMAG`` / ``DCMPLX`` / ``DCONJG`` in its
    # INTRINSIC list (used to convert ``DOUBLE COMPLEX`` ↔ ``DOUBLE
    # PRECISION``), while S/C-half doesn't need them (kind4 default
    # handles the same conversions implicitly). After kind16 migration
    # both halves use the kind-promoted versions; the asymmetric
    # INTRINSIC declarations remain as cosmetic text drift. Strip the
    # type-conversion subset on every target to converge.
    type_conv_intrinsics = {
        'REAL', 'CMPLX', 'DCMPLX', 'DBLE',
        'AIMAG', 'CONJG', 'DCONJG', 'DIMAG',
    }
    if not target_mode.is_kind_based:
        # Multifloats also strips the full generic-overload set so
        # INTRINSIC declarations of names that ``USE multifloats``
        # provides as generic interfaces don't clash (gfortran: "Cannot
        # change attributes of USE-associated symbol").
        strip_set = type_conv_intrinsics | {
            # Unary real -> real
            'ABS', 'SQRT', 'SIN', 'COS', 'TAN', 'EXP', 'LOG', 'LOG10',
            'ATAN', 'ASIN', 'ACOS', 'AINT', 'ANINT', 'SINH', 'COSH',
            'TANH', 'ASINH', 'ACOSH', 'ATANH', 'ERF', 'ERFC',
            'ERFC_SCALED', 'GAMMA', 'LOG_GAMMA', 'BESSEL_J0',
            'BESSEL_J1', 'BESSEL_Y0', 'BESSEL_Y1', 'FRACTION',
            'RRSPACING', 'SPACING', 'EPSILON', 'HUGE', 'TINY',
            # Binary real
            'SIGN', 'MOD', 'ATAN2', 'DIM', 'MODULO', 'HYPOT', 'NEAREST',
            # Variadic
            'MAX', 'MIN',
            # Extra type-conversion
            'INT', 'NINT', 'CEILING', 'FLOOR',
        }
    else:
        strip_set = type_conv_intrinsics

    def clean_intrinsic(m):
        indent, sep, funcs_str, newline = m.group(1), m.group(2), m.group(3), m.group(4)
        # Skip continuation lines: trailing ``,`` means more on the next
        # physical line, and stripping a name on this line could leave
        # the continuation joined into a fused token (e.g.
        # ``MAXEXPONENT,\n     $MINEXPONENT`` → ``MAXEXPONENTMINEXPONENT``).
        # The full-list case is rare enough that conservative single-line
        # handling captures the asymmetry we care about (S/C vs D/Z
        # halves differ in the REAL/CMPLX type-conversion entries, which
        # are always on the first line) without touching multi-line
        # declarations.
        if funcs_str.rstrip().endswith(','):
            return m.group(0)
        kept = [f.strip() for f in funcs_str.split(',')
                if f.strip() and f.strip().upper() not in strip_set]
        return f"{indent}INTRINSIC{sep}{', '.join(kept)}{newline}" if kept else ""

    migrated = re.sub(
        r'(?im)^([ \t]*)INTRINSIC(\s*::\s*|\s+)([A-Za-z0-9_,\s]+?)(\r?\n|$)',
        clean_intrinsic, migrated,
    )
    return out_name, migrated


def _migrate_with_flang(source: str, ext: str, rename_map: dict[str, str], target_mode: TargetMode, facts,
                         source_kind: int | None = None) -> str:
    # Also include USE-statement module names + every variable's type
    # spec (where the type is precision-prefixed, e.g.
    # ``TYPE(DMUMPS_INTR_STRUC)``). The gfortran parser doesn't surface
    # type-of-variable as a call/external_name, so without this extra
    # source the rename_map gets filtered down to nothing for files
    # whose only precision-prefixed reference is via ``USE`` or
    # ``TYPE(...)``.
    use_names = {u.module_name for u in (facts.use_stmt_ranges or [])}
    type_refs: set[str] = set()
    if facts.variable_types:
        # variable_types maps name -> ts; the ts may itself be a
        # precision-prefixed derived type. Add any token of the form
        # ``DMUMPS_*`` from the source as a candidate (cheap regex).
        # This catches TYPE(DMUMPS_FOO) :: var declarations the parser
        # doesn't lift into a routine_def.
        for token in re.findall(r'\b([SDCZ]MUMPS_[A-Z0-9_]+)\b', source.upper()):
            type_refs.add(token)
    file_names = (
        {rd.name for rd in facts.routine_defs}
        | {cs.name for cs in facts.call_sites}
        | set(facts.external_names)
        | use_names
        | type_refs
    )
    file_rename_map = {k: v for k, v in rename_map.items() if k in file_names}
    has_float_types = any(td.type_spec in ('DoublePrecision', 'Real', 'Complex') for td in facts.type_decls)
    if not has_float_types:
        # gfortran's dump puts DOUBLE PRECISION components inside a
        # TYPE body in the parent type's symbol entry, not as
        # individual type_decls. Without this fall-back, files that
        # only declare DOUBLE PRECISION inside derived-type bodies
        # (e.g. MUMPS's *_intr_types.F) skip the promotion pass and
        # the migrated module still references DOUBLE PRECISION at
        # the field level — incompatible with callers that expect
        # the promoted REAL(KIND=...) target type.
        # The various ways an FP type can appear in a source file. Keep
        # this list aligned with `replace_type_decls` -- if it can
        # rewrite the form, this fallback should trigger so the rewrite
        # actually runs. Catches DOUBLE PRECISION, DOUBLE COMPLEX,
        # COMPLEX*N, REAL*N, COMPLEX(...), REAL(...) declarations.
        _fp_patterns = [
            r'\bDOUBLE\s+PRECISION\b',
            r'\bDOUBLE\s+COMPLEX\b',
            r'\bCOMPLEX\s*\*\s*\d+',
            r'\bREAL\s*\*\s*\d+',
            r'\bCOMPLEX\s*\(\s*(?:KIND\s*=\s*)?\d+\s*\)',
            r'\bREAL\s*\(\s*(?:KIND\s*=\s*)?\d+\s*\)',
        ]
        if any(re.search(p, source, re.IGNORECASE) for p in _fp_patterns):
            has_float_types = True
    has_real_literals = bool(facts.real_literals)
    if not has_real_literals:
        # gfortran's dump only emits `value:` lines for symbols whose
        # initializer is a single bare literal — `parameter :: ZERO =
        # 0.0D+0` lifts the literal into the symbol's `value:` and is
        # captured, but `parameter :: CZERO = (0.0E+0, 0.0E+0)` (a
        # complex constructor expression) is NOT captured. Without
        # this fallback, files whose only real literals live inside
        # complex/array PARAMETER initializers skip the literal-
        # promotion pass and `0.0E+0` survives the migration as a
        # single-precision REAL(KIND=4) value rather than the target
        # KIND. Catches `1.0D+0` / `1.0e0` / `0.0` style literals.
        if re.search(r'(?<![\[\d])\d+\.\d*[DdEe][+-]?\d+', source):
            has_real_literals = True
    if ext in ('.f90', '.f95', '.F90'): return _migrate_free_form_flang(source, file_rename_map, target_mode, has_float_types, source_kind=source_kind)
    return _migrate_fixed_form_flang(source, file_rename_map, target_mode, has_float_types, has_real_literals, source_kind=source_kind)


def _migrate_fixed_form_flang(source: str, rename_map: dict[str, str], target_mode: TargetMode, has_float_types: bool, has_real_literals: bool, source_kind: int | None = None) -> str:
    complex_names = _scan_complex_var_names(source) if not target_mode.is_kind_based else set()
    real_names = _scan_real_var_names(source) if not target_mode.is_kind_based else set()
    source = fix_misdeclared_statement_functions(source, source_kind=source_kind)
    source, removed_known = strip_known_constants_from_decls(source, target_mode)
    source, param_assignments, dropped_p = convert_parameter_stmts(source, target_mode)
    source, data_assignments, dropped_d = convert_data_stmts(source, target_mode)
    removed_known.update(dropped_p)
    removed_known.update(dropped_d)
    source = _dedup_intrinsic_stmts(source, target_mode)
    source = insert_use_multifloats(source, target_mode, extra_lines=param_assignments + data_assignments)
    # Join continuation lines into logical statements before running the
    # paren-walking passes (replace_intrinsic_calls,
    # replace_generic_conversions, etc.). Without this, a CMPLX call
    # split across a fixed-form continuation -- ``CMPLX(re,\n+ im)`` --
    # has its open paren on one line and close paren on another, the
    # paren-balancer runs out of input on the first line and bails
    # silently, and the missing ``, KIND=k`` truncates COMPLEX(KIND=16)
    # results to default COMPLEX(KIND=4) (single precision). The non-
    # flang `migrate_fixed_form` path already does this; the flang
    # variant matches it now.
    physical = source.splitlines(keepends=True)
    statements = _segment_fixed_form_statements(physical)
    result = []
    for kind, lines, terms, joined in statements:
        if kind == 'blank':
            result.append(lines[0] + terms[0])
            continue
        if kind == 'pp':
            result.append(lines[0] + terms[0])
            continue
        if kind == 'comment':
            s = replace_routine_names(lines[0], rename_map)
            if has_float_types:
                s = replace_type_decls(s, target_mode, complex_names=complex_names,
                                        source_kind=source_kind)
            result.append(s + terms[0])
            continue
        s = joined
        if has_float_types:
            s = replace_type_decls(s, target_mode, complex_names=complex_names,
                                    source_kind=source_kind)
            if not s:
                continue
            s = replace_standalone_real_complex(s, target_mode, source_kind=source_kind)
        if has_real_literals: s = replace_literals(s, target_mode, source_kind=source_kind)
        s = replace_intrinsic_calls(
            s, target_mode, real_names=real_names, complex_names=complex_names,
        )
        s = replace_intrinsic_decls(s, target_mode)
        s = replace_generic_conversions(s, target_mode, complex_names=complex_names)
        s = replace_routine_names(s, rename_map)
        s = replace_include_filenames(s, rename_map)
        s = replace_xerbla_strings(s, rename_map)
        s = replace_known_constants(s, target_mode, renames=removed_known)
        s = _rewrite_int_of_complex(s, complex_names)
        s = _rewrite_int_kind_on_real64x2(s, target_mode, real_names=real_names)
        s = _wrap_bare_complex_literals(s, target_mode)
        s = _unwrap_redundant_constructors(s, target_mode, real_names=real_names)
        if len(lines) > 1 and s == joined:
            # Multi-line, no transforms: emit physical lines verbatim
            # to avoid reformatting churn for unchanged statements.
            for line, term in zip(lines, terms):
                result.append(line + term)
            continue
        s = reformat_fixed_line(s)
        result.append(s + terms[-1])

    source = ''.join(result)
    if not target_mode.is_kind_based:
        source = re.sub(r'! !    integer, parameter :: wp = kind\(1\.d0\)',
                        '!    integer, parameter :: wp = kind(1.d0)', source)

    source = _dedup_intrinsic_stmts(source, target_mode)
    source = specialize_use_module(source, target_mode, fixed_form=True)
    return source


def _migrate_free_form_flang(source: str, rename_map: dict[str, str], target_mode: TargetMode, has_float_types: bool, source_kind: int | None = None) -> str:
    complex_names = _scan_complex_var_names(source) if not target_mode.is_kind_based else set()
    real_names = _scan_real_var_names(source) if not target_mode.is_kind_based else set()
    source = rewrite_la_constants_use(source, target_mode)
    source = fix_misdeclared_statement_functions(source, source_kind=source_kind)
    source, removed_known = strip_known_constants_from_decls(source, target_mode)
    source, param_assignments, dropped_p = convert_parameter_stmts(source, target_mode)
    source, data_assignments, dropped_d = convert_data_stmts(source, target_mode)
    removed_known.update(dropped_p)
    removed_known.update(dropped_d)
    source = _dedup_intrinsic_stmts(source, target_mode)
    source = insert_use_multifloats(source, target_mode, extra_lines=param_assignments + data_assignments)
    lines = source.splitlines(keepends=True)
    result = []
    for line in lines:
        stripped = line.rstrip('\n\r')
        nl = '\n' if line.endswith('\n') else ''
        stripped = _replace_kind_parameter(stripped, target_mode)
        if not stripped.lstrip().startswith('!'):
            stripped = _strip_iso_fortran_env_realN(stripped)
            if not stripped: continue
            stripped = replace_intrinsic_calls(
                stripped, target_mode, real_names=real_names,
                complex_names=complex_names,
            )
            stripped = replace_intrinsic_decls(stripped, target_mode)
            stripped = replace_generic_conversions(
                stripped, target_mode, complex_names=complex_names,
            )
        stripped = replace_routine_names(stripped, rename_map)
        stripped = replace_include_filenames(stripped, rename_map)
        if stripped.lstrip().startswith('!'):
            stripped = replace_type_decls(stripped, target_mode, complex_names=complex_names,
                                           source_kind=source_kind)
        else:
            if not target_mode.is_kind_based:
                stripped = re.sub(r'REAL\s*\(\s*(?:KIND\s*=\s*)?' + _KIND_PARAM_NAMES + r'\s*\)', target_mode.real_type, stripped, flags=re.IGNORECASE)
                stripped = re.sub(r'COMPLEX\s*\(\s*(?:KIND\s*=\s*)?' + _KIND_PARAM_NAMES + r'\s*\)', target_mode.complex_type, stripped, flags=re.IGNORECASE)
            # Rewrite floating-point literals (D/E exponents and the
            # ``1.23_wp`` form). Without this, ``( 1.0_WP, 0.0_WP )``
            # complex constants in modern LAPACK files (DGEDMD/DGEDMDQ)
            # would be left untouched and gfortran would reject the
            # KIND parameter once ``wp`` itself has been stripped.
            if not target_mode.is_kind_based:
                stripped = replace_literals(stripped, target_mode, source_kind=source_kind)
            stripped = replace_known_constants(stripped, target_mode, renames=removed_known)
            stripped = _rewrite_int_of_complex(stripped, complex_names)
            stripped = _rewrite_int_kind_on_real64x2(stripped, target_mode, real_names=real_names)
            stripped = _wrap_bare_complex_literals(stripped, target_mode)
            stripped = _unwrap_redundant_constructors(stripped, target_mode, real_names=real_names)
        result.append(stripped + nl)

    source = ''.join(result)
    if not target_mode.is_kind_based:
        source = re.sub(r'(?i)!\s*!\s*integer\s*,\s*parameter\s*::\s*wp\s*=',
                        '!    integer, parameter :: wp =', source)

    source = _dedup_intrinsic_stmts(source, target_mode)
    source = specialize_use_module(source, target_mode, fixed_form=False)
    return source
