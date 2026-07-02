"""Type-declaration rewriting + declaration scanning (Cluster A).

Rewrites REAL/COMPLEX/DOUBLE type declarations to the target kind and scans
declared variable names that later clusters key off. Extracted verbatim from
``fortran_migrator.py``.
"""
import re
import functools

from ..target_mode import TargetMode
from .lex import is_continuation_line


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
