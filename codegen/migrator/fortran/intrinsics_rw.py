"""Intrinsic call / declaration / conversion rewriting (Cluster C).

Rewrites intrinsic calls and generic type conversions to their migrated
equivalents via the INTRINSIC_* maps. Extracted verbatim from
``fortran_migrator.py``.
"""
import re

from ..intrinsics import INTRINSIC_MAP, INTRINSIC_DECL_MAP
from ..target_mode import TargetMode
from .lex import is_continuation_line, split_top_level_commas


# ``REAL(`` / ``CMPLX(`` paren-call detector for ``replace_generic_conversions``.
_GENERIC_CONV_RE: dict[str, re.Pattern] = {
    _name: re.compile(rf'(?<=[=+\-*/,(.\0])\s*\b({_name})\s*\(', re.IGNORECASE)
    for _name in ('REAL', 'CMPLX')
}


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
                            args = split_top_level_commas(inner)
                            inner_first = args[0] if args else inner
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
    # Deliberately stricter than lex.is_continuation_line: comment
    # markers ('!', 'C', 'c', '*') and tab in col 6 must NOT be masked
    # with '\0' below, or the marker restore would corrupt them.
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
                    args = split_top_level_commas(inner)
                    inner_first = args[0] if args else inner
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


def _split_dcolon_sep(text: str) -> tuple[str, str]:
    """Split an optional leading ``::`` separator (with surrounding
    whitespace) off ``text``. Returns ``(sep, rest)`` where ``sep`` is
    the verbatim separator text ('' when absent)."""
    sep_m = re.match(r'\s*::\s*', text)
    if sep_m:
        return text[:sep_m.end()], text[sep_m.end():]
    return '', text


def _normalize_intrinsic_names(names_text: str,
                               target_mode: TargetMode | None) -> list[str]:
    """Comma-split an INTRINSIC name list, dedupe case-insensitively
    (order-preserving, first spelling wins) and, in wrap_constructor
    mode, drop names the target module overloads — gfortran refuses a
    USE-associated generic name when an INTRINSIC declaration in the
    same scope binds it to the standard intrinsic of incompatible
    signature."""
    names = [n.strip() for n in names_text.split(',') if n.strip()]
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
    return deduped


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
        sep, stripped = _split_dcolon_sep(stripped)
        deduped = _normalize_intrinsic_names(stripped, target_mode)
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
            if is_continuation_line(next_line):
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
            if is_continuation_line(sl):
                name_part += ' ' + sl[6:]
            else:
                stripped = sl.lstrip()
                if stripped.startswith('&'): stripped = stripped[1:]
                name_part += ' ' + stripped

        name_part = name_part.replace('&', ' ')
        # Drop an optional ``::`` separator following the keyword.
        _, name_part = _split_dcolon_sep(name_part)
        deduped = _normalize_intrinsic_names(name_part, target_mode)

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
    parts = split_top_level_commas(inner)

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


# Type-conversion intrinsics frequently drift asymmetrically between
# co-family halves: D/Z-half upstream tends to declare ``REAL`` /
# ``CMPLX`` / ``DBLE`` / ``DIMAG`` / ``DCMPLX`` / ``DCONJG`` in its
# INTRINSIC list (used to convert ``DOUBLE COMPLEX`` <-> ``DOUBLE
# PRECISION``), while S/C-half doesn't need them (kind4 default
# handles the same conversions implicitly). After kind16 migration
# both halves use the kind-promoted versions; the asymmetric
# INTRINSIC declarations remain as cosmetic text drift.
# :func:`strip_overloaded_intrinsics` strips this subset on every
# target to converge.
_TYPE_CONV_INTRINSICS = frozenset({
    'REAL', 'CMPLX', 'DCMPLX', 'DBLE',
    'AIMAG', 'CONJG', 'DCONJG', 'DIMAG',
})

# Multifloats also strips the full generic-overload set so INTRINSIC
# declarations of names that ``USE multifloats`` provides as generic
# interfaces don't clash (gfortran: "Cannot change attributes of
# USE-associated symbol").
_MF_OVERLOADED_INTRINSICS = _TYPE_CONV_INTRINSICS | {
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


def strip_overloaded_intrinsics(line: str, target_mode: TargetMode) -> str:
    """Strip target-supplied names from an INTRINSIC declaration statement.

    Runs on a single joined logical line (via ``_apply_local_passes``),
    so a multi-line INTRINSIC list is cleaned as a whole. Kind-based
    targets strip only the type-conversion subset; multifloats strips
    the full generic-overload set. Returns ``''`` when the kept list is
    empty so the statement is dropped.
    """
    strip_set = (_TYPE_CONV_INTRINSICS if target_mode.is_kind_based
                 else _MF_OVERLOADED_INTRINSICS)

    def clean_intrinsic(m):
        indent, sep, funcs_str, newline = m.group(1), m.group(2), m.group(3), m.group(4)
        kept = [f.strip() for f in funcs_str.split(',')
                if f.strip() and f.strip().upper() not in strip_set]
        return f"{indent}INTRINSIC{sep}{', '.join(kept)}{newline}" if kept else ""

    return re.sub(
        r'(?im)^([ \t]*)INTRINSIC(\s*::\s*|\s+)([A-Za-z0-9_,\s]+?)(\r?\n|$)',
        clean_intrinsic, line,
    )
