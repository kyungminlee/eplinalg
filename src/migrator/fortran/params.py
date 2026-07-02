"""PARAMETER / DATA statement conversion (Cluster H).

Retypes constants declared in PARAMETER and DATA statements to the target
kind. Extracted verbatim from ``fortran_migrator.py``.
"""
import re

from ..target_mode import TargetMode
from .lex import _is_fp_value, _join_continued_lines, _scope_indices
from .decls import _scan_complex_var_names


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
