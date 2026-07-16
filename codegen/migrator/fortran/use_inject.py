"""USE-module injection / specialization (Cluster I).

Builds the ``USE ..., ONLY:`` clauses the migrated bodies need and inserts /
specializes them at the right scope. Extracted verbatim from
``fortran_migrator.py``.
"""
import re

from ..target_mode import TargetMode
from .lex import (
    _END_PROC_RE, _INTERFACE_BEGIN_RE, _INTERFACE_END_RE, _PROC_HEADER_RE,
    _count_open_parens, _looks_like_statement_function,
    _scan_local_declared_names, _scan_referenced_identifiers,
    is_continuation_line,
)
from .keepkind import _KK_DBLE_SENTINEL
from .decls import _scan_complex_var_names
from .renames import _INCLUDE_RE


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
    for h in headers:
        # Determine the end of this procedure.
        next_end = next((e for e in ends if e >= h), len(lines) - 1)
        proc_lines = lines[h:next_end + 1]
        only_clause = _build_use_only_clause(proc_lines, target_mode)
        if only_clause is None:
            # Scope references no module public name — drop the injected
            # bare ``USE <module>`` line(s) entirely. Leaving them would
            # import the operator/assignment generics into a scope that
            # can never use them, gratuitously diverging this file's
            # object from the (identical) kind10/kind16 builds.
            for k in range(h, next_end + 1):
                if use_mod_re.match(out[k]):
                    out[k] = ''
            continue
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


def _build_use_only_clause(proc_lines: list[str],
                           target_mode: TargetMode) -> str | None:
    """Compute the ``, only:`` clause for a module-based USE statement.

    Returns the empty string if the target does not use a module. Returns
    ``None`` when the scope references no module public name at all — the
    injected bare ``USE`` should then be removed entirely (its operator
    generics are dead weight). Otherwise returns
    ``", only: name1, name2, ..., operator(+), ..."`` listing the module
    public names referenced by ``proc_lines`` (minus any name that the
    procedure declares as a local variable, so that the local declaration
    is not shadowed by the use-associated generic interface).
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
    # Predict the ``real`` calls that ``_rewrite_int_of_complex`` will
    # inject when ``INT(zvar)`` / ``NINT(zvar)`` appear in a body that
    # declares a complex variable. The rewrite emits ``real(zvar)`` and
    # the only-clause must import the multifloats ``real`` generic for
    # gfortran to dispatch it. Mirror the predictive ``dble`` path:
    # add ``real`` whenever an INT/NINT call co-occurs with any complex
    # variable name in the procedure.
    if (target_mode.intrinsic_mode == 'wrap_constructor'
            and re.search(r'\b(?:int|nint)\s*\(', body_text, re.IGNORECASE)
            and _scan_complex_var_names(body_text)):
        referenced.add('real')
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
    # A scope that references no module public name (no target type, no
    # module constant, no wrapped generic like ``dble``/``real``) cannot
    # invoke any operator/assignment overload either — an overloaded
    # ``+``/``=`` only fires on a target-typed operand, and every such
    # operand is a local/dummy whose declaration names the target type,
    # which would have put that type into ``selected``. So the operator
    # generics are dead weight here. Signal removal of the bare USE line
    # with ``None`` (distinct from ``''``, "leave as-is"). If the scope
    # does need the overloads via some path the scanner missed, gfortran
    # rejects the build — a hard gate, never a silent miscompile.
    if not selected:
        return None
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
                        # in codegen/recipes/*/keep-kind.manifest so they survive
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
