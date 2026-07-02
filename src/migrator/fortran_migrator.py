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

# Lexical + scanning primitives now live in ``fortran/lex.py`` (bottom layer
# of the split). Re-imported here so this module's namespace — and every
# caller that reaches these names through it — is unchanged.
from .fortran.lex import (  # noqa: F401  (re-export)
    is_comment_line, _iter_outside_strings, is_continuation_line,
    _find_inline_bang, _count_open_parens, _build_split_mask,
    _strip_strings_and_comments, _looks_like_statement_function,
    _join_continued_lines, _is_fp_value, _scope_indices,
    _scan_local_declared_names, _scan_referenced_identifiers,
    _STRING_RE, _IDENT_RE, _DECL_LINE_RE, _STMT_FN_RE, _FP_VALUE_RE,
    _FORTRAN_OP_RE, _INT_CALL_RE, _NINT_CALL_RE,
    _PROC_HEADER_RE, _END_PROC_RE, _PROC_HEADER_RE_SCOPE,
    _INTERFACE_BEGIN_RE, _INTERFACE_END_RE,
)
# --- re-exported rewriter submodules (facade); see fortran/ subpackage ---
# REEXPORT-ANCHOR
from .fortran.params import (  # noqa: F401  (re-export)
    convert_data_stmts, convert_parameter_stmts,
)
from .fortran.literals import (  # noqa: F401  (re-export)
    replace_literals, _rewrite_int_kind_on_real64x2, _rewrite_int_of_complex,
)
from .fortran.decls import (  # noqa: F401  (re-export)
    _TD_DECL_TAIL, _TD_DBL_PREC, _TD_DBL_CMPLX, _TD_CMPLX_STAR16, _TD_REAL_STAR8, _TD_CMPLX_STAR8, _TD_REAL_STAR4, _TD_REAL_KIND_E0, _TD_CMPLX_KIND_E0, _TD_REAL_KIND_D0, _TD_CMPLX_KIND_D0, _TD_REAL_KIND_WP, _TD_CMPLX_KIND_WP, _TD_GATE_RE, _TD_REAL_KIND_4, _TD_CMPLX_KIND_4, _TD_REAL_KIND_8, _TD_CMPLX_KIND_8, replace_type_decls, _filter_known_decl_re, _filter_known_constants_from_decl, replace_standalone_real_complex, _DECL_START_RE, _scan_typed_var_names, _COMPLEX_DECL_RE, _REAL_DECL_RE, _scan_complex_var_names, _scan_real_var_names, strip_known_constants_from_decls, fix_misdeclared_statement_functions,
)
from .fortran.fixedform import (  # noqa: F401  (re-export)
    reformat_fixed_line, _strip_inline_comment, _segment_fixed_form_statements,
)
from .fortran.renames import (  # noqa: F401  (re-export)
    _ROUTINE_NAME_TOK_RE, _get_rename_pattern, replace_routine_names, _INCLUDE_RE, _RENAME_PATTERN_CACHE, replace_include_filenames, _XERBLA_STR_RE, replace_xerbla_strings, _known_constants_pattern, _all_known_constant_renames, _DECL_KEYWORD_RE, replace_known_constants,
)
from .fortran.intrinsics_rw import (  # noqa: F401  (re-export)
    _INTRINSIC_CALL_GATE_RE, replace_intrinsic_calls, replace_generic_conversions, _wrap_complex_args, replace_intrinsic_decls, _dedup_intrinsic_stmts, _GENERIC_CONV_RE, _INTRINSIC_CALL_RE, _INTRINSIC_CALL_RE_REPL,
)
from .fortran.complex_rw import (  # noqa: F401  (re-export)
    _wrap_bare_complex_literals, _unwrap_redundant_constructors, _unwrap_ctor_re,
)
from .fortran.mpi import (  # noqa: F401  (re-export)
    _MPI_DOUBLE_COMPLEX_RE, _MPI_DOUBLE_PRECISION_RE, _MPI_COMPLEX_RE, _MPI_REAL_RE, _INCLUDE_PREFIXED_H_RE, _MPI_REDUCE_CALL_RE, _MPI_OP_RE, _STANDARD_MPI_HANDLES, _rewrite_prefixed_includes, _rewrite_mpi_datatypes, _rewrite_mpi_sum, _custom_mpi_tokens, insert_use_multifloats_mpi_f,
)
from .fortran.keepkind import (  # noqa: F401  (re-export)
    _KK_SENTINEL, _KK_DBLE_SENTINEL, _KK_DCMPLX_SENTINEL, _apply_keep_kind_sentinel, _restore_keep_kind_sentinel, _strip_roundup_lwork,
)


# ---------------------------------------------------------------------------
# Type declaration replacement
# ---------------------------------------------------------------------------





















# ---------------------------------------------------------------------------
# Literal constant replacement
# ---------------------------------------------------------------------------



# ---------------------------------------------------------------------------
# Intrinsic function replacement
# ---------------------------------------------------------------------------


































# Note: an experimental ``_force_int_assignment`` pass was prototyped
# here that wraps the RHS of ``INT_VAR = ...`` with ``INT(...)`` when
# the RHS mentions a known float64x2 variable. It was removed because
# the heuristic ("any token in real_names") misclassifies the case
# where a float64x2 variable is *passed* to an integer-returning
# function (e.g. ``JP = J - 1 + IUAMAX(M-J+1, A(J,J), 1)`` where A is
# float64x2 but IUAMAX returns INTEGER). Reliable handling needs
# semantic facts (the migrated function's return type), which is
# Phase 1.5 work.




































# Module public names (type names, constants, generics, operator generics)
# are now loaded from the target YAML via TargetMode fields:
#   target_mode.module_type_names
#   target_mode.module_constant_names
#   target_mode.module_generic_names
#   target_mode.module_public_names  (union of the above three)
#   target_mode.module_operator_generics









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
