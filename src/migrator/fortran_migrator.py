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
from .fortran.use_inject import (  # noqa: F401  (re-export)
    _wrap_use_clause, insert_use_multifloats, _build_use_only_clause, specialize_use_module,
)
from .fortran.la_constants import (  # noqa: F401  (re-export)
    _KIND_PARAM_NAMES, _KIND_PARAM_RE, _replace_kind_parameter, _ISO_USE_ONLY_RE, _strip_iso_fortran_env_realN, _LA_CONSTANTS_REAL_NAMES, _LA_CONSTANTS_COMPLEX_NAMES, _la_constants_rename_map, rewrite_la_constants_use,
)
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
from .fortran.logical_lines import (  # noqa: F401  (re-export)
    segment_free_form, reformat_free_line,
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
    _MPI_DOUBLE_COMPLEX_RE, _MPI_DOUBLE_PRECISION_RE, _MPI_COMPLEX_RE, _MPI_REAL_RE, _INCLUDE_PREFIXED_H_RE, _MPI_REDUCE_HEAD_RE, _sub_mpi_reduce_calls, _MPI_OP_RE, _STANDARD_MPI_HANDLES, _rewrite_prefixed_includes, _rewrite_mpi_datatypes, _rewrite_mpi_sum, _custom_mpi_tokens, insert_use_multifloats_mpi_f,
)
from .fortran.keepkind import (  # noqa: F401  (re-export)
    _KK_SENTINEL, _KK_DBLE_SENTINEL, _KK_DCMPLX_SENTINEL, _apply_keep_kind_sentinel, _restore_keep_kind_sentinel, _strip_roundup_lwork,
)


# ---------------------------------------------------------------------------
# Orchestrators
#
# The per-concern rewrite passes now live in the ``fortran/`` subpackage and
# are re-exported above; what remains here are the top-level drivers that
# sequence them. Two design notes that outlived their original call sites:
#
# * An experimental ``_force_int_assignment`` pass (wrap the RHS of
#   ``INT_VAR = ...`` with ``INT(...)`` when the RHS mentions a known
#   float64x2 variable) was prototyped and removed: the heuristic ("any
#   token in real_names") misclassifies a float64x2 variable *passed* to an
#   integer-returning function (e.g. ``JP = J - 1 + IUAMAX(M-J+1, A(J,J), 1)``
#   where A is float64x2 but IUAMAX returns INTEGER). Reliable handling needs
#   semantic facts (the migrated return type) — Phase 1.5 work.
# * Module public names (type names, constants, generics, operator generics)
#   are loaded from the target YAML via TargetMode fields: module_type_names,
#   module_constant_names, module_generic_names, module_public_names (their
#   union), and module_operator_generics.
# ---------------------------------------------------------------------------


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


def _apply_local_passes(source: str, passes, *, fixed_form: bool) -> str:
    """Run ``passes`` over each statement's joined logical line so a pattern
    the source split across continuations is matched as one string, but
    re-emit statements the passes leave unchanged as their *original* physical
    lines. This gives the join regime's matching power without the diff churn
    of reflowing every multi-line statement in the file — a whole-file
    join -> re-wrap round-trip reformats ~386/436 MUMPS files cosmetically,
    this touches only the statements a pass rewrote. A pass returning ``''``
    drops the statement (e.g. an emptied INTRINSIC)."""
    physical = source.splitlines(keepends=True)
    segs = (_segment_fixed_form_statements(physical) if fixed_form
            else segment_free_form(physical))
    reflow = reformat_fixed_line if fixed_form else reformat_free_line
    out: list[str] = []
    for kind, lines, terms, joined in segs:
        if kind != 'code':
            for line, term in zip(lines, terms):
                out.append(line + term)
            continue
        new = joined
        for p in passes:
            new = p(new)
        if new == joined:
            for line, term in zip(lines, terms):
                out.append(line + term)
        elif new == '':
            continue
        else:
            term = terms[0] if terms else '\n'
            out.append(reflow(new) + term)
    return ''.join(out)


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

    # The MPI datatype/reduce-op rewrites match patterns that can span
    # continuation lines (an MPI_REDUCE call, its op and datatype args), so run
    # them through the continuation-collapsed, verbatim-preserving regime:
    # each statement is matched as one joined logical line, but statements the
    # passes don't touch keep their original physical formatting. Previously
    # each pass re-walked continuations itself and periodically missed cases
    # (e.g. an MPI_REDUCE op whose argument nested two paren levels).
    fixed_form = ext in ('.f', '.for', '.h')
    migrated = _apply_local_passes(
        migrated,
        [lambda s: _rewrite_mpi_datatypes(s, target_mode),
         lambda s: _rewrite_mpi_sum(s, target_mode)],
        fixed_form=fixed_form,
    )
    # The remaining post-passes insert a whole USE line / rewrite single
    # (never-continued) include lines / strip a LWORK round-up — none match
    # across continuations, so they stay as whole-file passes. insert_use must
    # run after _rewrite_mpi_sum (it keys off the MPI_QQ_* ops that pass emits).
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
        kept = [f.strip() for f in funcs_str.split(',')
                if f.strip() and f.strip().upper() not in strip_set]
        return f"{indent}INTRINSIC{sep}{', '.join(kept)}{newline}" if kept else ""

    def clean_intrinsic_line(s: str) -> str:
        # Runs on a single joined logical line, so a multi-line INTRINSIC list
        # (previously punted on to avoid fusing ``MAXEXPONENT,\n$MINEXPONENT``
        # into one token) is now cleaned as a whole. Returns ``''`` when the
        # kept list is empty so the statement is dropped.
        return re.sub(
            r'(?im)^([ \t]*)INTRINSIC(\s*::\s*|\s+)([A-Za-z0-9_,\s]+?)(\r?\n|$)',
            clean_intrinsic, s,
        )

    # Runs in the verbatim-preserving joined regime: multi-line INTRINSIC
    # lists are collapsed before the strip, statements it doesn't alter keep
    # their original physical formatting, and re-wrapping happens once here.
    migrated = _apply_local_passes(migrated, [clean_intrinsic_line], fixed_form=fixed_form)
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
