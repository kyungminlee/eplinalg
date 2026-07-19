"""General-purpose Fortran type migration engine.

Handles fixed-form (.f) and free-form (.f90) Fortran source files.
Works for any library following the S/D/C/Z precision prefix convention
(BLAS, LAPACK, ScaLAPACK, MUMPS, etc.).

Uses a hybrid strategy:
  1. Parse with a Fortran compiler (Flang or GFortran, selected via
     ``--parser``) to get a structured parse tree
  2. Extract facts: type declarations, routine names, call sites,
     literals, intrinsics, XERBLA strings
  3. Apply transformations as source-level replacements, preserving
     all original formatting, comments, and preprocessor directives

Runs regex-only scanning when no ``--parser`` is given.

Transformations:
  1. Type declarations (DOUBLE PRECISION → REAL(KIND=k), etc.)
  2. Standalone REAL/COMPLEX in declaration context
  3. Floating-point literals (D-exponent → E_kind suffix)
  4. Type-specific intrinsic function calls
  5. INTRINSIC declaration statements
  6. Routine name prefixes (using rename map)
  7. XERBLA string arguments

The per-concern rewrite passes live in the ``fortran/`` subpackage; this
module hosts only the top-level drivers that sequence them and imports
just the names those drivers use.
"""

import dataclasses
import re
from pathlib import Path

from .target_mode import TargetMode

from .fortran.use_inject import insert_use_multifloats, specialize_use_module
from .fortran.la_constants import (
    _KIND_PARAM_NAMES, _replace_kind_parameter, _strip_iso_fortran_env_realN,
    rewrite_la_constants_use,
)
from .fortran.params import (
    convert_data_stmts, convert_parameter_stmts, nuke_multifloats_params,
)
from .fortran.literals import (
    replace_literals, _rewrite_int_kind_on_real64x2, _rewrite_int_of_complex,
)
from .fortran.decls import (
    replace_type_decls, replace_standalone_real_complex,
    _scan_complex_var_names, _scan_real_var_names,
    strip_known_constants_from_decls, fix_misdeclared_statement_functions,
)
from .fortran.fixedform import reformat_fixed_line, _segment_fixed_form_statements
from .fortran.logical_lines import segment_free_form, reformat_free_line
from .fortran.renames import (
    replace_routine_names, replace_include_filenames, replace_xerbla_strings,
    replace_known_constants,
)
from .fortran.intrinsics_rw import (
    replace_intrinsic_calls, replace_generic_conversions,
    replace_intrinsic_decls, _dedup_intrinsic_stmts,
    strip_overloaded_intrinsics,
)
from .fortran.complex_rw import (
    _wrap_bare_complex_literals, _unwrap_redundant_constructors,
)
from .fortran.mpi import (
    _rewrite_prefixed_includes, _rewrite_mpi_datatypes, _rewrite_mpi_sum,
    insert_use_multifloats_mpi_f,
)
from .fortran.keepkind import (
    _apply_keep_kind_sentinel, _restore_keep_kind_sentinel,
    _strip_roundup_lwork,
)
from .fortran.io_narrow import (
    narrow_multifloats_io,
    narrow_multifloats_io_open,
    narrow_io_continuation,
    is_fixed_io_continuation,
)


# ---------------------------------------------------------------------------
# Orchestrators
#
# Two design notes that outlived their original call sites:
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


@dataclasses.dataclass(frozen=True)
class MigrationContext:
    """Per-file-invariant inputs shared by the migration drivers.

    ``rename_map`` maps upper-cased source routine/module names to their
    migrated names; ``source_kind`` is the source half's native kind (4
    or 8, see :func:`_source_kind_from_filename`) or None.

    ``comp_real`` / ``comp_complex`` are the global derived-type component
    oracle (real64x2 / cmplx64x2 struct fields harvested across the whole
    staged tree); formatted-output narrowing consults them for
    ``%``-qualified references (e.g. ``id%DKEEP(160)``) whose type is
    declared in a struct module this file only ``USE``s.
    """
    rename_map: dict[str, str]
    target_mode: TargetMode
    source_kind: int | None = None
    comp_real: frozenset[str] = frozenset()
    comp_complex: frozenset[str] = frozenset()


def _apply_intrinsic_passes(s: str, target_mode: TargetMode,
                            real_names: set[str],
                            complex_names: set[str]) -> str:
    """Intrinsic-rewrite chain shared by the fixed- and free-form drivers."""
    s = replace_intrinsic_calls(
        s, target_mode, real_names=real_names, complex_names=complex_names,
    )
    s = replace_intrinsic_decls(s, target_mode)
    s = replace_generic_conversions(s, target_mode, complex_names=complex_names)
    return s


def _apply_value_passes(s: str, target_mode: TargetMode,
                        removed_known: dict[str, str],
                        real_names: set[str],
                        complex_names: set[str]) -> str:
    """Constant/value-rewrite chain shared by the fixed- and free-form
    drivers (known-constant renames, INT-of-complex, complex literal
    wrapping/unwrapping)."""
    s = replace_known_constants(s, target_mode, renames=removed_known)
    s = _rewrite_int_of_complex(s, complex_names)
    s = _rewrite_int_kind_on_real64x2(s, target_mode, real_names=real_names)
    s = _wrap_bare_complex_literals(s, target_mode)
    s = _unwrap_redundant_constructors(s, target_mode, real_names=real_names)
    return s


def migrate_fixed_form(source: str, ctx: MigrationContext, *,
                        has_float_types: bool = True,
                        has_real_literals: bool = True) -> str:
    """Migrate one fixed-form source string.

    ``has_float_types`` / ``has_real_literals`` let the parser-guided
    caller (:func:`_migrate_with_facts`) skip the type-declaration /
    literal-promotion passes when the parse facts show the file contains
    no floating-point declarations / real literals. The parser-less path
    has no such evidence and runs every pass (the defaults).

    Note the fixed-form-only passes with no free-form counterpart:
    :func:`replace_standalone_real_complex` and
    :func:`replace_xerbla_strings` — a known, intentional divergence
    from :func:`migrate_free_form`.
    """
    rename_map, target_mode = ctx.rename_map, ctx.target_mode
    source_kind = ctx.source_kind
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
    # Tracks a formatted WRITE/PRINT whose output list spills across
    # cpp-broken continuation fragments (see io_narrow docstring): the
    # segmenter splits the statement at each ``#if``/``#endif``, so the
    # fragment carrying the real64x2 item has no WRITE keyword and must be
    # narrowed by continuation. Reset the moment a non-continuation
    # statement begins.
    io_list_open = False
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
        # 'code' — apply all per-line transforms to the joined logical
        # line so paren-walking passes (e.g. replace_intrinsic_calls)
        # can match calls that span fixed-form continuations. Single-
        # physical-line statements have ``joined == lines[0]`` and pass
        # through identically.
        s = joined
        if has_float_types:
            s = replace_type_decls(s, target_mode, complex_names=complex_names,
                                    source_kind=source_kind)
            if not s:
                continue
            s = replace_standalone_real_complex(s, target_mode, source_kind=source_kind)
        if has_real_literals:
            s = replace_literals(s, target_mode, source_kind=source_kind)
        s = _apply_intrinsic_passes(s, target_mode, real_names, complex_names)
        s = replace_routine_names(s, rename_map)
        s = replace_include_filenames(s, rename_map)
        s = replace_xerbla_strings(s, rename_map)
        s = _apply_value_passes(s, target_mode, removed_known,
                                real_names, complex_names)
        if is_fixed_io_continuation(lines[0]):
            # Headless continuation fragment of a cpp-broken statement.
            # Narrow it only while a formatted output list is still open;
            # leave its open/closed state to the statement that owns it.
            if io_list_open:
                s = narrow_io_continuation(s, real_names, complex_names,
                                           ctx.comp_real, ctx.comp_complex)
        else:
            # A new logical statement begins: close any open output list,
            # then re-open if this statement is itself a formatted WRITE.
            io_list_open = False
            s, io_list_open = narrow_multifloats_io_open(
                s, real_names, complex_names, ctx.comp_real, ctx.comp_complex)
        if len(lines) > 1 and s == joined:
            # Multi-line statement, no transforms applied — emit the
            # original physical lines verbatim to avoid reformat churn.
            for line, term in zip(lines, terms):
                result.append(line + term)
            continue
        s = reformat_fixed_line(s)
        # A reflowed multi-line statement replaces all of its physical
        # lines, so the surviving terminator is the LAST line's — this
        # preserves a missing newline at end-of-file.
        result.append(s + terms[-1])

    source = ''.join(result)
    if not target_mode.is_kind_based:
        source = re.sub(r'! !    integer, parameter :: wp = kind\(1\.d0\)',
                        '!    integer, parameter :: wp = kind(1.d0)', source)

    source = _dedup_intrinsic_stmts(source, target_mode)
    source = specialize_use_module(source, target_mode, fixed_form=True)
    return source


def migrate_free_form(source: str, ctx: MigrationContext, *,
                       nuke_mf_params: bool = True) -> str:
    """Migrate one free-form source string.

    ``nuke_mf_params`` controls the heuristic pass that comments out
    file-scope PARAMETER declarations whose names the multifloats module
    supplies as ``DD_*`` constants. The parser-guided caller
    (:func:`_migrate_with_facts`) has never applied that pass and
    disables it; the parser-less path keeps it (the default).
    """
    rename_map, target_mode = ctx.rename_map, ctx.target_mode
    source_kind = ctx.source_kind
    complex_names = _scan_complex_var_names(source) if not target_mode.is_kind_based else set()
    real_names = _scan_real_var_names(source) if not target_mode.is_kind_based else set()
    source = rewrite_la_constants_use(source, target_mode)
    source = fix_misdeclared_statement_functions(source, source_kind=source_kind)
    source, removed_known = strip_known_constants_from_decls(source, target_mode)
    if nuke_mf_params and not target_mode.is_kind_based:
        source = nuke_multifloats_params(source, removed_known)

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
            stripped = _apply_intrinsic_passes(stripped, target_mode,
                                               real_names, complex_names)
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
            stripped = _apply_value_passes(stripped, target_mode, removed_known,
                                           real_names, complex_names)
            stripped = narrow_multifloats_io(stripped, real_names, complex_names,
                                             ctx.comp_real, ctx.comp_complex)
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


def migrate_file_to_string(src_path: Path, rename_map: dict[str, str],
                           target_mode: TargetMode,
                           parser: str | None = None,
                           parser_cmd: str | None = None,
                           keep_kind_lines: frozenset[int] | None = None,
                           comp_real: frozenset[str] = frozenset(),
                           comp_complex: frozenset[str] = frozenset(),
                           ) -> tuple[str, str] | None:
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

    ctx = MigrationContext(
        rename_map=rename_map,
        target_mode=target_mode,
        source_kind=_source_kind_from_filename(src_path.name),
        comp_real=comp_real,
        comp_complex=comp_complex,
    )

    if facts is not None:
        migrated = _migrate_with_facts(source, ext, ctx, facts)
    elif ext in ('.f', '.for', '.h'):
        migrated = migrate_fixed_form(source, ctx)
    elif ext in ('.f90', '.f95', '.F90'):
        migrated = migrate_free_form(source, ctx)
    else:
        return None

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

    # Strip target-supplied names from INTRINSIC declarations (see
    # strip_overloaded_intrinsics for the co-family-drift / generic-
    # interface-clash rationale). Runs in the verbatim-preserving joined
    # regime: multi-line INTRINSIC lists are collapsed before the strip,
    # statements it doesn't alter keep their original physical
    # formatting, and re-wrapping happens once here.
    migrated = _apply_local_passes(
        migrated,
        [lambda s: strip_overloaded_intrinsics(s, target_mode)],
        fixed_form=fixed_form,
    )
    return out_name, migrated


def _migrate_with_facts(source: str, ext: str, ctx: MigrationContext,
                        facts) -> str:
    # Also include USE-statement module names + every variable's type
    # spec (where the type is precision-prefixed, e.g.
    # ``TYPE(DMUMPS_INTR_STRUC)``). The gfortran parser doesn't surface
    # type-of-variable as a call/external_name, so without this extra
    # source the rename_map gets filtered down to nothing for files
    # whose only precision-prefixed reference is via ``USE`` or
    # ``TYPE(...)``.
    use_names = set(facts.use_modules or ())
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
        | set(facts.call_sites)
        | set(facts.external_names)
        | use_names
        | type_refs
    )
    file_ctx = dataclasses.replace(
        ctx,
        rename_map={k: v for k, v in ctx.rename_map.items() if k in file_names},
    )
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
    if ext in ('.f90', '.f95', '.F90'):
        return migrate_free_form(source, file_ctx, nuke_mf_params=False)
    return migrate_fixed_form(source, file_ctx,
                              has_float_types=has_float_types,
                              has_real_literals=has_real_literals)
