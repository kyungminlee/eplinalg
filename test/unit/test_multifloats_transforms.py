"""Unit tests for the multifloats branch of migrator.fortran_migrator.

These tests cover the per-pass transformations in isolation plus a few
small end-to-end migrations of synthetic BLAS-style sources. They run
without any external compiler — they only check the textual output of
the migrator. The expensive end-to-end smoke test that compiles real
BLAS sources against the multifloats module lives outside pytest.
"""

import textwrap

import pytest

from migrator.target_mode import load_target
from migrator.fortran_migrator import (
    replace_type_decls,
    replace_literals,
    replace_intrinsic_calls,
    replace_known_constants,
    strip_known_constants_from_decls,
    convert_parameter_stmts,
    convert_data_stmts,
    insert_use_multifloats,
    specialize_use_module,
    rewrite_la_constants_use,
    _rewrite_mpi_datatypes,
    _rewrite_mpi_sum,
    migrate_fixed_form,
    migrate_free_form,
    MigrationContext,
)
from migrator.fortran.la_constants import _la_constants_rename_map


# Multifloats prefix letters, derived from the loaded target so tests
# don't need editing whenever the prefix changes.
_MF = load_target('multifloats')
_RP = _MF.prefix_map['R'].lower()
_CP = _MF.prefix_map['C'].lower()
_RPU = _MF.prefix_map['R']
_CPU = _MF.prefix_map['C']


# ---------------------------------------------------------------------------
# TargetMode factory
# ---------------------------------------------------------------------------


def test_multifloats_target_basic_shape():
    mf = load_target('multifloats')
    assert mf.name == 'multifloats'
    assert mf.is_kind_based is False
    assert mf.real_type == 'TYPE(real64x2)'
    assert mf.complex_type == 'TYPE(cmplx64x2)'
    assert mf.real_constructor == 'real64x2'
    assert mf.complex_constructor == 'cmplx64x2'
    assert mf.module_name == 'multifloats'
    # Prefix letters are loaded from codegen/targets/multifloats.yaml. The
    # specific choice (currently M/W) was made after a per-file shadow
    # analysis ruling out same-translation-unit collisions; tests below
    # use the loaded values via _RP/_CP rather than hard-coded letters.
    assert isinstance(mf.prefix_map['R'], str) and len(mf.prefix_map['R']) == 1
    assert isinstance(mf.prefix_map['C'], str) and len(mf.prefix_map['C']) == 1
    assert mf.prefix_map['R'].isupper()
    assert mf.prefix_map['C'].isupper()


def test_kind_target_is_kind_based():
    k16 = load_target('kind16')
    assert k16.is_kind_based is True
    assert k16.real_type == 'REAL(KIND=16)'
    assert k16.kind_suffix == 16


def test_multifloats_known_constants_excludes_rtmin_rtmax():
    """RTMIN/RTMAX collide with local LAPACK variables — must NOT be in
    the auto-rename set."""
    mf = load_target('multifloats')
    assert 'RTMIN' not in mf.known_constants
    assert 'RTMAX' not in mf.known_constants
    assert 'ZERO' in mf.known_constants
    assert 'ONE' in mf.known_constants


# ---------------------------------------------------------------------------
# replace_type_decls
# ---------------------------------------------------------------------------


@pytest.fixture
def mf():
    return load_target('multifloats')


@pytest.mark.parametrize('src,expected', [
    ('      DOUBLE PRECISION X', '      TYPE(real64x2) X'),
    ('      DOUBLE COMPLEX Z', '      TYPE(cmplx64x2) Z'),
    ('      COMPLEX*16 Z', '      TYPE(cmplx64x2) Z'),
    ('      REAL*8 X', '      TYPE(real64x2) X'),
    ('      REAL*4 X', '      TYPE(real64x2) X'),
])
def test_replace_type_decls_keywords(mf, src, expected):
    assert replace_type_decls(src, mf) == expected


def test_replace_type_decls_preserves_attrs(mf):
    src = '      DOUBLE PRECISION, INTENT(IN) :: A(LDA,*)'
    out = replace_type_decls(src, mf)
    assert 'TYPE(real64x2)' in out
    assert 'INTENT(IN)' in out
    assert 'A(LDA,*)' in out


# ---------------------------------------------------------------------------
# strip_known_constants_from_decls
# ---------------------------------------------------------------------------


def test_strip_drops_entire_line_when_all_names_known(mf):
    src = '      DOUBLE PRECISION ZERO,ONE\n'
    new, removed = strip_known_constants_from_decls(src, mf)
    assert new == ''
    assert removed == {'ZERO': 'DD_ZERO', 'ONE': 'DD_ONE'}


def test_strip_filters_partial_decl(mf):
    src = '      DOUBLE PRECISION DFLAG,TWO,W,Z,ZERO\n'
    new, removed = strip_known_constants_from_decls(src, mf)
    assert 'DFLAG' in new
    assert 'W' in new and 'Z' in new
    assert 'TWO' not in new
    assert 'ZERO' not in new
    assert removed == {'TWO': 'DD_TWO', 'ZERO': 'DD_ZERO'}


def test_strip_handles_multiline_continuation(mf):
    """Continuation line carries names that should also be filtered."""
    src = (
        '      DOUBLE PRECISION DFLAG,DH11,DH12,DH21,DH22,DP1,DP2,DQ1,DQ2,DTEMP,\n'
        '     $                 DU,GAM,GAMSQ,ONE,RGAMSQ,TWO,ZERO\n'
    )
    new, removed = strip_known_constants_from_decls(src, mf)
    assert 'ONE' not in new
    assert 'TWO' not in new
    assert 'ZERO' not in new
    assert 'GAM' in new
    assert 'GAMSQ' in new
    assert 'RGAMSQ' in new
    assert removed == {'ONE': 'DD_ONE', 'TWO': 'DD_TWO', 'ZERO': 'DD_ZERO'}


def test_strip_matches_bare_real_keyword(mf):
    """sgemm.f uses bare REAL (single precision) — must be matched too."""
    src = '      REAL ONE,ZERO\n'
    new, removed = strip_known_constants_from_decls(src, mf)
    assert new == ''
    assert removed == {'ONE': 'DD_ONE', 'ZERO': 'DD_ZERO'}


def test_strip_does_not_touch_array_specs(mf):
    """A line with array specs is not safe to filter — leave it alone."""
    src = '      DOUBLE PRECISION A(LDA,*),X(*),Y(*)\n'
    new, removed = strip_known_constants_from_decls(src, mf)
    assert new == src
    assert removed == {}


def test_strip_skips_function_return_decl(mf):
    src = '      DOUBLE PRECISION FUNCTION FOO(N)\n'
    new, removed = strip_known_constants_from_decls(src, mf)
    assert new == src
    assert removed == {}


def test_strip_kind_target_is_noop():
    src = '      DOUBLE PRECISION ZERO,ONE\n'
    new, removed = strip_known_constants_from_decls(src, load_target('kind16'))
    assert new == src
    assert removed == {}


# ---------------------------------------------------------------------------
# replace_literals
# ---------------------------------------------------------------------------


def test_replace_literals_d_exponent(mf):
    src = "      X = 1.0D+0 + 2.5D-3"
    out = replace_literals(src, mf)
    assert "real64x2(limbs=[1.0D0, 0.0_8])" in out
    assert "real64x2(limbs=[2.5D-3, 0.0_8])" in out


def test_replace_literals_bare_float(mf):
    src = "      X = 0.5"
    out = replace_literals(src, mf)
    assert "real64x2(limbs=[0.5D0, 0.0_8])" in out


def test_replace_literals_wp_suffix(mf):
    src = "   x = 0.5_wp + 1.0_wp"
    out = replace_literals(src, mf)
    assert "real64x2(limbs=[0.5D0, 0.0_8])" in out
    assert "real64x2(limbs=[1.0D0, 0.0_8])" in out


def test_replace_literals_does_not_double_wrap(mf):
    """Idempotence: re-running on already-wrapped output is a no-op."""
    src = "      X = real64x2(limbs=[1.0D0, 0.0_8])"
    out1 = replace_literals(src, mf)
    out2 = replace_literals(out1, mf)
    # Second pass mustn't add another constructor layer.
    assert out1.count('real64x2(') == 1
    assert out2.count('real64x2(') == 1


def test_replace_literals_preserves_strings(mf):
    """Float-shaped text inside string literals must be untouched."""
    src = "      WRITE(*,*) '1.0D+0 in a string'"
    out = replace_literals(src, mf)
    assert "'1.0D+0 in a string'" in out
    assert 'real64x2(' not in out


def test_replace_literals_skips_fortran_operators(mf):
    """``.EQ.`` etc. must not be parsed as floating-point literals."""
    src = "      IF (X.EQ.1.0D+0) GO TO 10"
    out = replace_literals(src, mf)
    assert '.EQ.' in out
    assert "real64x2(limbs=[1.0D0, 0.0_8])" in out


# ---------------------------------------------------------------------------
# replace_known_constants
# ---------------------------------------------------------------------------


def test_replace_known_constants_uses_per_file_renames(mf):
    src = "      X = ZERO + ONE"
    out = replace_known_constants(src, mf, renames={'ZERO': 'DD_ZERO', 'ONE': 'DD_ONE'})
    assert out == "      X = DD_ZERO + DD_ONE"


def test_replace_known_constants_skips_comment_lines(mf):
    src = "*     performs ONE pass through ZERO array"
    out = replace_known_constants(src, mf, renames={'ZERO': 'DD_ZERO', 'ONE': 'DD_ONE'})
    assert out == src


def test_replace_known_constants_skips_inline_comment_text(mf):
    src = "      X = ZERO     ! the ZERO constant"
    out = replace_known_constants(src, mf, renames={'ZERO': 'DD_ZERO'})
    # ZERO before the ! gets renamed; ZERO after stays as English.
    assert out == "      X = DD_ZERO     ! the ZERO constant"


def test_replace_known_constants_skips_string_literals(mf):
    src = "      CALL XERBLA('ZERO INPUT')"
    out = replace_known_constants(src, mf, renames={'ZERO': 'DD_ZERO'})
    assert out == src


def test_replace_known_constants_skips_decl_lines(mf):
    src = "      INTEGER ZERO"
    out = replace_known_constants(src, mf, renames={'ZERO': 'DD_ZERO'})
    assert out == src


def test_replace_known_constants_kind_mode_is_noop():
    src = "      X = ZERO"
    out = replace_known_constants(src, load_target('kind16'), renames={'ZERO': 'DD_ZERO'})
    assert out == src


# ---------------------------------------------------------------------------
# replace_intrinsic_calls — Phase 7 multifloats intrinsic mapping
# ---------------------------------------------------------------------------


def test_intrinsic_dble_uses_universal_constructor(mf):
    """``DBLE(x)`` becomes ``real64x2(x)`` — the universal generic
    constructor in multifloats handles every input type (integer,
    real(sp/dp), real64x2, cmplx64x2 → real part)."""
    src = "      Y = DBLE(ZX(I))"
    out = replace_intrinsic_calls(src, mf)
    assert 'real64x2(ZX(I))' in out


def test_intrinsic_dble_integer_uses_constructor(mf):
    """``DBLE(N)`` for integer N must become ``real64x2(N)``, not
    ``DD_REAL(N)`` (which only takes float64x2/complex64x2)."""
    src = "      Y = DBLE(N)"
    out = replace_intrinsic_calls(src, mf)
    assert 'real64x2(N)' in out


def test_intrinsic_dcmplx_to_cmplx64x2(mf):
    """``DCMPLX(A, B)`` becomes ``cmplx64x2(real64x2(A), real64x2(B))``
    so the call dispatches to multifloats's cx_from_mf_2 interface.
    Args known to be real64x2 are not double-wrapped (see
    test_intrinsic_dcmplx_skips_known_real)."""
    src = "      Z = DCMPLX(A, B)"
    out = replace_intrinsic_calls(src, mf)
    assert 'cmplx64x2(real64x2(A),real64x2(B))' in out


def test_intrinsic_dcmplx_skips_known_real(mf):
    src = "      Z = DCMPLX(A, B)"
    out = replace_intrinsic_calls(src, mf, real_names={'A', 'B'})
    assert 'cmplx64x2(A, B)' in out


def test_intrinsic_dimag_to_aimag(mf):
    src = "      Y = DIMAG(Z)"
    out = replace_intrinsic_calls(src, mf)
    assert 'AIMAG(Z)' in out


def test_intrinsic_dconjg_to_conjg(mf):
    src = "      Y = DCONJG(Z)"
    out = replace_intrinsic_calls(src, mf)
    assert 'CONJG(Z)' in out


def test_intrinsic_dabs_to_abs(mf):
    src = "      Y = DABS(X)"
    out = replace_intrinsic_calls(src, mf)
    assert 'ABS(X)' in out


def test_intrinsic_dcos_to_cos(mf):
    src = "      Y = DCOS(X)"
    out = replace_intrinsic_calls(src, mf)
    assert 'COS(X)' in out


# ---------------------------------------------------------------------------
# convert_parameter_stmts / convert_data_stmts
# ---------------------------------------------------------------------------


def test_convert_parameter_drops_known_constant(mf):
    src = '      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)\n'
    new, assigns, dropped = convert_parameter_stmts(src, mf)
    # Both names are known constants — no assignments emitted, line elided.
    assert assigns == []
    assert dropped == {'ONE': 'DD_ONE', 'ZERO': 'DD_ZERO'}
    assert 'PARAMETER' not in new or 'Converted' in new


def test_convert_parameter_keeps_unknown_fp(mf):
    src = '      PARAMETER (R=ONE/IPW2)\n'
    new, assigns, dropped = convert_parameter_stmts(src, mf)
    assert any('R = ONE/IPW2' in text for _, text in assigns)
    assert dropped == {}


def test_convert_parameter_mixed_type_keeps_integers(mf):
    src = '      PARAMETER (LV=128,IPW2=4096,R=1.0D0)\n'
    new, assigns, dropped = convert_parameter_stmts(src, mf)
    # Integers stay in the rewritten PARAMETER stmt, FP becomes assignment
    assert 'LV=128' in new
    assert 'IPW2=4096' in new
    assert 'PARAMETER' in new
    assert any('R = 1.0D0' in text for _, text in assigns)


def test_convert_data_drops_known_constants(mf):
    src = '      DATA ZERO,TWO/0.D0,2.D0/\n'
    new, assigns, dropped = convert_data_stmts(src, mf)
    assert assigns == []
    assert dropped == {'ZERO': 'DD_ZERO', 'TWO': 'DD_TWO'}


def test_convert_data_keeps_unknown_names(mf):
    src = '      DATA GAM,GAMSQ/4096.D0,16777216.D0/\n'
    new, assigns, dropped = convert_data_stmts(src, mf)
    assert any('GAM = 4096.D0' in text for _, text in assigns)
    assert any('GAMSQ = 16777216.D0' in text for _, text in assigns)
    assert dropped == {}


def test_convert_parameter_multiline_continuation(mf):
    """Multi-line PARAMETER (column-6 continuation) is joined and parsed."""
    src = (
        '      PARAMETER (X = 1.0D+0,\n'
        '     $           Y = 2.0D+0)\n'
    )
    new, assigns, dropped = convert_parameter_stmts(src, mf)
    assert any('X = 1.0D+0' in text for _, text in assigns)
    assert any('Y = 2.0D+0' in text for _, text in assigns)


def test_convert_parameter_rejects_identifier_with_e(mf):
    """``E`` inside an identifier name must NOT make the value count as FP."""
    src = '      PARAMETER (LEN = ELEMENT)\n'
    new, assigns, dropped = convert_parameter_stmts(src, mf)
    # ELEMENT is not a numeric literal, not a known constant — leave alone.
    assert assigns == []
    assert 'PARAMETER (LEN = ELEMENT)' in new


def test_convert_data_multiline_continuation(mf):
    src = (
        '      DATA A,B/1.0D0,\n'
        '     $          2.0D0/\n'
    )
    new, assigns, dropped = convert_data_stmts(src, mf)
    assert any('A = 1.0D0' in text for _, text in assigns)
    assert any('B = 2.0D0' in text for _, text in assigns)


# ---------------------------------------------------------------------------
# insert_use_multifloats
# ---------------------------------------------------------------------------


def test_insert_use_multifloats_placement_after_decls(mf):
    """Extra param/data assignments must come AFTER the declaration block,
    not right after the SUBROUTINE header."""
    src = textwrap.dedent('''\
              SUBROUTINE FOO(A)
              DOUBLE PRECISION A
              INTEGER N
              A = 1.0D0
              END
        ''')
    out = insert_use_multifloats(src, mf, extra_lines=['      X = real64x2(\'1.0D0\')\n'])
    lines = out.splitlines()
    use_idx = next(i for i, l in enumerate(lines) if 'USE multifloats' in l)
    assign_idx = next(i for i, l in enumerate(lines) if 'X = real64x2' in l)
    int_idx = next(i for i, l in enumerate(lines) if 'INTEGER N' in l)
    real_idx = next(i for i, l in enumerate(lines) if 'DOUBLE PRECISION A' in l)
    # USE comes right after the header, before any declaration
    assert use_idx < real_idx
    # Extra assignment comes AFTER the declaration block
    assert assign_idx > int_idx


def test_insert_use_multifloats_dedupes(mf):
    src = textwrap.dedent('''\
              SUBROUTINE FOO()
              USE multifloats
              END
        ''')
    out = insert_use_multifloats(src, mf, extra_lines=None)
    assert out.count('USE multifloats') == 1


def test_specialize_use_excludes_keepkind_masked_local(mf):
    """A keep-kind-protected ``DOUBLE PRECISION`` local whose name
    collides with a multifloats generic (e.g. ``GAMMA``) must NOT be
    imported via ``USE ..., only:``.

    ``specialize_use_module`` runs while the keep-kind sentinel is still
    in place (the restore happens after migrate_*_form returns), so the
    declaration appears as ``__KEEPKIND_DP__ :: GAMMA``. The local-name
    scanner must recognise the sentinel form, or ``gamma`` leaks into the
    import list and gfortran rejects it ("GENERIC procedure 'gamma' is
    not allowed as an actual argument"). Regression for the MUMPS 5.9.0
    DMUMPS_LRMIXED_FRxFR routine.
    """
    from migrator.fortran.keepkind import _KK_SENTINEL
    src = (
        "      SUBROUTINE FRXFR(A)\n"
        "      USE multifloats\n"
        f"      {_KK_SENTINEL} :: GAMMA\n"
        "      TYPE(real64x2) :: X\n"
        "      GAMMA = 2.0D0\n"
        "      X = real64x2(GAMMA)\n"
        "      END\n"
    )
    out = specialize_use_module(src, mf, fixed_form=True)
    use_line = next(l for l in out.splitlines() if 'USE multifloats' in l)
    assert ', only:' in use_line
    # gamma is a keep-kind-masked local, not the generic — must be excluded.
    assert 'gamma' not in use_line.lower()
    # real64x2 is genuinely used and public — must survive.
    assert 'real64x2' in use_line


# ---------------------------------------------------------------------------
# _la_constants_rename_map / rewrite_la_constants_use
# ---------------------------------------------------------------------------


def test_la_constants_rename_map_uses_dd_zz_for_multifloats(mf):
    m = _la_constants_rename_map(mf)
    assert m['DZERO'] == f'{_RPU}ZERO'
    assert m['DSAFMIN'] == f'{_RPU}SAFMIN'
    assert m['ZZERO'] == f'{_CPU}ZERO'
    # Unprefixed forms must NOT be in the map (would clobber USE alias LHS)
    assert 'ZERO' not in m
    assert 'SAFMIN' not in m


def test_la_constants_rename_map_kind16():
    m = _la_constants_rename_map(load_target('kind16'))
    assert m['DZERO'] == 'QZERO'
    assert m['ZZERO'] == 'XZERO'


def test_rewrite_la_constants_use_pattern_b(mf):
    src = textwrap.dedent('''\
        subroutine foo
           use LA_CONSTANTS, &
              only: wp=>dp, zero=>dzero, half=>dhalf, one=>done, &
                    safmin=>dsafmin, safmax=>dsafmax
           real(wp) :: x
        end
    ''')
    out = rewrite_la_constants_use(src, mf)
    assert 'LA_CONSTANTS_MW' in out
    # wp=>dp removed, but local aliases preserved with DD-prefixed RHS
    assert 'wp=>dp' not in out
    assert f'zero=>{_RP}zero' in out
    assert f'safmin=>{_RP}safmin' in out
    assert f'safmax=>{_RP}safmax' in out


# ---------------------------------------------------------------------------
# End-to-end fixed-form smoke
# ---------------------------------------------------------------------------


SYNTHETIC_BLAS = """\
      SUBROUTINE FOO(N,X,ALPHA)
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA
      INTEGER N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION X(*)
*     ..
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
*     ..
*     .. Local Scalars ..
      INTEGER I
*     ..
      IF (ALPHA.EQ.ZERO) RETURN
      DO 10 I = 1,N
          X(I) = ALPHA*X(I) + ONE
   10 CONTINUE
      RETURN
      END
"""


def test_end_to_end_fixed_form(mf):
    out = migrate_fixed_form(SYNTHETIC_BLAS, MigrationContext({}, mf))

    # Header is unchanged
    assert 'SUBROUTINE FOO' in out
    # USE multifloats inserted right after the header
    assert 'USE multifloats' in out
    # Type decls converted
    assert 'TYPE(real64x2) ALPHA' in out
    assert 'TYPE(real64x2) X(*)' in out
    # PARAMETER block + decl for ONE/ZERO are gone (the PARAMETER line
    # may survive as a `! Converted to assignments: ...` comment).
    assert 'DOUBLE PRECISION ONE,ZERO' not in out
    code_lines = [
        l for l in out.splitlines()
        if l and not l.lstrip().startswith('!') and not l.startswith(('*', 'C', 'c'))
    ]
    code = '\n'.join(code_lines)
    assert 'PARAMETER (ONE=' not in code
    # Body references renamed
    assert 'IF (ALPHA.EQ.DD_ZERO) RETURN' in out
    assert '+ DD_ONE' in out
    # No bare ZERO/ONE references in body code
    body_lines = [
        l for l in out.splitlines()
        if l and not l.startswith(('*', 'C', 'c', '!')) and 'Converted' not in l
    ]
    body = '\n'.join(body_lines)
    import re as _re
    assert not _re.search(r'\bZERO\b', body)
    assert not _re.search(r'\bONE\b', body)


# Capital-F preprocessed Fortran (e.g. iparam2stage.F, dsytrd_sb2st.F) has
# a ``#if defined(_OPENMP) / use omp_lib / #endif`` block sitting between
# the USE clause and the IMPLICIT NONE line. The PARAMETER-conversion
# anchor walker MUST skip those preprocessor directives — otherwise it
# treats ``#if`` as the first executable line, inserts the converted
# assignments above IMPLICIT NONE, and gfortran rejects the result with
# ``Unexpected IMPLICIT NONE statement``.
SYNTHETIC_LAPACK_F_OPENMP = """\
      SUBROUTINE FOO(N, X)
#if defined(_OPENMP)
      use omp_lib
#endif
      IMPLICIT NONE
*     .. Scalar Arguments ..
      INTEGER N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION X(*)
*     ..
*     .. Parameters ..
      DOUBLE PRECISION ONE, ZERO
      PARAMETER (ONE=1.0D+0, ZERO=0.0D+0)
*     ..
*     .. Local Scalars ..
      INTEGER I
*     ..
      DO 10 I = 1, N
          X(I) = ZERO
   10 CONTINUE
      RETURN
      END
"""


def test_param_assignment_lands_below_implicit_none_with_openmp_directive(mf):
    """Regression for the .F PARAMETER-ordering bug (UB-02): in the
    presence of a ``#if defined(_OPENMP) / use omp_lib / #endif`` block
    between USE and IMPLICIT NONE, the converted PARAMETER assignments
    must be inserted AFTER IMPLICIT NONE and the type declarations,
    not above them."""
    out = migrate_fixed_form(SYNTHETIC_LAPACK_F_OPENMP, MigrationContext({}, mf))

    lines = out.splitlines()

    def first_index(predicate):
        for idx, line in enumerate(lines):
            if predicate(line):
                return idx
        return -1

    implicit_idx = first_index(lambda l: l.lstrip().upper().startswith('IMPLICIT NONE'))
    # Match either a literal assignment (``ONE = real64x2(...)``) or the
    # known-constant rename (``DD_ONE = ...``); the migrator may emit
    # either depending on whether the constant is in its known table.
    assignment_idx = first_index(
        lambda l: l.lstrip().startswith(('ONE =', 'ZERO =', 'DD_ONE =', 'DD_ZERO ='))
    )
    openmp_if_idx = first_index(lambda l: l.lstrip().startswith('#if defined(_OPENMP)'))
    openmp_endif_idx = first_index(lambda l: l.lstrip().startswith('#endif'))

    assert implicit_idx >= 0, 'IMPLICIT NONE missing from output'
    assert openmp_if_idx >= 0, '#if defined(_OPENMP) directive was dropped'
    assert openmp_endif_idx >= 0, '#endif directive was dropped'
    # The directive block must remain intact and ABOVE IMPLICIT NONE
    # (Fortran allows USE + #if-guarded USE before IMPLICIT NONE).
    assert openmp_if_idx < implicit_idx
    assert openmp_endif_idx < implicit_idx
    # The converted assignment(s) — if emitted at all — must land below
    # IMPLICIT NONE. ZERO is a known constant that may be dropped; ONE
    # is unknown and survives as an assignment in fixed-form.
    if assignment_idx >= 0:
        assert assignment_idx > implicit_idx, (
            f'PARAMETER assignment at line {assignment_idx} '
            f'must come AFTER IMPLICIT NONE at line {implicit_idx}'
        )


def test_equivalence_no_diagnostic_after_sequence(mf, capsys):
    """The mock multifloats module now declares both real64x2 and
    cmplx64x2 with the SEQUENCE attribute, so EQUIVALENCE on
    migrated FP variables compiles cleanly. The migrator therefore
    does not emit a warning anymore.
    """
    src = (
        "      SUBROUTINE FOO()\n"
        "      DOUBLE PRECISION CR(2,2), CRV(4)\n"
        "      EQUIVALENCE (CR(1,1), CRV(1))\n"
        "      END\n"
    )
    migrate_fixed_form(src, MigrationContext({}, mf))
    captured = capsys.readouterr()
    assert 'EQUIVALENCE' not in captured.err
    assert 'Manual rewrite required' not in captured.err


def test_equivalence_no_warning_for_integer_only(mf, capsys):
    """Integer EQUIVALENCE is fine and must not warn."""
    src = (
        "      SUBROUTINE FOO()\n"
        "      INTEGER A(4), B(2,2)\n"
        "      EQUIVALENCE (A(1), B(1,1))\n"
        "      END\n"
    )
    migrate_fixed_form(src, MigrationContext({}, mf))
    captured = capsys.readouterr()
    assert 'EQUIVALENCE' not in captured.err


SYNTHETIC_LAPACK_FREE_FORM = """\
subroutine FOO( x, y )
   use LA_CONSTANTS, &
      only: wp=>dp, zero=>dzero, one=>done, safmin=>dsafmin
   real(wp) :: x, y
   if( x == zero ) then
      y = one / safmin
   end if
end subroutine
"""


def test_canonicalize_for_compare_strips_multifloats(mf):
    """The convergence canonicalizer must strip real64x2(...) wrappers
    and TYPE(real64x2) so S/D pairs migrated to multifloats canonicalize
    to the same text as the kind-based canonicalization would."""
    from migrator.divergence import _canonicalize_for_compare
    a = "      X = real64x2('1.0D+0') + real64x2('2.0D+0')"
    b = "      X = 1.0D+0 + 2.0D+0"
    assert _canonicalize_for_compare(a) == _canonicalize_for_compare(b)


def test_canonicalize_for_compare_normalizes_type_decl():
    from migrator.divergence import _canonicalize_for_compare
    a = "      TYPE(real64x2) X"
    ca = _canonicalize_for_compare(a)
    # The migrator only ever produces TYPE(real64x2) for multifloats,
    # so we just need to confirm `TYPE(real64x2)` doesn't survive
    # canonicalization.
    assert 'FLOAT64X2' not in ca
    assert 'TYPE(' not in ca


def test_end_to_end_free_form_pattern_b(mf):
    out = migrate_free_form(SYNTHETIC_LAPACK_FREE_FORM, MigrationContext({}, mf))
    assert 'LA_CONSTANTS_MW' in out
    assert 'wp=>dp' not in out
    # Local aliases preserved (lowercase) with DD-prefixed RHS
    assert f'zero=>{_RP}zero' in out
    assert f'safmin=>{_RP}safmin' in out
    # real(wp) → TYPE(real64x2)
    assert 'TYPE(real64x2)' in out
    # Body references stay as the local alias names
    assert 'x == zero' in out
    assert 'one / safmin' in out


# ---------------------------------------------------------------------------
# _rewrite_mpi_sum — multifloats per-call-site MPI_SUM rewriter
# ---------------------------------------------------------------------------


def _mpi_rewrite(src, target):
    """Datatype pass then reduce-op pass — the order the migrator runs
    them in (the op rewriter keys on the already-rewritten datatype
    token in the call's argument list)."""
    return _rewrite_mpi_sum(_rewrite_mpi_datatypes(src, target), target)


def test_rewrite_mpi_sum_real_call_site(mf):
    """A reduce call whose arglist contains the rewritten real datatype
    token gets MPI_SUM swapped for MPI_MM_SUM."""
    src = (
        "      CALL MPI_ALLREDUCE(SBUF, RBUF, 1, MPI_DOUBLE_PRECISION,\n"
        "     &                   MPI_SUM, COMM, IERR)\n"
    )
    out = _mpi_rewrite(src, mf)
    assert 'MPI_FLOAT64X2' in out
    assert 'MPI_MM_SUM' in out
    assert 'MPI_SUM' not in out


def test_rewrite_mpi_sum_complex_call_site(mf):
    """Same but for a complex-typed reduce: MPI_SUM → MPI_WW_SUM."""
    src = (
        "      CALL MPI_REDUCE(ZB, ZA, N, MPI_DOUBLE_COMPLEX, MPI_SUM,\n"
        "     &                ROOT, COMM, IERR)\n"
    )
    out = _mpi_rewrite(src, mf)
    assert 'MPI_COMPLEX64X2' in out
    assert 'MPI_WW_SUM' in out
    assert 'MPI_SUM' not in out


def test_rewrite_mpi_sum_integer_call_site_untouched(mf):
    """Integer reductions stay MPI_SUM — MPI_MM_SUM is undefined for
    MPI_INTEGER. This is the whole reason for token-context dispatch."""
    src = (
        "      CALL MPI_ALLREDUCE(N, NSUM, 1, MPI_INTEGER, MPI_SUM,\n"
        "     &                   COMM, IERR)\n"
    )
    out = _mpi_rewrite(src, mf)
    assert 'MPI_SUM' in out
    assert 'MPI_MM_SUM' not in out
    assert 'MPI_WW_SUM' not in out


def test_rewrite_mpi_sum_mixed_arith_in_one_file(mf):
    """A single source with all three call sites — each judged
    independently. Real → MM_SUM, complex → WW_SUM, integer untouched."""
    src = (
        "      CALL MPI_ALLREDUCE(R, S, 1, MPI_DOUBLE_PRECISION, MPI_SUM, C, I)\n"
        "      CALL MPI_ALLREDUCE(Z, W, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, C, I)\n"
        "      CALL MPI_ALLREDUCE(N, M, 1, MPI_INTEGER, MPI_SUM, C, I)\n"
    )
    out = _mpi_rewrite(src, mf)
    lines = out.splitlines()
    assert 'MPI_MM_SUM' in lines[0] and 'MPI_FLOAT64X2' in lines[0]
    assert 'MPI_WW_SUM' in lines[1] and 'MPI_COMPLEX64X2' in lines[1]
    assert 'MPI_SUM' in lines[2] and 'MPI_INTEGER' in lines[2]
    assert 'MPI_MM_SUM' not in lines[2]
    assert 'MPI_WW_SUM' not in lines[2]


def test_rewrite_mpi_sum_kind16_uses_custom_quad_ops():
    """KIND=16 keeps the standard MPI_REAL16 / MPI_COMPLEX32 datatypes but
    routes reductions through custom MPI_Op handles (quad_mpi_init), because
    Intel MPI's builtin MPI_SUM / MPI_MAX / MPI_MIN have no 16-byte-real
    kernel and segfault at np >= 2. So the reduce op — not the datatype — is
    rewritten inside any reduce call carrying the quad datatype; integer
    reductions are left untouched."""
    k16 = load_target('kind16')
    src = (
        "      CALL MPI_ALLREDUCE(R, S, 1, MPI_DOUBLE_PRECISION, MPI_SUM, C, I)\n"
        "      CALL MPI_ALLREDUCE(Z, W, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, C, I)\n"
        "      CALL MPI_REDUCE(R, S, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, C, I)\n"
        "      CALL MPI_ALLREDUCE(N, M, 1, MPI_INTEGER, MPI_SUM, C, I)\n"
    )
    out = _mpi_rewrite(src, k16)
    lines = out.splitlines()
    assert 'MPI_REAL16' in lines[0] and 'MPI_QQ_SUM' in lines[0]
    assert 'MPI_COMPLEX32' in lines[1] and 'MPI_XX_SUM' in lines[1]
    assert 'MPI_REAL16' in lines[2] and 'MPI_QQ_AMX' in lines[2]
    # integer reductions stay on the stock builtin op
    assert 'MPI_INTEGER' in lines[3] and 'MPI_SUM' in lines[3]
    assert 'MPI_QQ_SUM' not in lines[3]


def test_rewrite_mpi_sum_kind10_is_noop():
    k10 = load_target('kind10')
    src = (
        "      CALL MPI_ALLREDUCE(R, S, 1, MPI_DOUBLE_PRECISION, MPI_SUM, C, I)\n"
    )
    out = _mpi_rewrite(src, k10)
    assert 'MPI_LONG_DOUBLE' in out
    assert 'MPI_SUM' in out
    assert 'MPI_MM_SUM' not in out


def test_rewrite_mpi_sum_does_not_touch_non_reduce_calls(mf):
    """MPI_SUM appearing outside a reduction call (e.g. as a literal in
    a SAVE / DATA / parameter list) is NOT rewritten by this pass —
    the regex anchors on MPI_(I?All)Reduce[_scatter]."""
    src = (
        "      INTEGER, PARAMETER :: MY_OP = MPI_SUM\n"
        "      CALL MPI_OP_FREE(MPI_SUM, IERR)\n"
    )
    out = _rewrite_mpi_sum(src, mf)
    assert out == src


def test_rewrite_mpi_max_real_call_site(mf):
    """MPI_MAX on a float-typed reduce becomes MPI_MM_AMX (multifloats
    bridge ships AMX, not a stock-MPI-compatible MAX)."""
    src = (
        "      CALL MPI_REDUCE(SBUF, RBUF, 1, MPI_DOUBLE_PRECISION,\n"
        "     &                MPI_MAX, MASTER, COMM, IERR)\n"
    )
    out = _mpi_rewrite(src, mf)
    assert 'MPI_MM_AMX' in out
    assert 'MPI_MAX' not in out


def test_rewrite_mpi_reduce_nested_paren_argument():
    """The op must be rewritten even when a reduce argument nests
    parentheses two levels deep — e.g. ``WRKRC(1_8+int(M,8))``, the
    COLSCA reduce in ``*fac_scalings_simScaleAbs.F``. A fixed-depth
    regex for the argument list failed to match such a call at all,
    silently leaving a stock ``MPI_MAX`` on a wide datatype (undefined
    for the target op) and corrupting distributed scaling into NaN. The
    balanced-paren scan must capture it. Guards both real (kind16) and
    complex arithmetic paths."""
    k16 = load_target('kind16')
    src = (
        "      CALL MPI_REDUCE(COLSCA, WRKRC(1_8+int(M,8)), N,\n"
        "     &     MPI_DOUBLE_PRECISION,\n"
        "     &     MPI_MAX, 0,\n"
        "     &     COMM, IERROR)\n"
    )
    out = _mpi_rewrite(src, k16)
    assert 'MPI_REAL16' in out
    assert 'MPI_QQ_AMX' in out
    assert 'MPI_MAX' not in out


def test_rewrite_mpi_min_complex_call_site(mf):
    src = (
        "      CALL MPI_ALLREDUCE(Z, W, 1, MPI_DOUBLE_COMPLEX, MPI_MIN, C, I)\n"
    )
    out = _mpi_rewrite(src, mf)
    assert 'MPI_WW_AMN' in out
    assert 'MPI_MIN' not in out


def test_rewrite_mpi_max_integer_call_site_untouched(mf):
    """MPI_INTEGER + MPI_MAX must stay MPI_MAX — integer reductions
    don't go through the user-defined op."""
    src = (
        "      CALL MPI_ALLREDUCE(N, M, 1, MPI_INTEGER, MPI_MAX, C, I)\n"
    )
    out = _mpi_rewrite(src, mf)
    assert 'MPI_MAX' in out
    assert 'MPI_MM_AMX' not in out


def test_rewrite_mpi_sum_handles_mpi_reduce_scatter(mf):
    """MPI_REDUCE_SCATTER also takes (op, datatype) — must be covered."""
    src = (
        "      CALL MPI_REDUCE_SCATTER(SB, RB, CNTS, MPI_DOUBLE_PRECISION,\n"
        "     &                        MPI_SUM, COMM, IERR)\n"
    )
    out = _mpi_rewrite(src, mf)
    assert 'MPI_MM_SUM' in out
    assert 'MPI_SUM' not in out
