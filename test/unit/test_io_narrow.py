"""Unit tests for the formatted-output narrowing pass
(``migrator.fortran.io_narrow``).

The pass wraps direct real64x2 / cmplx64x2 references in a formatted
WRITE/PRINT output list as ``dble(x)`` / ``cmplx(dble(x%re), ...)`` so a
single real edit descriptor consumes a single value (see the module
docstring for the derived-type-expansion crash it fixes). These tests
exercise the single-statement path directly and the cpp-interrupted
continuation case through ``migrate_fixed_form``.
"""

from migrator.target_mode import load_target
from migrator.fortran_migrator import migrate_fixed_form, MigrationContext
from migrator.fortran.io_narrow import (
    narrow_multifloats_io,
    narrow_multifloats_io_open,
    narrow_io_continuation,
    is_fixed_io_continuation,
)

REAL = {'RINFOG', 'RINFO'}
CPLX = {'CX'}


# --- single-statement narrowing -------------------------------------------

def test_formatted_write_narrows_real_reference():
    line = '      WRITE(MPG, 99992) KEEP(1), RINFOG(1)'
    assert narrow_multifloats_io(line, REAL, CPLX) == \
        '      WRITE(MPG, 99992) KEEP(1), dble(RINFOG(1))'


def test_formatted_write_leaves_integer_items():
    line = '      WRITE(6, 10) KEEP(1), INFOG(4), ITMP'
    assert narrow_multifloats_io(line, REAL, CPLX) == line


def test_whole_array_reference_narrowed():
    line = '      WRITE(6, 10) RINFO'
    assert narrow_multifloats_io(line, REAL, CPLX) == \
        '      WRITE(6, 10) dble(RINFO)'


def test_complex_reference_expands_to_cmplx():
    line = '      WRITE(6, 10) CX(2)'
    assert narrow_multifloats_io(line, REAL, CPLX) == \
        '      WRITE(6, 10) cmplx(dble(CX(2)%re), dble(CX(2)%im), kind=8)'


def test_print_statement_narrowed():
    line = '      PRINT 10, RINFOG(1)'
    assert narrow_multifloats_io(line, REAL, CPLX) == \
        '      PRINT 10, dble(RINFOG(1))'


def test_list_directed_write_narrowed():
    line = '      WRITE(6, *) RINFOG(1)'
    assert narrow_multifloats_io(line, REAL, CPLX) == \
        '      WRITE(6, *) dble(RINFOG(1))'


def test_logical_if_prefix_narrowed():
    line = '      IF (K.GT.1) WRITE(6, 10) RINFOG(2)'
    assert narrow_multifloats_io(line, REAL, CPLX) == \
        '      IF (K.GT.1) WRITE(6, 10) dble(RINFOG(2))'


def test_unformatted_write_untouched():
    # Unit only — no format spec — is unformatted; never narrow.
    line = '      WRITE(21) RINFOG(1)'
    assert narrow_multifloats_io(line, REAL, CPLX) == line


def test_iostat_keyword_write_untouched():
    line = '      WRITE(21, IOSTAT=IERR) RINFOG(1)'
    assert narrow_multifloats_io(line, REAL, CPLX) == line


def test_read_never_narrowed():
    line = '      READ(6, 10) RINFOG(1)'
    assert narrow_multifloats_io(line, REAL, CPLX) == line


def test_expression_item_untouched():
    # An expression is not a pure reference: narrowing it could feed a
    # non-real value to a real descriptor. Leave it alone.
    line = '      WRITE(6, 10) RINFOG(1) + RINFO(2)'
    assert narrow_multifloats_io(line, REAL, CPLX) == line


def test_intrinsic_call_item_untouched():
    line = '      WRITE(6, 10) INT(RINFOG(1))'
    assert narrow_multifloats_io(line, REAL, CPLX) == line


def test_string_literal_not_confused_for_reference():
    line = "      WRITE(6, '(A)') 'RINFOG(1)'"
    assert narrow_multifloats_io(line, REAL, CPLX) == line


def test_noop_when_name_oracle_empty():
    line = '      WRITE(6, 10) RINFOG(1)'
    assert narrow_multifloats_io(line, set(), set()) == line


def test_open_flag_true_for_formatted_write():
    line = '      WRITE(6, 10) KEEP(1),'
    _, opened = narrow_multifloats_io_open(line, REAL, CPLX)
    assert opened is True


def test_open_flag_false_for_non_output():
    line = '      X = RINFOG(1)'
    _, opened = narrow_multifloats_io_open(line, REAL, CPLX)
    assert opened is False


# --- global derived-type component oracle ---------------------------------

# Component names known only through the struct definition in a USE'd module
# (e.g. id%DKEEP(160) in dsol_driver.F). The per-file oracle can't see them,
# so they arrive via the comp_real / comp_complex sets and narrow ONLY when
# reached through a ``%``-qualified reference.
COMP_REAL = {'DKEEP', 'RINFOG'}
COMP_CPLX = {'LPS'}


def test_qualified_component_reference_narrowed():
    line = '      WRITE( MPG, 434 ) id%DKEEP(160)'
    assert narrow_multifloats_io(line, set(), set(), COMP_REAL, COMP_CPLX) == \
        '      WRITE( MPG, 434 ) dble(id%DKEEP(160))'


def test_bare_component_name_not_narrowed():
    # A bare DKEEP is a different (local) entity — the component oracle must
    # only fire on a %-qualified reference, never on the bare name.
    line = '      WRITE(6, 10) DKEEP(160)'
    assert narrow_multifloats_io(line, set(), set(), COMP_REAL, COMP_CPLX) == line


def test_qualified_complex_component_expands_to_cmplx():
    line = '      WRITE(6, 10) id%LPS(2)'
    assert narrow_multifloats_io(line, set(), set(), COMP_REAL, COMP_CPLX) == \
        '      WRITE(6, 10) cmplx(dble(id%LPS(2)%re), dble(id%LPS(2)%im), kind=8)'


def test_qualified_component_not_in_oracle_untouched():
    line = '      WRITE(6, 10) id%RHS(1)'
    assert narrow_multifloats_io(line, set(), set(), COMP_REAL, COMP_CPLX) == line


def test_component_oracle_noop_when_all_empty():
    line = '      WRITE(6, 10) id%DKEEP(160)'
    assert narrow_multifloats_io(line, set(), set()) == line


def test_component_oracle_via_continuation():
    frag = '     &  KEEP(1), id%RINFOG(1)'
    assert narrow_io_continuation(frag, set(), set(), COMP_REAL, COMP_CPLX) == \
        '     &  KEEP(1), dble(id%RINFOG(1))'


# --- fixed-form continuation helpers --------------------------------------

def test_is_fixed_io_continuation():
    assert is_fixed_io_continuation('     &  KEEP(1), RINFOG(1)')
    assert is_fixed_io_continuation('     +  KEEP(1)')
    assert not is_fixed_io_continuation('      WRITE(6, 10) RINFOG(1)')
    assert not is_fixed_io_continuation('     0  RINFOG(1)')  # col-6 '0'


def test_continuation_narrowing_preserves_marker():
    frag = '     &  KEEP(56), KEEP(61), RINFOG(1)'
    assert narrow_io_continuation(frag, REAL, CPLX) == \
        '     &  KEEP(56), KEEP(61), dble(RINFOG(1))'


def test_continuation_not_narrowed_when_no_real():
    frag = '     &  KEEP(56), KEEP(61)'
    assert narrow_io_continuation(frag, REAL, CPLX) == frag


# --- end-to-end: cpp-interrupted WRITE continuation -----------------------

_CPP_WRITE_SRC = """\
      SUBROUTINE DIAG
      DOUBLE PRECISION, INTENT(IN) :: RINFOG(40)
      INTEGER, INTENT(IN) :: KEEP(500)
      WRITE(MPG, 99992) KEEP(1), KEEP(2),
     &  KEEP(3),
#if defined(FOO)
     &  KEEP(4),
#endif
     &  KEEP(56), KEEP(61), RINFOG(1)
99992 FORMAT('x =', I16, ' y =', 1PD10.3)
      END SUBROUTINE
"""


def test_cpp_interrupted_write_continuation_narrowed():
    tm = load_target('multifloats')
    out = migrate_fixed_form(_CPP_WRITE_SRC, MigrationContext({}, tm, source_kind=8))
    # The real64x2 item on the headless continuation fragment is narrowed,
    assert 'dble(RINFOG(1))' in out
    # the interleaved cpp directives survive intact,
    assert '#if defined(FOO)' in out
    assert '#endif' in out
    # and the USE clause imports the dble generic.
    assert 'dble' in out.split('IMPLICIT NONE')[0]


def test_cpp_interrupted_write_does_not_leak_state():
    # A following non-output statement that happens to reference a
    # continuation-shaped fragment must not be narrowed once the WRITE
    # list has closed.
    src = _CPP_WRITE_SRC.replace(
        "99992 FORMAT('x =', I16, ' y =', 1PD10.3)",
        "      RINFOG(2) = RINFOG(1)\n"
        "99992 FORMAT('x =', I16, ' y =', 1PD10.3)",
    )
    tm = load_target('multifloats')
    out = migrate_fixed_form(src, MigrationContext({}, tm, source_kind=8))
    # The assignment RHS/LHS must stay bare — only the WRITE item narrowed.
    assert 'dble(RINFOG(1))' in out
    assert out.count('dble(RINFOG') == 1
