"""Unit tests for the derived-type component-type oracle
(``migrator.fortran.decls.scan_type_component_names``) and the pipeline
helper that unions it across the staged tree
(``migrator.pipeline._build_component_oracle``).

The oracle exists so formatted-output narrowing can recognise a struct
component reached as ``id%DKEEP(160)`` in a file that only ``USE``s the
struct-definition module — the per-file variable oracle cannot see it.
"""

from migrator.fortran.decls import scan_type_component_names


_STRUC = """\
      TYPE DMUMPS_STRUC
        INTEGER :: N, NNZ
        DOUBLE PRECISION :: CNTL(15)
        DOUBLE PRECISION :: RINFOG(40)
        DOUBLE PRECISION :: DKEEP(230)
        DOUBLE PRECISION, DIMENSION(:), POINTER :: A
        INTEGER, DIMENSION(:), POINTER :: IRN
      END TYPE DMUMPS_STRUC
"""

# The c/z struct: the same data-array component A is COMPLEX here.
_STRUC_C = """\
      TYPE CMUMPS_STRUC
        DOUBLE PRECISION :: CNTL(15)
        DOUBLE PRECISION :: RINFOG(40)
        COMPLEX, DIMENSION(:), POINTER :: A
      END TYPE CMUMPS_STRUC
"""


def test_scan_harvests_real_components():
    real, cplx = scan_type_component_names(_STRUC)
    assert {'CNTL', 'RINFOG', 'DKEEP', 'A'} <= real
    assert cplx == set()
    # Integer fields are not real-typed data.
    assert 'N' not in real and 'IRN' not in real


def test_scan_harvests_complex_components():
    real, cplx = scan_type_component_names(_STRUC_C)
    assert 'A' in cplx
    assert {'CNTL', 'RINFOG'} <= real


def test_scan_empty_when_no_derived_type():
    src = """\
      SUBROUTINE FOO(X)
      DOUBLE PRECISION :: X
      END SUBROUTINE
"""
    assert scan_type_component_names(src) == (set(), set())


def test_type_variable_decl_not_a_definition():
    # ``TYPE(...)`` is a variable declaration, not a definition block —
    # its "body" must not be harvested.
    src = """\
      SUBROUTINE FOO(ID)
      TYPE(DMUMPS_STRUC) :: ID
      DOUBLE PRECISION :: LOCAL
      END SUBROUTINE
"""
    assert scan_type_component_names(src) == (set(), set())


def test_union_drops_real_intersect_complex():
    # Mirror the pipeline reduction: A is real in the d-struct and complex
    # in the c-struct, so it is ambiguous and dropped from both sets; the
    # printed statistics survive in comp_real.
    r = set()
    c = set()
    for src in (_STRUC, _STRUC_C):
        rr, cc = scan_type_component_names(src)
        r |= rr
        c |= cc
    ambiguous = r & c
    comp_real = r - ambiguous
    comp_complex = c - ambiguous
    assert 'A' in ambiguous
    assert 'A' not in comp_real and 'A' not in comp_complex
    assert {'CNTL', 'RINFOG', 'DKEEP'} <= comp_real
