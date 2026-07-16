"""Tests for the empirical precision-prefix classifier.

The classifier walks a symbol set and discovers ``S/D/C/Z`` precision
families purely from name shape. These tests pin the classifier on a
representative subset of BLAS / LAPACK / ScaLAPACK symbols and the
known cross-type / doubled-prefix edge cases.
"""

from migrator.prefix_classifier import (
    CHAR_TYPE,
    classify_symbols,
    build_rename_map,
)
from migrator.target_mode import load_target


def test_char_type_table():
    """S/D map to real, C/Z map to complex; nothing else is recognized."""
    assert CHAR_TYPE == {'S': 'R', 'D': 'R', 'C': 'C', 'Z': 'C'}


def test_blas_gemm_family():
    cls = classify_symbols({'SGEMM', 'DGEMM', 'CGEMM', 'ZGEMM'})
    real_fam = cls.get_family('DGEMM')
    cplx_fam = cls.get_family('ZGEMM')
    assert real_fam is not None and cplx_fam is not None
    # Real and complex are *different* families even though both
    # have a single prefix slot at position 0.
    assert real_fam is not cplx_fam
    # SGEMM is the S-half of DGEMM's family; ZGEMM is the Z-half.
    assert real_fam is cls.get_family('SGEMM')
    assert cplx_fam is cls.get_family('CGEMM')


def test_independent_symbols():
    """Symbols with no precision-letter in any candidate position
    should land in ``independent``, not in any family."""
    cls = classify_symbols({'XERBLA', 'LSAME', 'IEEECK', 'SGEMM', 'DGEMM'})
    assert 'XERBLA' in cls.independent
    assert 'LSAME' in cls.independent
    assert cls.get_family('XERBLA') is None


def test_scalapack_p_prefix():
    """P-prefixed ScaLAPACK symbols pair on the *second* character."""
    cls = classify_symbols({'PSGESV', 'PDGESV', 'PCGESV', 'PZGESV'})
    fam = cls.get_family('PDGESV')
    assert fam is not None
    assert fam is cls.get_family('PSGESV')
    # The complex sibling forms a separate family.
    assert cls.get_family('PZGESV') is not fam


def test_cross_type_csrot_zdrot():
    """``CSROT`` / ``ZDROT`` pair via a 2-slot ``#C#R ROT`` pattern."""
    cls = classify_symbols({'CSROT', 'ZDROT'})
    fam = cls.get_family('CSROT')
    assert fam is not None
    assert fam is cls.get_family('ZDROT')
    # Two slots, one complex (pos 0) and one real (pos 1).
    assert sorted(s.position for s in fam.slots) == [0, 1]


def test_doubled_prefix_zzdotc_ccdotc():
    """ScaLAPACK's ZZDOTC/CCDOTC use both leading letters as the prefix
    slot — only the 2-position pattern can pair them."""
    cls = classify_symbols({'CCDOTC', 'ZZDOTC'})
    fam = cls.get_family('CCDOTC')
    assert fam is not None
    assert fam is cls.get_family('ZZDOTC')
    assert sorted(s.position for s in fam.slots) == [0, 1]


def test_dsytrs_dsytrd_not_grouped_via_trailing_letter():
    """DSYTRS and DSYTRD share a trailing S/D but should NOT be grouped
    on it — the first-position prefix family wins via the position
    tiebreaker."""
    cls = classify_symbols({
        'SSYTRS', 'DSYTRS',  # legitimate -TRS family (slot at pos 0)
        'SSYTRD', 'DSYTRD',  # legitimate -TRD family (slot at pos 0)
    })
    trs = cls.get_family('DSYTRS')
    trd = cls.get_family('DSYTRD')
    assert trs is not None and trd is not None
    assert trs is not trd
    # Each family's slot is at position 0, not the trailing S/D.
    for fam in (trs, trd):
        assert [s.position for s in fam.slots] == [0]


def test_build_rename_map_kind16():
    """End-to-end: kind16 maps DGEMM→QGEMM, ZGEMM→XGEMM, leaves S/C
    halves and independents alone."""
    target = load_target('kind16')
    rename = build_rename_map(
        {'SGEMM', 'DGEMM', 'CGEMM', 'ZGEMM', 'XERBLA'},
        target,
    )
    # Every family member maps to the same target name. D→Q, Z→X
    # under kind16; S and C halves co-converge to those same targets
    # because they share the family. (The convergence check happens
    # at migration time; same target name is by design.)
    assert rename['DGEMM'] == 'QGEMM'
    assert rename['SGEMM'] == 'QGEMM'
    assert rename['ZGEMM'] == 'XGEMM'
    assert rename['CGEMM'] == 'XGEMM'
    # Independents are absent from the rename map.
    assert 'XERBLA' not in rename


def test_no_collision_for_well_formed_input():
    """The collision diagnostic must not raise on a clean BLAS subset."""
    target = load_target('kind16')
    # No exception expected.
    build_rename_map({'SGEMM', 'DGEMM', 'CGEMM', 'ZGEMM'}, target)
