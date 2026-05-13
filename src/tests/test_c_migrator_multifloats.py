"""Unit tests for the multifloats branch of migrator.c_migrator.

These cover the c_migrator changes that unblock multifloats:
- _build_sub_vars produces the right template variables for both
  KIND and multifloats targets;
- the BLACS regex substitution rules emit DD/ZZ-prefixed names plus
  the new MPI_SUM->MPI_DD_SUM/MPI_ZZ_SUM rule;
- _patch_bdef_header for multifloats emits an #include of the libmfc
  header instead of the legacy __float128 typedef block;
- the override_files post-pass copies hand-written replacement files
  on top of clones.
"""

import textwrap
from pathlib import Path

import pytest

from migrator.target_mode import load_target
from migrator.c_migrator import (
    _build_sub_vars,
    _apply_aliases_to_original,
    _apply_overrides,
    _build_rename_regex,
    _make_rename_substituter,
    _patch_bdef_header,
    migrate_c_file_to_string,
)


# Multifloats prefix letters, derived from the loaded target (so tests
# don't have to be edited every time the prefix is changed).
_MF = load_target('multifloats')
_RP = _MF.prefix_map['R'].lower()       # real-prefix, lower
_CP = _MF.prefix_map['C'].lower()       # complex-prefix, lower
_RPU = _MF.prefix_map['R']              # real-prefix, upper
_CPU = _MF.prefix_map['C']              # complex-prefix, upper


# ---------------------------------------------------------------------------
# _build_sub_vars
# ---------------------------------------------------------------------------


def test_build_sub_vars_multifloats_struct_types():
    mf = load_target('multifloats')
    v = _build_sub_vars(mf)
    assert v['REAL_TYPE'] == 'float64x2'
    assert v['COMPLEX_TYPE'] == 'complex64x2'
    assert v['C_REAL_TYPE'] == 'float64x2'
    assert v['MPI_REAL'] == 'MPI_FLOAT64X2'
    assert v['MPI_COMPLEX'] == 'MPI_COMPLEX64X2'
    assert v['MPI_SUM_REAL'] == 'MPI_DD_SUM'
    assert v['MPI_SUM_COMPLEX'] == 'MPI_ZZ_SUM'
    assert v['RP'] == _RP
    assert v['CP'] == _CP
    assert v['RPU'] == _RPU
    assert v['CPU'] == _CPU


def test_build_sub_vars_kind16_unchanged():
    """KIND targets keep the legacy values; the new MPI_SUM_* keys
    expand to plain 'MPI_SUM' so the substitution rule is a no-op."""
    k16 = load_target('kind16')
    v = _build_sub_vars(k16)
    assert v['REAL_TYPE'] == 'QREAL'
    assert v['COMPLEX_TYPE'] == 'QCOMPLEX'
    assert v['C_REAL_TYPE'] == '__float128'
    assert v['MPI_REAL'] == 'MPI_REAL16'
    assert v['MPI_COMPLEX'] == 'MPI_COMPLEX32'
    assert v['MPI_SUM_REAL'] == 'MPI_SUM'
    assert v['MPI_SUM_COMPLEX'] == 'MPI_SUM'
    assert v['RP'] == 'q'
    assert v['CP'] == 'x'


def test_build_sub_vars_kind10_unchanged():
    k10 = load_target('kind10')
    v = _build_sub_vars(k10)
    assert v['C_REAL_TYPE'] == 'long double'
    assert v['MPI_REAL'] == 'MPI_LONG_DOUBLE'
    assert v['MPI_SUM_REAL'] == 'MPI_SUM'


# ---------------------------------------------------------------------------
# Per-file BLACS substitution (REAL clone)
# ---------------------------------------------------------------------------


_BLACS_REAL_SAMPLE = '''\
#include "Bdef.h"
void BI_dvvsum(Int N, double *v1, double *v2) {
   for (int k = 0; k < N; k++) v1[k] += v2[k];
}

void Cdgebs2d_aux(double *buf, Int n) {
   MPI_Reduce(buf, buf, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}

void BI_dMPI_sum(double *a, double *b, int *len, MPI_Datatype *dt) {
   for (int i = 0; i < *len; i++) b[i] += a[i];
}
'''


def test_blacs_real_clone_multifloats(tmp_path):
    src = tmp_path / 'BI_dvvsum.c'
    src.write_text(_BLACS_REAL_SAMPLE)
    target = load_target('multifloats')
    result = migrate_c_file_to_string(src, target)
    assert result is not None
    new_name, text = result
    assert new_name == f'BI_{_RP}vvsum.c'
    assert 'float64x2 *v1' in text
    assert 'float64x2 *v2' in text
    assert 'MPI_FLOAT64X2' in text
    assert 'MPI_DOUBLE' not in text
    # MPI_SUM must be rewritten to the user-defined op
    assert 'MPI_DD_SUM' in text
    assert 'MPI_SUM,' not in text  # the old constant is gone
    # Routine name renames cover both Cd* and BI_d* prefixes
    assert f'C{_RP}gebs2d_aux' in text
    assert f'BI_{_RP}MPI_sum' in text


def test_blacs_real_clone_kind16_still_uses_mpi_sum(tmp_path):
    """Regression: the new MPI_SUM->{MPI_SUM_REAL} rule must be a
    textual no-op for KIND targets."""
    src = tmp_path / 'BI_dvvsum.c'
    src.write_text(_BLACS_REAL_SAMPLE)
    target = load_target('kind16')
    result = migrate_c_file_to_string(src, target)
    assert result is not None
    new_name, text = result
    assert new_name == 'BI_qvvsum.c'
    assert 'QREAL *v1' in text
    assert 'MPI_REAL16' in text
    # KIND target still uses stock MPI_SUM
    assert 'MPI_SUM' in text
    assert 'MPI_DD_SUM' not in text


# ---------------------------------------------------------------------------
# Per-file BLACS substitution (COMPLEX clone)
# ---------------------------------------------------------------------------


_BLACS_COMPLEX_SAMPLE = '''\
#include "Bdef.h"
void BI_zvvsum(Int N, DCOMPLEX *v1, DCOMPLEX *v2) {
   N *= 2;
   double *r1 = (double*)v1;
   double *r2 = (double*)v2;
   for (int k = 0; k < N; k++) r1[k] += r2[k];
   MPI_Reduce(v1, v2, N, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);
}

void Czgebs2d_aux(DCOMPLEX *buf, Int n) {
   (void)buf;
}
'''


def test_blacs_complex_clone_multifloats(tmp_path):
    src = tmp_path / 'BI_zvvsum.c'
    src.write_text(_BLACS_COMPLEX_SAMPLE)
    target = load_target('multifloats')
    result = migrate_c_file_to_string(src, target)
    assert result is not None
    new_name, text = result
    assert new_name == f'BI_{_CP}vvsum.c'
    assert 'complex64x2 *v1' in text
    assert 'float64x2 *r1' in text
    assert 'MPI_COMPLEX64X2' in text
    # Complex files use the zz reduction op
    assert 'MPI_ZZ_SUM' in text
    assert 'MPI_SUM,' not in text
    assert f'C{_CP}gebs2d_aux' in text


# ---------------------------------------------------------------------------
# _patch_bdef_header
# ---------------------------------------------------------------------------


_BDEF_SKELETON = '''\
typedef struct {double r, i;} DCOMPLEX;
typedef struct {float r, i;} SCOMPLEX;

#define BI_zvmcopy(m, n, A, lda, buff) \\
        BI_dvmcopy(2*(m), (n), (double *) (A), 2*(lda), (double *) (buff))
#define zgamn2d_   zgamn2d
#define zgamn2d_   ZGAMN2D
'''


def test_patch_bdef_header_multifloats_includes_mfc_header(tmp_path):
    bdef = tmp_path / 'Bdef.h'
    bdef.write_text(_BDEF_SKELETON)
    mf = load_target('multifloats')
    template_vars = _build_sub_vars(mf)
    _patch_bdef_header(bdef, mf, template_vars)
    out = bdef.read_text()
    # Multifloats branch uses libmfc's plain-C header for type defs
    assert '#include "multifloats_bridge.h"' in out
    # Should NOT emit the legacy __float128 typedef
    assert '__float128' not in out
    assert 'typedef __float128' not in out
    # Forward decls for the BI_*mvcopy/vmcopy helpers still present
    assert f'BI_{_RP}mvcopy' in out
    assert f'BI_{_RP}vmcopy' in out
    # Name-mangling defines for the dd/zz prefixes
    assert f'{_RP}gamn2d_' in out
    assert f'{_CP}gamn2d_' in out


def test_patch_bdef_header_kind16_emits_legacy_typedef(tmp_path):
    """Regression: KIND=16 target keeps the original __float128 block."""
    bdef = tmp_path / 'Bdef.h'
    bdef.write_text(_BDEF_SKELETON)
    k16 = load_target('kind16')
    template_vars = _build_sub_vars(k16)
    _patch_bdef_header(bdef, k16, template_vars)
    out = bdef.read_text()
    assert 'typedef __float128 QREAL' in out
    assert '#include "multifloats_bridge.h"' not in out


# ---------------------------------------------------------------------------
# _apply_overrides
# ---------------------------------------------------------------------------


def test_apply_overrides_overwrites_clone(tmp_path):
    output_dir = tmp_path / 'output'
    output_dir.mkdir()
    src_dir = tmp_path / 'src'
    src_dir.mkdir()

    name = f'BI_{_RP}vvsum.c'
    # Pre-existing clone (e.g. produced by the regex pass) -- will be
    # overwritten.
    (output_dir / name).write_text('// clone produced by regex\n')

    # Hand-written override.
    override_text = f'// hand-written override\nvoid BI_{_RP}vvsum(...) {{}}\n'
    (src_dir / name).write_text(override_text)

    applied = _apply_overrides(
        output_dir,
        [(src_dir / name, name)],
    )
    assert applied == [name]
    assert (output_dir / name).read_text() == override_text


def test_apply_overrides_missing_source_raises(tmp_path):
    output_dir = tmp_path / 'output'
    output_dir.mkdir()
    with pytest.raises(FileNotFoundError):
        _apply_overrides(
            output_dir,
            [(tmp_path / 'does_not_exist.c', f'BI_{_RP}vvsum.c')],
        )


# ---------------------------------------------------------------------------
# _make_rename_substituter -- multi-char prefix expansion
# ---------------------------------------------------------------------------


def test_rename_substituter_equal_length_prefix():
    """KIND-style single-char swap: D->Q, lengths match."""
    rename_map = {'PB_CDTYPESET': 'PB_CQTYPESET', 'DGEMM': 'QGEMM'}
    pattern, combined = _build_rename_regex(rename_map)
    sub = _make_rename_substituter(pattern, combined)

    # Lowercase identifier preserves case
    assert pattern.sub(sub, 'PB_Cdtypeset') == 'PB_Cqtypeset'
    # Trailing underscore form
    assert pattern.sub(sub, 'dgemm_') == 'qgemm_'
    # Uppercase form
    assert pattern.sub(sub, 'DGEMM_') == 'QGEMM_'


def test_rename_substituter_multichar_expansion():
    """Multifloats-style expansion: D->DD, target one char longer.

    Regression for the zip(src, new_lower) bug that truncated the
    trailing 't' in PB_Cdtypeset -> PB_Cddtypese (missing 't').
    """
    rename_map = {'PB_CDTYPESET': 'PB_CDDTYPESET', 'DGEMM': 'DDGEMM'}
    pattern, combined = _build_rename_regex(rename_map)
    sub = _make_rename_substituter(pattern, combined)

    # Mixed-case PascalCase identifier: head case preserved, inserted
    # 'd' takes the case of the original precision letter (lowercase).
    assert pattern.sub(sub, 'PB_Cdtypeset') == 'PB_Cddtypeset'
    # Trailing underscore form
    assert pattern.sub(sub, 'dgemm_') == 'ddgemm_'
    # Uppercase form: inserted char is uppercase too
    assert pattern.sub(sub, 'DGEMM_') == 'DDGEMM_'


def test_rename_substituter_does_not_match_inside_longer_identifier():
    """\\b boundaries prevent DGER inside PDGER from being renamed."""
    rename_map = {'DGER': 'DDGER'}
    pattern, combined = _build_rename_regex(rename_map)
    sub = _make_rename_substituter(pattern, combined)

    assert pattern.sub(sub, 'pdger_') == 'pdger_'    # untouched
    assert pattern.sub(sub, 'dger_') == 'ddger_'


# ---------------------------------------------------------------------------
# _apply_aliases_to_original — copy-original alias / pointer-cast pass (UB-04)
# ---------------------------------------------------------------------------


def test_aliases_to_original_renames_local_decl():
    """Regression for UB-04: a precision-independent dispatcher like
    PB_Ctzher2k declares ``cmplx16 Calph16;`` as a local. The
    copy-originals path must apply ``c_type_aliases`` so the local
    declaration becomes ``cmplxQ Calph16;`` on the kind16 target —
    otherwise PB_Cconjg writes only 16 of 32 bytes into the buffer
    and the conjugated alpha is corrupted."""
    src = textwrap.dedent("""\
        void PB_Ctzher2k(...) {
            cmplx16 Calph16;
            int i;
        }
    """)
    template_vars = {'REAL_TYPE': 'quad', 'RPU': 'Q'}
    aliases = [{'from': ['cmplx', 'cmplx16'], 'to': 'cmplx{RPU}'}]
    out = _apply_aliases_to_original(src, template_vars, aliases, None)
    assert 'cmplxQ Calph16;' in out
    assert 'cmplx16' not in out


def test_aliases_to_original_pointer_cast_replacement():
    """Regression for UB-04: PB_Cconjg.c uses ``((double*)CALPHA)`` and
    ``((float*)CALPHA)`` to step through scalar buffers in 8-byte and
    4-byte strides. The copy-originals path must rewrite the casts so
    a kind16 cmplxQ buffer (32 bytes total, 16 bytes per real component)
    is stepped in 16-byte strides."""
    src = textwrap.dedent("""\
        case DCPLX:
            ((double*)(CALPHA))[REAL_PART] =  ((double*)(ALPHA))[REAL_PART];
            ((double*)(CALPHA))[IMAG_PART] = -((double*)(ALPHA))[IMAG_PART];
            break;
        case SCPLX:
            ((float*)(CALPHA))[REAL_PART]  =  ((float*)(ALPHA))[REAL_PART];
            ((float*)(CALPHA))[IMAG_PART]  = -((float*)(ALPHA))[IMAG_PART];
            break;
    """)
    template_vars = {'REAL_TYPE': 'quad', 'RPU': 'Q'}
    cast_aliases = [{'from': ['(double*)', '(float*)'], 'to': '({REAL_TYPE}*)'}]
    out = _apply_aliases_to_original(src, template_vars, None, cast_aliases)
    assert '(quad*)' in out
    assert '(double*)' not in out
    assert '(float*)' not in out
    # The bare ``DCPLX``/``SCPLX`` dispatch labels must NOT be touched —
    # they are enum constants, not types. The function deliberately
    # avoids the broad ``double``/``float`` substitution for this reason.
    assert 'case DCPLX:' in out
    assert 'case SCPLX:' in out


def test_aliases_to_original_does_not_promote_bare_double_keyword():
    """The copy-originals pass must NOT do the broad ``double`` →
    ``REAL_TYPE`` substitution. A precision-independent dispatcher may
    have a literal ``double`` type that is part of the dispatch logic
    (e.g. ``case DREAL: double_field = ...``); the broad sub would
    wrongly turn it into ``quad`` on kind16, breaking SCPLX/DCPLX
    discrimination."""
    src = 'double tmp; float other;\n'
    template_vars = {'REAL_TYPE': 'quad', 'RPU': 'Q'}
    out = _apply_aliases_to_original(src, template_vars, None, None)
    # No aliases configured → text unchanged.
    assert out == src


def test_aliases_to_original_no_aliases_is_identity():
    src = 'cmplx16 X; ((double*)P)[0] = 0.0;\n'
    out = _apply_aliases_to_original(src, {}, None, None)
    assert out == src
