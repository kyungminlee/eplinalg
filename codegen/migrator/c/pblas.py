"""PBLAS-specific C migration knowledge: header typedef patching and
the multifloats struct-type idiom rewrites for PBLAS entry points.
"""

import re
from pathlib import Path


def patch_pblas_header(pblas_path: Path,
                       template_vars: dict[str, str]) -> None:
    """Insert extended-precision typedefs into pblas.h for KIND targets.

    The migrator replaces ``double``/``complex16`` with the target type
    names (e.g. QREAL, QCOMPLEX) in .c files, but pblas.h keeps the
    original types.  We add typedefs so both old and new names compile.
    """
    real_type = template_vars['REAL_TYPE']
    complex_type = template_vars['COMPLEX_TYPE']
    c_type = template_vars['C_REAL_TYPE']

    block = (f'typedef {c_type} {real_type};\n'
             f'typedef struct {{ {real_type} re, im; }} {complex_type};\n')

    text = pblas_path.read_text(errors='replace')
    if real_type in text and 'typedef' in text.split(real_type)[0]:
        return  # already patched
    # Insert just before the first "typedef struct" line
    marker = 'typedef struct'
    idx = text.find(marker)
    if idx >= 0:
        text = text[:idx] + block + '\n' + text[idx:]
        pblas_path.write_text(text)


# Cost-estimate local variable names used by PBLAS Level-3 entry points
# (pdgemm/pdsymm/pdsyrk/...) for algorithm selection. These are pure
# heuristic doubles that must NOT be promoted to the multifloats struct
# type, otherwise (double) casts and `*=` arithmetic in the cost model
# stop compiling. Survey of pblas/SRC/p[dz]*.c shows just two declaration
# lines containing these names; recognising them by name is sufficient.
#
# Note: ``tmp\d+`` deliberately excluded — PBLAS pairs them with
# ABest/ACest/BCest on the same line (so the line still matches via
# those anchors), while XBLAS uses bare ``double tmp1;`` accumulators
# that DO need to be promoted to float64x2 for working-precision math.
PBLAS_COST_LOCAL = (
    r'ABest|ACest|BCest|ABestL|ABestR|Best'
)


_MF_PBLAS_PART = r'(\w+)\s*\[\s*(REAL_PART|IMAG_PART)\s*\]'


def apply_multifloats_pblas_subs(text: str) -> str:
    """Rewrite PBLAS scalar quick-return checks for the float64x2_t
    struct type. C operators ``==`` / ``!=`` are not defined on
    structs, so we replace ``ALPHA[REAL_PART] == ZERO`` with the
    inline macro ``MF_IS_ZERO(ALPHA[REAL_PART])`` (and likewise for
    ``ONE`` / ``IMAG_PART`` / negation).
    """
    # ZERO comparisons
    text = re.sub(
        _MF_PBLAS_PART + r'\s*==\s*ZERO\b',
        r'MF_IS_ZERO(\1[\2])', text)
    text = re.sub(
        _MF_PBLAS_PART + r'\s*!=\s*ZERO\b',
        r'(!MF_IS_ZERO(\1[\2]))', text)
    # ONE comparisons
    text = re.sub(
        _MF_PBLAS_PART + r'\s*==\s*ONE\b',
        r'MF_IS_ONE(\1[\2])', text)
    text = re.sub(
        _MF_PBLAS_PART + r'\s*!=\s*ONE\b',
        r'(!MF_IS_ONE(\1[\2]))', text)
    # Integer truncation. PBLAS packs integer indices into work
    # buffers: ``(Int)(work[1])`` or ``(Int)(work[1][REAL_PART])``
    # becomes invalid on struct types. Rewrite to ``mf_to_int(...)``
    # (defined in both multifloats_c.h and the C++ bridge header).
    text = re.sub(
        r'\(Int\)\(\s*(\w+(?:\[\s*\w+\s*\])+)\s*\)',
        r'mf_to_int(\1)', text)
    return text
