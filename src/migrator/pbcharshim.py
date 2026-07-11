"""Hidden-CHARACTER-length ABI bridge for migrated PBLAS typesets.

The type-generic PBLAS drivers reach the BLAS through ``PBTYP_T`` function
pointers (``GEMM_T``, ``HEMV_T``, ...) whose C prototypes use

    typedef char * F_CHAR_T;   #define C2F_CHAR(a) (a)

i.e. they pass NO hidden Fortran CHARACTER-length arguments and reserve no
stack for them.  A migrated gfortran facade (``ygemm_``, ``ytrsm_``, ...) is a
genuine Fortran routine: per the SysV AMD64 ABI it may freely use its incoming
hidden-length slots as scratch, and gfortran's codegen for these routines does
exactly that.  Because the C caller never provisioned those slots, the writes
clobber the caller's live locals -> wild pointers -> SIGSEGV / garbage MPI
ranks.  Reference/MKL BLAS merely never happen to touch those slots, which is
why the s/c/d/z path is unaffected.

The fix is a per-(field, leaf) C trampoline with the driver's exact no-length
signature that forwards to the facade with the correct number of hidden
CHARACTER lengths (all 1 -- every BLAS char argument is ``CHARACTER*1``)
appended.  The facade's scratch stores then land in the trampoline's own
outgoing-argument area instead of the caller's frame.

This module emits ``PBcharshim.h`` (the reusable substrate: length-carrying
function-pointer typedefs + per-shape generator macros) next to the migrated
typesets and rewrites each extended typeset so every ``F<slot>`` char-callback
is assigned a trampoline instead of the raw leaf symbol.  Confined to the
extended (e/y/q/x, and by hand m/w) typesets; the MKL s/c/d/z path is
untouched.
"""
from __future__ import annotations

import re
from pathlib import Path

# The verified header substrate, shipped alongside this module.
_PBCHARSHIM_H_PATH = Path(__file__).with_name('PBcharshim.h')


def pbcharshim_header_text() -> str:
    """Return the exact bytes of the reusable ``PBcharshim.h`` substrate."""
    return _PBCHARSHIM_H_PATH.read_text()


# The 27 ``PBTYP_T`` fields whose callee takes one or more Fortran CHARACTER
# arguments, mapped to (generator-macro suffix, PBTYP_T field cast typedef).
# Every other field (Fmmadd, Fvvdotc, Fgerc, Fgeru, Faxpy, ...) has no
# CHARACTER dummy and must NOT be wrapped.  Shapes are shared across
# callbacks with the same layout (verified against pblas.h): GEMV/AGEMV,
# SYMV/ASYMV/HEMV/AHEMV, TRMV/TRSV, SYR/HER, SYR2/HER2, SYMM/HEMM,
# SYRK/HERK, SYR2K/HER2K, TRMM/TRSM, and the three TZSCAL aliases.
FIELD_SHIM: dict[str, tuple[str, str]] = {
    'Ftzpad':    ('TZPAD',    'TZPAD_T'),
    'Ftzpadcpy': ('TZPADCPY', 'TZPADCPY_T'),
    'Ftzscal':   ('TZSCAL',   'TZSCAL_T'),
    'Fhescal':   ('TZSCAL',   'TZSCAL_T'),
    'Ftzcnjg':   ('TZSCAL',   'TZSCAL_T'),
    'Fgemv':     ('GEMV',     'GEMV_T'),
    'Fagemv':    ('GEMV',     'AGEMV_T'),
    'Fsymv':     ('SYMV',     'SYMV_T'),
    'Fhemv':     ('SYMV',     'HEMV_T'),
    'Fasymv':    ('SYMV',     'ASYMV_T'),
    'Fahemv':    ('SYMV',     'AHEMV_T'),
    'Ftrmv':     ('TRMV',     'TRMV_T'),
    'Ftrsv':     ('TRMV',     'TRSV_T'),
    'Fatrmv':    ('ATRMV',    'ATRMV_T'),
    'Fsyr':      ('SYR',      'SYR_T'),
    'Fher':      ('SYR',      'HER_T'),
    'Fsyr2':     ('SYR2',     'SYR2_T'),
    'Fher2':     ('SYR2',     'HER2_T'),
    'Fgemm':     ('GEMM',     'GEMM_T'),
    'Fsymm':     ('SYMM',     'SYMM_T'),
    'Fhemm':     ('SYMM',     'HEMM_T'),
    'Fsyrk':     ('SYRK',     'SYRK_T'),
    'Fherk':     ('SYRK',     'HERK_T'),
    'Fsyr2k':    ('SYR2K',    'SYR2K_T'),
    'Fher2k':    ('SYR2K',    'HER2K_T'),
    'Ftrmm':     ('TRMM',     'TRMM_T'),
    'Ftrsm':     ('TRMM',     'TRSM_T'),
}

# ``TypeStruct.<Field> = <leaf>_ ;`` — the raw callback assignment produced by
# the stock typeset (after rename + type substitution).  RHS is a bare Fortran
# leaf symbol (``ygemm_``, ``etzscal_``, ...); aliased real slots reuse a leaf
# (``Fhescal = etzscal_``), which is fine — each gets its own trampoline.
_ASSIGN_RE = re.compile(
    r'(?P<pre>TypeStruct\.)(?P<field>F\w+)(?P<mid>\s*=\s*)'
    r'(?P<rhs>[A-Za-z_]\w*)(?P<post>\s*;)'
)

# Insert the shim block after the final top-of-file ``#include "..."`` line.
_INCLUDE_RE = re.compile(r'^#include\s*"[^"]+"\s*$', re.MULTILINE)


def is_typeset_stem(stem: str) -> bool:
    """True for ``pb_c<letter>typeset`` clone stems (PBLAS only)."""
    return bool(re.fullmatch(r'pb_c[a-z]+typeset', stem))


def transform_typeset(text: str) -> tuple[str, bool]:
    """Rewrite a migrated PBLAS typeset to route char-callbacks through
    hidden-length trampolines.

    Returns ``(new_text, changed)``.  Idempotent: a text already carrying the
    ``PBcharshim.h`` include is returned unchanged.
    """
    if 'PBcharshim.h' in text:
        return text, False

    shims: list[str] = []  # (macro, tag, rhs) rendered EP_MK_* lines
    seen_tags: set[str] = set()

    def _repl(m: re.Match) -> str:
        field = m.group('field')
        info = FIELD_SHIM.get(field)
        if info is None:
            return m.group(0)  # non-char field: leave untouched
        macro, cast = info
        rhs = m.group('rhs')
        tag = field[1:].lower()  # Fgemm -> gemm, Fsyr2k -> syr2k
        if tag not in seen_tags:
            seen_tags.add(tag)
            shims.append(f'EP_MK_{macro}({tag}, {rhs})')
        return (f"{m.group('pre')}{field}{m.group('mid')}"
                f"({cast}) EP_SHIM({tag}){m.group('post')}")

    new_text = _ASSIGN_RE.sub(_repl, text)
    if not shims:
        return text, False

    block = (
        '\n/*\n'
        '*  Hidden-CHARACTER-length ABI trampolines.  Each F<slot> callback\n'
        '*  below is reached by the type-generic drivers through a C prototype\n'
        '*  that passes NO Fortran CHARACTER lengths; the migrated gfortran\n'
        '*  leaf expects them and spills scratch into the (unprovisioned)\n'
        '*  hidden-length slots, clobbering the caller frame.  These\n'
        '*  trampolines carry the driver\'s no-length signature inward and\n'
        '*  append the lengths (all 1) outward.  See PBcharshim.h.\n'
        '*/\n'
        '#include "PBcharshim.h"\n\n'
        + '\n'.join(shims) + '\n'
    )

    includes = list(_INCLUDE_RE.finditer(new_text))
    if includes:
        pos = includes[-1].end()
        new_text = new_text[:pos] + '\n' + block + new_text[pos:]
    else:
        new_text = block + new_text

    return new_text, True
