"""MPI datatype / reduce-op rewriting for extended-precision reductions (Cluster L).

Rewrites the standard MPI datatype and MPI_SUM handles to the custom ops the
migrated archives need (MPI_QQ_* / MPI_DD_* / ...) and injects the matching
``USE ..._mpi_f`` clause. Extracted verbatim from ``fortran_migrator.py``.
"""
import re

from ..target_mode import TargetMode
from .lex import (
    _END_PROC_RE, _PROC_HEADER_RE, _count_open_parens, is_continuation_line,
)


# Fortran-side MPI datatype name rewriter. In an `s*`/`c*` source MUMPS
# uses ``MPI_REAL`` and ``MPI_COMPLEX``; in `d*`/`z*` it uses
# ``MPI_DOUBLE_PRECISION`` and ``MPI_DOUBLE_COMPLEX``. After migration
# both halves must refer to the target's wider datatype (e.g.
# ``MPI_REAL16``/``MPI_COMPLEX32`` for kind16). Without this pass the two
# halves' outputs disagree on every MPI call, even though they are
# semantically identical. Word boundaries keep the rewrite idempotent —
# ``MPI_REAL16`` already has no ``\b`` after ``REAL`` so it is not
# rematched.
_MPI_DOUBLE_COMPLEX_RE = re.compile(r'\bMPI_DOUBLE_COMPLEX\b')


_MPI_DOUBLE_PRECISION_RE = re.compile(r'\bMPI_DOUBLE_PRECISION\b')


_MPI_COMPLEX_RE = re.compile(r'\bMPI_COMPLEX\b')


_MPI_REAL_RE = re.compile(r'\bMPI_REAL\b')


# Fortran ``INCLUDE 'foo.h'`` whose filename starts with an arithmetic
# letter. MUMPS uses this for per-arithmetic C-interop headers
# (``INCLUDE 'dmumps_struc.h'`` in dmumps_struc_def.F). After migration
# the file on disk is renamed by target_filename's first-char fallback
# (dmumps_struc.h → qmumps_struc.h); the INCLUDE string must be
# rewritten to match, otherwise the compiler can't find it.
_INCLUDE_PREFIXED_H_RE = re.compile(
    r"(INCLUDE\s*['\"])([SsDdCcZz])([\w]*\.h)(['\"])",
    re.IGNORECASE,
)


def _rewrite_prefixed_includes(source: str, target_mode: TargetMode) -> str:
    from ..prefix_classifier import CHAR_TYPE
    pmap = target_mode.prefix_map

    def sub(m: re.Match) -> str:
        head, prefix, rest, tail = m.group(1), m.group(2), m.group(3), m.group(4)
        family = CHAR_TYPE.get(prefix.upper())
        new = pmap.get(family) if family else None
        if not new:
            return m.group(0)
        first = new if prefix.isupper() else new.lower()
        return head + first + rest + tail

    return _INCLUDE_PREFIXED_H_RE.sub(sub, source)


def _rewrite_mpi_datatypes(source: str, target_mode: TargetMode) -> str:
    mpi_real = target_mode.c_mpi_real
    mpi_complex = target_mode.c_mpi_complex
    if not mpi_real and not mpi_complex:
        return source
    if mpi_complex:
        source = _MPI_DOUBLE_COMPLEX_RE.sub(mpi_complex, source)
        source = _MPI_COMPLEX_RE.sub(mpi_complex, source)
    if mpi_real:
        source = _MPI_DOUBLE_PRECISION_RE.sub(mpi_real, source)
        source = _MPI_REAL_RE.sub(mpi_real, source)
    return source


# Token-context MPI reduction-op rewriter for multifloats. The target
# ships user-defined ops (MPI_DD_SUM / MPI_DD_AMX / MPI_DD_AMN, plus
# the ZZ-prefixed complex variants) that operate on the user-defined
# float64x2 / complex64x2 datatypes; the stock MPI_SUM / MPI_MAX /
# MPI_MIN are undefined for those datatypes and the runtime aborts.
# The C migrator handles this by file-arithmetic context (real /
# complex sibling files use separate rule lists) but Fortran sources
# from MUMPS mix integer reductions (MPI_INTEGER + MPI_SUM, which
# MUST stay MPI_SUM) with floating-point reductions in the same
# translation unit, so a blunt file-prefix rewrite would corrupt the
# integer calls. Instead we match each MPI_(I?All)Reduce[_scatter]
# call as a whole, run *after* _rewrite_mpi_datatypes, and only
# rewrite the op inside a call whose argument list contains the
# rewritten floating datatype token (e.g. MPI_FLOAT64X2). Integer-
# typed calls are left alone by construction. For KIND targets the
# c_mpi_(sum|max|min)_* fields expand to 'MPI_(SUM|MAX|MIN)' or are
# unset, and the pass is a no-op.
_MPI_REDUCE_CALL_RE = re.compile(
    r'\bMPI_(?:I?ALL)?REDUCE(?:_SCATTER)?\b\s*'
    r'\((?:[^()]|\([^()]*\))*\)',
    re.IGNORECASE,
)


_MPI_OP_RE = {
    'MPI_SUM': re.compile(r'\bMPI_SUM\b'),
    'MPI_MAX': re.compile(r'\bMPI_MAX\b'),
    'MPI_MIN': re.compile(r'\bMPI_MIN\b'),
}


def _rewrite_mpi_sum(source: str, target_mode: TargetMode) -> str:
    real_tok = target_mode.c_mpi_real
    complex_tok = target_mode.c_mpi_complex
    rules: list[tuple[re.Pattern, str | None, str | None]] = []
    for stock, mf_real, mf_complex in (
        ('MPI_SUM', target_mode.c_mpi_sum_real, target_mode.c_mpi_sum_complex),
        ('MPI_MAX', target_mode.c_mpi_max_real, target_mode.c_mpi_max_complex),
        ('MPI_MIN', target_mode.c_mpi_min_real, target_mode.c_mpi_min_complex),
    ):
        r = mf_real if mf_real and mf_real != stock else None
        c = mf_complex if mf_complex and mf_complex != stock else None
        if r or c:
            rules.append((_MPI_OP_RE[stock], r, c))
    if not rules:
        return source

    def rewrite(m: re.Match) -> str:
        call = m.group(0)
        is_complex = bool(complex_tok) and complex_tok in call
        is_real = bool(real_tok) and real_tok in call
        if not (is_complex or is_real):
            return call
        for pat, r, c in rules:
            if is_complex and c:
                call = pat.sub(c, call)
            elif is_real and r:
                call = pat.sub(r, call)
        return call

    return _MPI_REDUCE_CALL_RE.sub(rewrite, source)


# Predefined MPI handle names that mpif.h / the MPI module already
# declares. A target may keep one of these as its datatype (kind16
# transports quad reals as the standard ``MPI_REAL16``) while routing
# only the *reduction ops* through custom handles — those standard
# names must NOT trigger a support-module ``USE`` (they are already in
# scope) and are never declared by the support module.
_STANDARD_MPI_HANDLES = frozenset({
    'MPI_SUM', 'MPI_MAX', 'MPI_MIN',
    'MPI_REAL16', 'MPI_COMPLEX32',
    'MPI_LONG_DOUBLE',
    'MPI_DOUBLE_PRECISION', 'MPI_DOUBLE_COMPLEX',
    'MPI_REAL', 'MPI_COMPLEX', 'MPI_REAL8', 'MPI_COMPLEX16',
})


def _custom_mpi_tokens(target_mode: TargetMode) -> tuple[str, ...]:
    """The custom (non-predefined) MPI handle names the target's support
    module declares — i.e. the ``c_mpi_*`` datatype/op values that are
    not standard MPI handles. For multifloats these are the custom
    datatypes (``MPI_FLOAT64X2`` ...) plus all six reduction ops; for
    kind16 the datatypes stay standard so only the six quad ops
    (``MPI_QQ_SUM`` ...) qualify."""
    seen: dict[str, None] = {}
    for v in (target_mode.c_mpi_real, target_mode.c_mpi_complex,
              target_mode.c_mpi_sum_real, target_mode.c_mpi_sum_complex,
              target_mode.c_mpi_max_real, target_mode.c_mpi_max_complex,
              target_mode.c_mpi_min_real, target_mode.c_mpi_min_complex):
        if v and v not in _STANDARD_MPI_HANDLES:
            seen[v] = None
    return tuple(seen)


def insert_use_multifloats_mpi_f(source: str, target_mode: TargetMode) -> str:
    """Inject ``USE <c_mpi_module>`` after each procedure header in
    sources that reference the target's custom MPI handle names. Runs
    after ``_rewrite_mpi_datatypes`` / ``_rewrite_mpi_sum`` so the
    rewritten tokens are visible. Migrated MUMPS Fortran calls MPI
    directly with names like ``MPI_FLOAT64X2`` (multifloats) or
    ``MPI_QQ_SUM`` (kind16); without a Fortran module declaring them,
    gfortran rejects every site with "has no IMPLICIT type". The
    support module (``multifloats_mpi_f`` / ``quad_mpi_f``) bind-to-C
    declares them as default-INTEGER handles populated at MPI init
    time. It exports only ``MPI_*``-prefixed names plus its ``*_init``
    routine, so no ONLY clause is needed — collision risk is nil and
    the bare USE keeps fixed-form lines under column 72.

    Targets without a ``mpi_module`` (genuine d/z, kind10) declare no
    custom handles and are a no-op."""
    module = target_mode.c_mpi_module
    if not module:
        return source
    tokens = _custom_mpi_tokens(target_mode)
    if not tokens:
        return source
    any_token_re = re.compile(r'\b(?:' + '|'.join(map(re.escape, tokens)) + r')\b')
    if not any_token_re.search(source):
        return source

    lines = source.splitlines(keepends=True)
    proc_header_re = _PROC_HEADER_RE
    end_proc_re = _END_PROC_RE

    result: list[str] = []
    i = 0
    while i < len(lines):
        line = lines[i]
        result.append(line)
        m = proc_header_re.match(line)
        if not m:
            i += 1
            continue
        # Walk past the procedure header (continuations + CPP blocks),
        # mirroring insert_use_multifloats so paren-balanced argument
        # lists that span ``#if defined(metis)`` blocks are handled.
        j = i + 1
        paren_depth = _count_open_parens(line)
        prev_has_amp = line.rstrip().rstrip('\n').endswith('&')
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
        # Find this procedure's matching END via depth tracking. Many
        # MUMPS files put an INTERFACE block inside the outer SUBROUTINE
        # whose internal SUBROUTINE signatures look syntactically like
        # procedure headers; their END SUBROUTINE would falsely
        # terminate a naive linear scan, missing the MPI_FLOAT64X2
        # sites further down in the outer body. Counting nested opens
        # vs closes lets us land on the actual matching END at depth 0.
        depth = 1
        k = j
        while k < len(lines):
            if proc_header_re.match(lines[k]):
                depth += 1
            elif end_proc_re.match(lines[k]):
                depth -= 1
                if depth == 0:
                    break
            k += 1
        body_text = ''.join(lines[j:k + 1])
        if any_token_re.search(body_text):
            use_stmt = f"USE {module}"
            already_has = any(use_stmt in lines[kk]
                              for kk in range(j, min(j + 30, len(lines))))
            if not already_has:
                indent = m.group(1)
                ind = indent if indent.strip() else '      '
                result.append(f"{ind}{use_stmt}\n")
        i = j
    return ''.join(result)
