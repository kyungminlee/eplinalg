"""LA_CONSTANTS USE rewriting + kind-parameter helpers (Cluster J tail).

Rewrites ``USE LA_CONSTANTS`` clauses and KIND parameter declarations to the
migrated names/values. Extracted verbatim from ``fortran_migrator.py``.
"""
import re

from ..target_mode import TargetMode
from .renames import replace_routine_names


_KIND_PARAM_NAMES = r'(?:wp|sp|dp)'


_KIND_PARAM_RE = re.compile(rf'(integer\s*,\s*parameter\s*::\s*{_KIND_PARAM_NAMES}\s*=\s*)(?:kind\s*\(\s*1\.[de]0\s*\)|real(?:32|64|128))', re.IGNORECASE)


def _replace_kind_parameter(line: str, target_mode: TargetMode) -> str:
    if target_mode.is_kind_based: return _KIND_PARAM_RE.sub(rf'\g<1>{target_mode.kind_suffix}', line)
    return ('! ' + line) if _KIND_PARAM_RE.search(line) else line


_ISO_USE_ONLY_RE = re.compile(r'^(?P<lead>\s*)USE\s*,\s*INTRINSIC\s*::\s*ISO_FORTRAN_ENV\s*,\s*ONLY\s*:\s*(?P<names>[^\n!]*?)\s*(?P<tail>!.*)?$', re.IGNORECASE)


def _strip_iso_fortran_env_realN(line: str) -> str:
    m = _ISO_USE_ONLY_RE.match(line)
    if not m: return line
    names = [n.strip() for n in m.group('names').split(',') if n.strip()]
    kept = [n for n in names if not re.fullmatch(r'real(?:32|64|128)', n, re.IGNORECASE)]
    if not kept: return ''
    tail = (' ' + m.group('tail')) if m.group('tail') else ''
    return f'{m.group("lead")}USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: {", ".join(kept)}{tail}'


def rewrite_la_constants_use(source: str, target_mode: TargetMode) -> str:
    """Rewrite ``USE LA_CONSTANTS`` clauses for the chosen target.

    KIND mode (extended precision): the LAPACK la_constants module is
    cloned to a per-target module (``LA_CONSTANTS_EY`` for KIND=10,
    ``LA_CONSTANTS_QX`` for KIND=16) named via ``la_constants_suffix``,
    so we just rename the module reference and rename each constant to
    its prefixed equivalent (E*/Y* or Q*/X*).

    Multifloats mode: there is no ``la_constants_mw`` module — instead,
    we rewrite the import to point at the real ``multifloats`` module
    and rename each LAPACK constant (``dzero``, ``dsafmin``, ...) to its
    multifloats equivalent (``MF_ZERO``, ``MF_SAFMIN``, ...). The
    ``wp=>dp`` rename entry is dropped because ``wp`` is no longer
    meaningful once the type becomes ``TYPE(float64x2)``.
    """
    const_renames = _la_constants_rename_map(target_mode)
    lines, result, in_use_stmt = source.split('\n'), [], False
    suffix = target_mode.la_constants_suffix
    target_module_upper = f'LA_CONSTANTS{suffix}'
    target_module_lower = f'la_constants{suffix.lower()}'
    target_xisnan_upper = f'LA_XISNAN{suffix}'
    target_xisnan_lower = f'la_xisnan{suffix.lower()}'

    for line in lines:
        upper = line.upper().lstrip()
        if re.search(r'\bUSE\s+LA_XISNAN\b', upper) and target_xisnan_upper not in upper:
            line = re.sub(
                r'(?i)\bLA_XISNAN\b',
                lambda m: target_xisnan_lower if m.group().islower() else target_xisnan_upper,
                line,
            )
        if re.search(r'\bUSE\s+LA_CONSTANTS\b', upper) and target_module_upper not in upper:
            in_use_stmt = True
            line = re.sub(
                r'(?i)\bLA_CONSTANTS\b',
                lambda m: target_module_lower if m.group().islower() else target_module_upper,
                line,
            )
        if in_use_stmt:
            line = replace_routine_names(line, const_renames)
            if target_mode.is_kind_based:
                # Rename ``wp=>dp`` / ``wp=>sp`` to ``wp=>qp`` (kind16) or
                # ``wp=>ep`` (kind10).  The target kind parameter is the
                # real prefix lowercased + "p".
                target_kp = target_mode.prefix_map['R'].lower() + 'p'
                for kindname in ('dp', 'sp'):
                    line = re.sub(
                        rf'(?i)\b{kindname}\b',
                        lambda m, kp=target_kp: kp if m.group().islower() else kp.upper(),
                        line,
                    )
            else:
                # Strip ``wp=>dp`` (D-source) and ``wp=>sp`` (S-source)
                # entries — both become meaningless after the migrator
                # collapses both halves to float64x2.
                for kindname in ('dp', 'sp'):
                    line = re.sub(rf',\s*wp\s*=>\s*{kindname}\s*,', ',', line, flags=re.IGNORECASE)
                    line = re.sub(rf',\s*wp\s*=>\s*{kindname}\s*(?=[!&]|$)', '', line, flags=re.IGNORECASE)
                    line = re.sub(rf'(ONLY\s*:\s*)wp\s*=>\s*{kindname}\s*,', r'\1', line, flags=re.IGNORECASE)
                    line = re.sub(rf'(ONLY\s*:\s*)wp\s*=>\s*{kindname}\s*(?=[!&]|$)', r'\1', line, flags=re.IGNORECASE)
            if not line.rstrip().endswith('&'):
                in_use_stmt = False
        result.append(line)
    return '\n'.join(result)


_LA_CONSTANTS_REAL_NAMES = (
    'ZERO', 'HALF', 'ONE', 'TWO', 'THREE', 'FOUR', 'EIGHT', 'TEN',
    'PREFIX', 'ULP', 'EPS',
    'SAFMIN', 'SAFMAX', 'SMLNUM', 'BIGNUM',
    'RTMIN', 'RTMAX',
    'TSML', 'TBIG', 'SSML', 'SBIG',
)


_LA_CONSTANTS_COMPLEX_NAMES = ('ZERO', 'HALF', 'ONE', 'PREFIX')


def _la_constants_rename_map(target_mode: TargetMode) -> dict[str, str]:
    """Build a rename map for the RHS of LA_CONSTANTS USE-clause aliases.

    Maps the LAPACK la_constants names ``DZERO``, ``DSAFMIN``, ``ZZERO``,
    etc. to the equivalent names exported by the target la_constants
    auxiliary module:

      KIND=10  → ``ezero``, ``ysafmin`` (la_constants_ey)
      KIND=16  → ``qzero``, ``xsafmin`` (la_constants_qx)
      multifloats → ``ddzero``, ``zzsafmin`` (la_constants_mw)

    Only the prefixed names are mapped — the LHS aliases ``zero``,
    ``half`` etc. are intentionally left untouched so the body of the
    routine continues to reference them through the local alias.
    """
    # Build S/D → target_real_prefix, C/Z → target_complex_prefix map
    # from the target's prefix_map (which maps R→prefix, C→prefix).
    real_pfx = target_mode.prefix_map.get('R', 'Q')
    cmplx_pfx = target_mode.prefix_map.get('C', 'X')
    pmap = {'S': real_pfx, 'D': real_pfx, 'C': cmplx_pfx, 'Z': cmplx_pfx}

    renames: dict[str, str] = {}
    for p in ('S', 'D'):
        for base in _LA_CONSTANTS_REAL_NAMES:
            renames[p + base] = pmap[p] + base
    for p in ('C', 'Z'):
        for base in _LA_CONSTANTS_COMPLEX_NAMES:
            renames[p + base] = pmap[p] + base
    return renames
