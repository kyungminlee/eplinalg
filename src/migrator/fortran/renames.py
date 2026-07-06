"""Routine / include / xerbla / known-constant renaming (Clusters D+E).

Applies the recipe's symbol renames to routine names, INCLUDE filenames,
XERBLA string arguments, and known named constants. Extracted verbatim from
``fortran_migrator.py``.
"""
import re
import functools
from pathlib import Path

from ..target_mode import TargetMode
from .lex import is_comment_line


@functools.cache
def _known_constants_pattern(keys: frozenset) -> re.Pattern:
    """Cache the comma-name alternation in ``replace_known_constants``.

    The ``renames`` dict is invariant for the duration of a file
    (built once before the per-line loop), so the compiled pattern
    can be cached. The function is called per source line; rebuilding
    the alternation on every call was the dominant cost.
    """
    names_alt = '|'.join(re.escape(n) for n in sorted(keys, key=len, reverse=True))
    return re.compile(rf'(?<![A-Za-z0-9_])({names_alt})(?![A-Za-z0-9_])', re.IGNORECASE)


_ROUTINE_NAME_TOK_RE = re.compile(r'[A-Za-z_]\w*')


def replace_routine_names(line: str, rename_map: dict[str, str]) -> str:
    """Replace routine names using the rename map (case-preserving).

    Tokenize-then-lookup: scan ``line`` for identifier-shaped tokens
    and dict-look each one up. This sidesteps the multi-thousand
    alternation regex that Python's ``re`` engine evaluates by
    backtracking through every alternative at every position, which
    profiling showed dominated migration runtime (~60% of total time
    for LAPACK kind16). One small DFA-friendly tokenizer plus O(1)
    dict lookups is dramatically faster (70× on a representative
    LAPACK rename_map) and semantically equivalent: identifier tokens
    captured by ``[A-Za-z_]\\w*`` match exactly the ``\\b...\\b``-
    bounded forms the alternation regex used to find.
    """
    _, upper_map = _get_rename_pattern(rename_map)
    if not upper_map:
        return line

    def case_replace(m):
        tok = m.group(0)
        new = upper_map.get(tok.upper())
        if new is None:
            return tok
        if tok.isupper():
            return new.upper()
        if tok.islower():
            return new.lower()
        return new.capitalize()

    return _ROUTINE_NAME_TOK_RE.sub(case_replace, line)


_INCLUDE_RE = re.compile(
    r'''^(?P<lead>\s*(?:INCLUDE|include|Include))\s+(?P<q>['"])(?P<name>[^'"]+)(?P=q)(?P<tail>.*)$''',
)


def replace_include_filenames(line: str, rename_map: dict[str, str]) -> str:
    """Rewrite ``INCLUDE 'xxx.h'`` when ``xxx`` is in the rename_map.

    Needed for precision-prefixed headers (e.g. MUMPS's ``dmumps_struc.h``,
    which the ``mumps_struc`` recipe migrates to ``ddmumps_struc.h`` in
    multifloats mode). The including ``.F`` file's literal filename
    string isn't touched by :func:`replace_routine_names` because it's
    inside a quoted Fortran string, so it needs a dedicated rewrite.
    """
    # Fast reject: ``INCLUDE`` statements are rare (≈1 per file) so we
    # avoid the per-call rename_map scan until we know the line matches.
    m = _INCLUDE_RE.match(line)
    if not m:
        return line
    name = m.group('name')
    stem = Path(name).stem
    ext = Path(name).suffix
    upper_stem = stem.upper()
    # rename_map keys are already uppercase (built by build_rename_map);
    # fall back to a case-insensitive lookup if not.
    new_stem = rename_map.get(upper_stem) or rename_map.get(stem)
    if new_stem is None:
        return line
    new_name = (new_stem.upper() if stem.isupper() else
                (new_stem.lower() if stem.islower() else new_stem)) + ext
    return f'{m.group("lead")} {m.group("q")}{new_name}{m.group("q")}{m.group("tail")}'


# key: id(rename_map) -> (pattern, upper_map, source_dict). The source
# dict is stored in the value so it stays alive for as long as its cache
# entry does — this is what makes keying on id() safe (see below).
_RENAME_PATTERN_CACHE: dict[int, tuple[re.Pattern, dict[str, str], dict[str, str]]] = {}


def _get_rename_pattern(rename_map: dict[str, str]) -> tuple[re.Pattern, dict[str, str]]:
    # Cache key is id(rename_map). id() alone is NOT a safe key across
    # different dicts: CPython reuses the address of a garbage-collected
    # object, so a stale entry could be returned for an unrelated dict
    # that happens to land at the same address. This is not hypothetical
    # here — callers pass short-lived *per-file* ``file_rename_map``
    # objects (built fresh in ``_migrate_with_flang`` for every source
    # file), not just the one stable run-wide map. When such a dict was
    # freed and its id reused for the next file's map (a same-length map
    # collided under the old ``(id, len)`` key), the cache silently
    # handed back the previous file's ``upper_map`` and the current
    # file's own routine name was dropped from the rename — e.g. a file
    # renamed to ``qgemv.f`` kept ``SUBROUTINE DGEMV`` in its body,
    # exporting ``dgemv_`` and breaking the link. Two guards:
    #   1. The value keeps a reference to the source dict, so a live
    #      cache entry pins its id and no other live dict can reuse it.
    #   2. On a hit we still verify identity (``is``); a mismatch forces
    #      a recompute. Belt-and-suspenders against any future eviction.
    key = id(rename_map)
    cached = _RENAME_PATTERN_CACHE.get(key)
    if cached is not None and cached[2] is rename_map:
        return cached[0], cached[1]
    upper_map = {k.upper(): v for k, v in rename_map.items()}
    names = sorted(upper_map.keys(), key=len, reverse=True)
    pattern = re.compile(r'\b(' + '|'.join(re.escape(n) for n in names) + r')\b', re.IGNORECASE) if names else re.compile(r'(?!x)x')
    # Bound memory: per-file maps accumulate one entry each over a full
    # stage. Clear wholesale past a generous cap (entries are cheap to
    # rebuild; the hot path is the stable run-wide map, re-cached on the
    # next call).
    if len(_RENAME_PATTERN_CACHE) > 4096:
        _RENAME_PATTERN_CACHE.clear()
    _RENAME_PATTERN_CACHE[key] = (pattern, upper_map, rename_map)
    return pattern, upper_map


def _all_known_constant_renames(target_mode: TargetMode) -> dict[str, str]:
    """Combined uppercase rename map for known + la_constants names."""
    out: dict[str, str] = {}
    for k, v in target_mode.known_constants.items():
        out[k.upper()] = v
    for k, v in target_mode.la_constants_map.items():
        out[k.upper()] = v
    return out


_DECL_KEYWORD_RE = re.compile(
    r'^\s*(?:TYPE\s*\(|INTEGER|REAL|COMPLEX|LOGICAL|CHARACTER|'
    r'DOUBLE\s+PRECISION|DOUBLE\s+COMPLEX|PARAMETER\b|DATA\b|'
    r'IMPLICIT\b|DIMENSION\b|EXTERNAL\b|INTRINSIC\b|SAVE\b|'
    r'COMMON\b|EQUIVALENCE\b|USE\b|'
    # Keep-kind sentinels — the pre-pass masks ``DOUBLE PRECISION`` /
    # ``dble`` / ``dcmplx`` on keep-kind lines, and these decl-keyword
    # detection points must treat the masked tokens the same as the
    # un-masked originals so that subsequent passes (e.g.
    # ``replace_known_constants``) don't try to rewrite identifiers in
    # what is logically still a declaration line.
    r'__KEEPKIND_DP__\b)',
    re.IGNORECASE,
)


def replace_known_constants(
    line: str,
    target_mode: TargetMode,
    renames: dict[str, str] | None = None,
) -> str:
    """Substitute names in ``renames`` in the code portion of ``line`` only.

    ``renames`` is the per-file set of names that were removed from a
    local declaration by :func:`strip_known_constants_from_decls`. They
    are exactly the names that need module-import substitution. Names
    bound via a ``USE LA_CONSTANTS_MW, ONLY: zero=>dzero`` alias are
    NOT in this set and therefore are not touched here.

    When ``renames`` is ``None`` the function falls back to the union
    of ``target_mode.known_constants`` and ``target_mode.la_constants_map``;
    that legacy mode is preserved for callers that haven't yet been
    threaded through the new per-file API.

    Skips: comment lines (fixed-form column-1 ``C/c/*/!`` and bare
    ``!``-prefixed free-form), inline comment text after ``!``, and
    declaration / PARAMETER / DATA / USE statement lines (those need
    structural handling, not regex). String literals are also masked
    out so an English ``ZERO`` inside an XERBLA-style message stays
    untouched.
    """
    if target_mode.is_kind_based:
        return line
    if not line or is_comment_line(line):
        return line
    if _DECL_KEYWORD_RE.match(line):
        return line

    if renames is None:
        renames = _all_known_constant_renames(target_mode)
    if not renames:
        return line
    renames = {k.upper(): v for k, v in renames.items()}

    # An assignment of the form ``ZERO = <expr>`` typically came from
    # ``convert_parameter_stmts`` keeping a name as a runtime local
    # because it carries complex semantics (e.g. ``COMPLEX(kind=8) ZERO;
    # PARAMETER(ZERO = (0.0D0, 0.0D0))``). Renaming the LHS to
    # ``DD_ZERO`` (a multifloats public constant) would turn the
    # assignment into a write to an imported PARAMETER, which gfortran
    # rejects as "Named constant in variable definition context". Skip
    # the rename for the LHS name on such lines (use sites on the
    # right-hand side are still a problem at the procedure-scope level
    # but are out of scope for this per-line fix).
    assn_m = re.match(r'^\s*([A-Za-z_]\w*)\s*=(?!=)', line)
    if assn_m and assn_m.group(1).upper() in renames:
        renames = {
            k: v for k, v in renames.items()
            if k != assn_m.group(1).upper()
        }
        if not renames:
            return line

    # Split off any inline comment so we don't substitute inside it.
    # We must respect string literals when locating the '!'.
    code_end = len(line)
    in_str, qch = False, ''
    i = 0
    while i < len(line):
        ch = line[i]
        if in_str:
            if ch == qch:
                if i + 1 < len(line) and line[i + 1] == qch:
                    i += 2  # doubled-quote escape stays inside the string
                    continue
                in_str = False
        elif ch in ("'", '"'):
            in_str, qch = True, ch
        elif ch == '!':
            code_end = i
            break
        i += 1
    code, tail = line[:code_end], line[code_end:]

    # Mask out string literal interiors in code segment.
    parts = re.split(r"('(?:[^']|'')*'|\"(?:[^\"]|\"\")*\")", code)
    pattern = _known_constants_pattern(frozenset(renames.keys()))

    def _sub(m):
        return renames[m.group(1).upper()]

    for idx in range(0, len(parts), 2):
        parts[idx] = pattern.sub(_sub, parts[idx])
    return ''.join(parts) + tail


_XERBLA_STR_RE = re.compile(r"'([A-Za-z][A-Za-z0-9_]*)( ?)'")


def replace_xerbla_strings(line: str, rename_map: dict[str, str]) -> str:
    """Replace routine names inside XERBLA string arguments.

    XERBLA's first argument is a Fortran string literal naming the
    failing routine — e.g. ``CALL XERBLA('DGEMM ', 1)``. A routine
    rename must rewrite that literal too. Previously this function
    looped over the whole rename_map and called ``str.replace`` twice
    per entry per line, which is O(N * lines) — catastrophic once the
    map hits MUMPS scale (6k+ entries). The regex version below is
    O(lines) with a dict lookup per quoted identifier found.

    Wiring note: :func:`replace_routine_names` does not mask string
    literals, so by the time this pass runs (always after it) the
    quoted name has normally been renamed already and the lookup here
    misses. What remains is uppercasing the literal when the mapped
    name still resolves (identity-mapped entries). Only the fixed-form
    drivers call this; the free-form drivers deliberately rely on
    ``replace_routine_names`` alone — verified empirically that quoted
    names are renamed there too.
    """
    def sub(m):
        name = m.group(1)
        upper = name.upper()
        new = rename_map.get(upper)
        if new is None:
            return m.group(0)
        return f"'{new.upper()}{m.group(2)}'"
    return _XERBLA_STR_RE.sub(sub, line)
