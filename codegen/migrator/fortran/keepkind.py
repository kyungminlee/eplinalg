"""Keep-kind sentinel machinery + LWORK round-up stripping (Cluster K).

Protects source spans whose KIND must survive migration by swapping them for
private sentinels before the type rewrite and restoring them after. Self-
contained; extracted verbatim from ``fortran_migrator.py``.
"""
import re

from .lex import is_continuation_line


_KK_SENTINEL = '__KEEPKIND_DP__'


_KK_DBLE_SENTINEL = '__KEEPKIND_DBLE__'


_KK_DCMPLX_SENTINEL = '__KEEPKIND_DCMPLX__'


KEEPKIND_MARKERS = (_KK_SENTINEL, _KK_DBLE_SENTINEL, _KK_DCMPLX_SENTINEL)


def has_keepkind_marker(text: str) -> bool:
    """True when ``text`` carries any keep-kind sentinel (i.e. it is
    mid-pipeline between apply and restore and must not be rewritten)."""
    return any(marker in text for marker in KEEPKIND_MARKERS)


_KK_DP_RE = re.compile(r'DOUBLE\s+PRECISION', re.IGNORECASE)
_KK_DBLE_RE = re.compile(r'\bdble\s*\(', re.IGNORECASE)
_KK_DCMPLX_RE = re.compile(r'\bdcmplx\s*\(', re.IGNORECASE)


def _apply_keep_kind_sentinel(source: str, keep_kind_lines: frozenset[int]) -> str:
    """Replace DP-defining tokens with non-matching sentinels on each
    1-based line number in ``keep_kind_lines``. Protects the line from
    every migrator regex that rewrites those tokens; restored by
    :func:`_restore_keep_kind_sentinel` on the migrated output.

    Protected tokens: ``DOUBLE PRECISION`` (type declaration),
    ``dble(`` (convert-to-DP intrinsic), ``dcmplx(`` (convert-to-DC
    intrinsic). The last two are protected on call-site lines so that
    callers of verbatim (copy_files) DP-stable routines keep passing
    DP values instead of being rewritten to ``REAL(x, KIND=16)``.
    """
    lines = source.splitlines(keepends=True)
    for ln in keep_kind_lines:
        if 1 <= ln <= len(lines):
            t = _KK_DP_RE.sub(_KK_SENTINEL, lines[ln - 1])
            t = _KK_DBLE_RE.sub(_KK_DBLE_SENTINEL + '(', t)
            t = _KK_DCMPLX_RE.sub(_KK_DCMPLX_SENTINEL + '(', t)
            lines[ln - 1] = t
    return ''.join(lines)


def _restore_keep_kind_sentinel(source: str) -> str:
    source = source.replace(_KK_SENTINEL, 'DOUBLE PRECISION')
    source = source.replace(_KK_DBLE_SENTINEL, 'dble')
    source = source.replace(_KK_DCMPLX_SENTINEL, 'dcmplx')
    return source


_ROUNDUP_NAME_RE = re.compile(r'\b([A-Z])ROUNDUP_LWORK\b')

# Orphan-declaration matchers for pass 2 of _strip_roundup_lwork.
# REAL(KIND=N) <list> — strip ROUNDUP token from list.
_ROUNDUP_REAL_DECL_RE = re.compile(
    r'^(\s*REAL(?:\s*\(\s*KIND\s*=\s*\d+\s*\))?\s+)(.*)$', re.IGNORECASE)
_ROUNDUP_EXTERNAL_DECL_RE = re.compile(r'^(\s*EXTERNAL\s+)(.*)$', re.IGNORECASE)
# Modern F90 attribute-list form: ``REAL, EXTERNAL :: NAME``
# and ``DOUBLE PRECISION, EXTERNAL :: NAME``. Migrator-level
# strip needs to handle both to drop the orphan ROUNDUP_LWORK
# declaration after the call sites are stripped.
_ROUNDUP_ATTR_DECL_RE = re.compile(
    r'^(\s*(?:REAL|DOUBLE\s+PRECISION)'
    r'(?:\s*\(\s*KIND\s*=\s*\d+\s*\))?'
    r'\s*,\s*EXTERNAL\s*::\s*)(.*)$',
    re.IGNORECASE)


def _is_cont_free(prev: str) -> bool:
    # Free-form: previous line ended (before any trailing comment) in '&'.
    s = prev.rstrip('\n').rstrip()
    # strip a trailing '! ...' comment for the heuristic
    if '!' in s:
        s = s.split('!', 1)[0].rstrip()
    return s.endswith('&')


def _strip_roundup_lwork(source: str, target_mode) -> str:
    """At kind ≥ 8 the SROUNDUP_LWORK / DROUNDUP_LWORK wrapper is dead
    code: it guards a float→int round-trip bug that only bites when
    LWORK exceeds the float mantissa (2**24 at single precision). Any
    target with mantissa ≥ 32 bits (kind8 and above, plus multifloats)
    represents 32-bit INTEGER exactly, so the IF branch is unreachable
    and the function is the identity.

    Strip the call wrapper and its declarations so:
      - the S-derived re-migration matches the D-derived canonical
        (which already uses bare ``WORK(1)=LWKOPT`` upstream);
      - the dead ``qroundup_lwork.f`` symbol drops out of the link.

    Skipped at kind4 (single precision is exactly where the wrapper
    matters). Operates on the post-rename text, matching any single-
    letter prefix so it works regardless of target prefix (Q at
    kind16, M at multifloats, etc.).
    """
    if target_mode.is_kind_based and (target_mode.kind_suffix or 0) < 8:
        return source
    if not _ROUNDUP_NAME_RE.search(source):
        return source
    source = _strip_roundup_calls(source)
    return _strip_roundup_decls(source)


def _strip_roundup_calls(source: str) -> str:
    """Pass 1: strip call wrappers,
    ``<PREFIX>ROUNDUP_LWORK( <expr> )`` -> ``<expr>``."""
    out = []
    i = 0
    n = len(source)
    while i < n:
        m = _ROUNDUP_NAME_RE.search(source, i)
        if not m:
            out.append(source[i:])
            break
        out.append(source[i:m.start()])
        # Skip the function definition itself (FUNCTION QROUNDUP_LWORK(...))
        prefix = source[max(0, m.start() - 16):m.start()].upper()
        if re.search(r'\bFUNCTION\s*$', prefix):
            out.append(m.group(0))
            i = m.end()
            continue
        # Find following '(' (possibly with whitespace)
        j = m.end()
        while j < n and source[j] in ' \t':
            j += 1
        if j >= n or source[j] != '(':
            # Not a call (e.g. EXTERNAL listing) — keep token, handled below.
            out.append(m.group(0))
            i = m.end()
            continue
        # Walk balanced parens
        depth = 0
        k = j
        while k < n:
            c = source[k]
            if c == '(':
                depth += 1
            elif c == ')':
                depth -= 1
                if depth == 0:
                    break
            k += 1
        if k >= n:
            # Unbalanced — bail, keep original
            out.append(source[m.start():])
            break
        out.append(source[j + 1:k])  # the inner expression
        i = k + 1
    return ''.join(out)


def _strip_roundup_decls(source: str) -> str:
    """Pass 2: drop the now-orphan declarations. Logical-line-aware:
    joins fixed-form ($/+/&) and free-form (& at EOL) continuations,
    strips the dead token from any REAL/EXTERNAL list, then re-emits —
    also dropping lines whose entire list became empty."""
    physical = source.splitlines(keepends=True)

    def is_cont_fixed(line: str) -> bool:
        # Fixed-form: cols 1-5 blank and column 6 has a non-space,
        # non-zero, non-tab character ($/+/& conventionally).
        return is_continuation_line(line, tab_marker=False)

    out_lines: list[str] = []
    i = 0
    while i < len(physical):
        line = physical[i]
        # Comment lines pass through unchanged.
        if line[:1] in ('C', 'c', '*', '!'):
            out_lines.append(line); i += 1; continue
        # Build logical-line group [i, j)
        j = i + 1
        while j < len(physical):
            nxt = physical[j]
            if nxt[:1] in ('C', 'c', '*', '!'):
                # In fixed form, comment lines may interleave continuations;
                # but if comment line, stop the group here.
                break
            if is_cont_fixed(nxt) or _is_cont_free(physical[j - 1]):
                j += 1
            else:
                break
        group = physical[i:j]
        # Reconstruct logical content (strip continuation glue for matching).
        logical = group[0]
        for k in range(1, len(group)):
            ln = group[k]
            # strip continuation marker (col 6 in fixed form) — we just
            # blank out positions 0..5 of continuation lines.
            if is_cont_fixed(ln):
                logical = logical.rstrip('\n').rstrip() + ' ' + ln[6:].lstrip()
            else:
                # free-form: drop trailing & on prior, drop leading & on this
                logical = logical.rstrip('\n').rstrip()
                if logical.endswith('&'):
                    logical = logical[:-1].rstrip()
                cur = ln.lstrip()
                if cur.startswith('&'):
                    cur = cur[1:].lstrip()
                logical = logical + ' ' + cur
        # Match on the joined logical line.
        joined = logical.rstrip('\n').rstrip()

        m_real = _ROUNDUP_REAL_DECL_RE.match(joined)
        m_ext = _ROUNDUP_EXTERNAL_DECL_RE.match(joined)
        m_attr = _ROUNDUP_ATTR_DECL_RE.match(joined)
        is_attr = bool(m_attr)
        is_real = m_real and not m_ext and not is_attr
        is_ext  = bool(m_ext) and not is_attr
        if is_ext or is_real or is_attr:
            m = m_attr if is_attr else (m_ext if is_ext else m_real)
            head, names = m.group(1), m.group(2)
            tokens = [t.strip() for t in names.split(',')]
            keep = [t for t in tokens
                    if t and not re.fullmatch(
                        r'[A-Z]ROUNDUP_LWORK', t, re.IGNORECASE)]
            if len(keep) == len(tokens):
                # Token not present — emit group unchanged.
                out_lines.extend(group)
            elif not keep:
                # Drop entire logical-line group.
                pass
            else:
                # Re-emit as a single line (give up on preserving the
                # original continuation layout — the result is still
                # legal Fortran in either fixed or free form, since
                # EXTERNAL/REAL declaration lines are not column-7
                # constrained).
                # Keep leading whitespace from original line for cosmetics.
                lead = re.match(r'^(\s*)', group[0]).group(1)
                new_logical = head + ', '.join(keep) + '\n'
                # Preserve original leading indent from first phys line.
                # head may not have leading whitespace; ensure it does:
                if not new_logical.startswith(lead):
                    new_logical = lead + new_logical.lstrip()
                out_lines.append(new_logical)
        else:
            out_lines.extend(group)
        i = j
    return ''.join(out_lines)
