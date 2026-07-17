"""Free-form logical-line segmentation + re-wrapping.

Free-form counterpart to ``fixedform.py``'s
``_segment_fixed_form_statements`` / ``reformat_fixed_line``. Lets the
per-line transform passes operate on a *joined* logical line so a pattern
that a source split across ``&`` continuations is still matched as one
string, then re-wraps the migrated result under a column limit.

The build compiles with ``-ffree-line-length-none`` so wrapping is for
readability / diff-parity, not correctness — a joined line that never
gets re-split still compiles. We wrap anyway to keep migrated free-form
sources readable and comparable to their siblings.

Segment tuple shape mirrors the fixed-form one:
``(kind, lines, terms, joined)`` with ``kind`` in
``'blank' | 'comment' | 'pp' | 'code'``.
"""
from .lex import _ends_in_string, _find_inline_bang, _iter_outside_strings

__all__ = [
    'segment_free_form', 'reformat_free_line',
    '_ends_in_string', '_code_and_cont', '_split_term', '_strip_leading_amp',
]


def _split_term(raw: str) -> tuple[str, str]:
    """Split a physical line into (text, line_terminator)."""
    if raw.endswith('\r\n'):
        return raw[:-2], '\r\n'
    if raw.endswith('\n') or raw.endswith('\r'):
        return raw[:-1], raw[-1]
    return raw, ''


def _code_and_cont(text: str) -> tuple[str, bool, bool]:
    """Return ``(code, continues, in_string_at_amp)``.

    ``code`` is ``text`` with any inline ``!`` comment removed and, if the
    statement continues, the trailing ``&`` stripped (right-trimmed).
    ``continues`` is True when the last non-blank code character is ``&``.
    ``in_string_at_amp`` is True when that ``&`` sits inside an open
    character literal (a string split across the continuation).
    """
    bang = _find_inline_bang(text)
    code = text[:bang]
    stripped = code.rstrip()
    if stripped.endswith('&'):
        before = stripped[:-1]
        in_str = _ends_in_string(before)
        # Trailing whitespace before the ``&`` is significant only inside a
        # string literal; drop it otherwise so token-boundary joins are clean.
        if not in_str:
            before = before.rstrip()
        return before, True, in_str
    # No continuation: keep the code (comment already dropped only for the
    # continues case; here we return the original text so single-line
    # statements retain any inline comment verbatim).
    return text, False, False


def _strip_leading_amp(text: str) -> tuple[str, bool]:
    """If the first non-blank char is ``&``, return (content-after-&, True);
    otherwise (text, False)."""
    ls = text.lstrip()
    if ls.startswith('&'):
        return ls[1:], True
    return text, False


def segment_free_form(
    physical: list[str],
) -> list[tuple[str, list[str], list[str], str]]:
    """Group physical free-form lines into logical statements.

    Each entry is ``(kind, lines, terms, joined)``. For ``'code'``
    statements spanning ``&`` continuations, ``joined`` fuses the pieces:
    a continuation whose next line begins with ``&`` is joined tightly (a
    token split across the break — no space), otherwise with a single
    space (the break was at a token boundary). String literals split
    across a continuation are joined tightly and verbatim. Full-line ``!``
    comments interleaved between continuation lines are absorbed into the
    statement's ``lines`` (kept verbatim on output) but contribute nothing
    to ``joined``.
    """
    out: list[tuple[str, list[str], list[str], str]] = []
    i = 0
    n = len(physical)
    while i < n:
        raw = physical[i]
        text, term = _split_term(raw)
        st = text.strip()
        if not st:
            out.append(('blank', [text], [term], text))
            i += 1
            continue
        if st.startswith('#'):
            out.append(('pp', [text], [term], text))
            i += 1
            continue
        if st.startswith('!'):
            out.append(('comment', [text], [term], text))
            i += 1
            continue

        # Code head. Absorb continuation lines while the current piece ends
        # with ``&``.
        lines = [text]
        terms = [term]
        code, cont, in_str = _code_and_cont(text)
        joined = code
        j = i + 1
        while cont and j < n:
            nxt = physical[j]
            ntext, nterm = _split_term(nxt)
            nst = ntext.strip()
            if not nst:
                # Blank line: not a legal continuation body. Stop absorbing;
                # leave it for the next iteration.
                break
            if nst.startswith('!'):
                # Comment line between continuations — allowed. Keep it with
                # the statement but do not fold it into ``joined``.
                lines.append(ntext)
                terms.append(nterm)
                j += 1
                continue
            if nst.startswith('#'):
                # A preprocessor directive cannot be folded into a statement;
                # stop here and let it segment on its own.
                break
            content, had_leading = _strip_leading_amp(ntext)
            ncode, ncont, _ = _code_and_cont(content)
            if in_str:
                # String continues across the break — verbatim, tight.
                joined = joined + ncode
            elif had_leading:
                # Explicit resume marker: a token was split — fuse tightly.
                joined = joined.rstrip() + ncode.lstrip()
            else:
                # Break at a token boundary — a space stands in for it.
                joined = joined.rstrip() + ' ' + ncode.strip()
            lines.append(ntext)
            terms.append(nterm)
            cont = ncont
            in_str = _ends_in_string(joined) if ncont else False
            j += 1
        out.append(('code', lines, terms, joined))
        i = j
    return out


def reformat_free_line(line: str, width: int = 132, indent_extra: int = 4) -> str:
    """Re-wrap a long free-form logical line under ``width`` columns using
    trailing-``&`` continuations. Comment / preprocessor lines and lines
    already within ``width`` are returned unchanged. Breaks only at
    ``,``/`` `` outside string literals so no token or string is split;
    each continuation after the first carries a leading ``&`` as well so a
    break inside (say) a long unbroken token region stays legal."""
    if len(line) <= width:
        return line
    ls = line.lstrip()
    if ls.startswith('!') or ls.startswith('#'):
        return line
    # An inline comment within the limit means the significant code already
    # fits; leave the line intact rather than risk splitting across the ``!``.
    bang = _find_inline_bang(line)
    if bang < width:
        return line

    indent = len(line) - len(ls)
    cont_prefix = ' ' * (indent + indent_extra) + '&'
    # Safe-split mask: positions not inside a string literal.
    safe = [False] * len(line)
    for idx, _ in _iter_outside_strings(line):
        safe[idx] = True

    n = len(line)
    chunks: list[tuple[bool, str]] = []  # (is_first, raw content slice)
    pos = 0
    first = True
    while pos < n:
        prefix_len = 0 if first else len(cont_prefix)
        # If the remainder fits on this physical line (no trailing ``&``
        # needed), emit it as the final chunk.
        if n - pos <= width - prefix_len:
            chunks.append((first, line[pos:]))
            break
        # Otherwise reserve two columns for a trailing `` &`` and find the
        # last safe break within budget.
        hi = pos + (width - prefix_len - 2)
        split = -1
        lo = pos + 1
        for k in range(min(hi, n) - 1, lo, -1):
            if line[k] in (',', ' ') and safe[k]:
                split = k + 1
                break
        if split < 0:
            # No safe break in the window — extend outward rather than split
            # a token or string literal.
            for k in range(hi, n):
                if line[k] in (',', ' ') and safe[k]:
                    split = k + 1
                    break
        if split < 0:
            chunks.append((first, line[pos:]))  # unbreakable tail
            break
        chunks.append((first, line[pos:split]))
        pos = split
        first = False
    if len(chunks) == 1:
        return line
    out: list[str] = []
    for idx, (is_first, content) in enumerate(chunks):
        is_last = idx == len(chunks) - 1
        text = content.rstrip() if is_first else cont_prefix + content.strip()
        if not is_last:
            text = text + ' &'
        out.append(text)
    # If an unbreakable segment (a long token or string literal) left any
    # physical line still over the limit, wrapping bought nothing — emit the
    # original single long line instead (legal under -ffree-line-length-none).
    if any(len(p) > width for p in out):
        return line
    return '\n'.join(out)
