"""Fixed-form line reformatting + statement segmentation.

Re-splits migrated fixed-form lines that overflow column 72 and segments
fixed-form statements at their continuation boundaries. Extracted verbatim
from ``fortran_migrator.py``.
"""
from .lex import _build_split_mask, _find_inline_bang, is_comment_line


def reformat_fixed_line(line: str, cont_char: str = '+') -> str:
    # Preprocessor directives (``#if``, ``#include``, ``#define`` ...) are
    # not bound by fixed-form column 72 and must not be split into
    # continuation lines — doing so produces a truncated directive on
    # the first line (e.g. ``#if A || B ||`` dangling) and a second line
    # the preprocessor doesn't understand. Leave them alone regardless
    # of length.
    if line.lstrip().startswith('#'):
        return line
    if len(line) <= 72 or is_comment_line(line) or (len(line) > 6 and line[6:].lstrip().startswith('!')):
        return line
    # If an inline ``!`` comment sits within the first 72 columns, keep
    # the whole line intact — fixed-form Fortran ignores columns past
    # 72, and we must NOT split across the comment (the text after
    # ``!`` would otherwise land on a continuation line as code).
    # Scan for the first ``!`` outside a string literal. If it sits
    # within col 72 we must keep the line intact — splitting it would
    # land the comment text on a continuation chunk.
    bang = _find_inline_bang(line)
    if bang < len(line) and bang < 72:
        return line
    prefix, body = line[:6] if len(line) >= 6 else line.ljust(6), line[6:]
    safe = _build_split_mask(body)
    chunks = []
    while len(body) > 66:
        split_pos = 66
        for i in range(65, max(35, 65 - 30), -1):
            if body[i] in (',', ' ') and safe[i]:
                split_pos = i + 1
                break
        else:
            for i in range(65, 0, -1):
                if safe[i]:
                    split_pos = i
                    break
        chunks.append(body[:split_pos])
        body, safe = body[split_pos:], safe[split_pos:]
    chunks.append(body)
    result_lines = [prefix + chunks[0]]
    for chunk in chunks[1:]: result_lines.append('     ' + cont_char + chunk)
    return '\n'.join(result_lines)


def _strip_inline_comment(text: str) -> str:
    """Strip an inline ``!`` comment from a Fortran code line, respecting
    string literals. The trailing whitespace before the ``!`` is also
    trimmed.

    Used by :func:`_segment_fixed_form_statements` when joining
    continuation lines: an inline ``!`` mid-statement would otherwise
    swallow every continuation that follows once
    :func:`reformat_fixed_line` re-splits the joined string at column
    66 and the ``!`` lands on a chunk that includes content from the
    next continuation line. The comment is irretrievably lost from the
    joined / reformatted output, which is the price of correctness.
    Single-physical-line statements with ``s == joined`` go through the
    no-transform path and keep their comments verbatim from ``lines``.
    """
    bang = _find_inline_bang(text)
    if bang < len(text):
        return text[:bang].rstrip()
    return text


def _segment_fixed_form_statements(
    physical: list[str],
) -> list[tuple[str, list[str], list[str], str]]:
    """Group physical fixed-form lines into logical statements.

    Each entry is ``(kind, lines, terminators, joined)`` where ``kind`` is
    ``'blank' | 'comment' | 'pp' | 'code'``. ``lines`` and ``terminators``
    are aligned slices of ``physical`` (text without newline / the
    original line terminator). For ``'code'`` statements with continuation
    lines, ``joined`` is the head plus each continuation's column-7+
    content concatenated with single spaces — with any inline ``!``
    comments stripped (see :func:`_strip_inline_comment`). This is what
    the per-line transform passes operate on, so paren-walkers can match
    across the physical line break. For other kinds, ``joined`` is the
    head text (with inline comments preserved — single-line statements
    don't risk the swallow-the-next-continuation failure mode).
    """
    out: list[tuple[str, list[str], list[str], str]] = []
    i = 0
    while i < len(physical):
        raw = physical[i]
        if raw.endswith('\r\n'):
            text, term = raw[:-2], '\r\n'
        elif raw.endswith('\n') or raw.endswith('\r'):
            text, term = raw[:-1], raw[-1]
        else:
            text, term = raw, ''
        if not text.strip():
            out.append(('blank', [text], [term], text))
            i += 1
            continue
        if text[0] in 'Cc*!':
            out.append(('comment', [text], [term], text))
            i += 1
            continue
        if text.lstrip().startswith('#'):
            out.append(('pp', [text], [term], text))
            i += 1
            continue
        # Code head — absorb any immediately-following continuation lines.
        lines = [text]
        terms = [term]
        joined = text
        j = i + 1
        while j < len(physical):
            nxt = physical[j]
            if nxt.endswith('\r\n'):
                ntext, nterm = nxt[:-2], '\r\n'
            elif nxt.endswith('\n') or nxt.endswith('\r'):
                ntext, nterm = nxt[:-1], nxt[-1]
            else:
                ntext, nterm = nxt, ''
            if (len(ntext) > 5 and ntext[:1] != '\t' and ntext[:5] == '     '
                    and ntext[5:6] not in (' ', '0', '\t')):
                lines.append(ntext)
                terms.append(nterm)
                # Strip the previous segment's inline ``!`` comment
                # before appending the next continuation, otherwise the
                # comment swallows every following chunk once
                # ``reformat_fixed_line`` re-splits the joined string.
                # On the first continuation, this also strips any inline
                # comment from the head line.
                joined = _strip_inline_comment(joined) + ' ' + ntext[6:]
                j += 1
            else:
                break
        # Strip any inline comment from the trailing continuation too
        # (so the rebuilt single-line / reformatted output isn't truncated
        # at the comment when reformat_fixed_line walks it).
        if len(lines) > 1:
            joined = _strip_inline_comment(joined)
        out.append(('code', lines, terms, joined))
        i = j
    return out
