"""Coverage for the canonical Fortran string-literal scanner and the
helpers built on top of it.

Before consolidation each scanner had its own char-walking loop and at
least two of them mishandled the doubled-quote escape (``''`` / ``""``).
These tests pin the corrected behavior of the shared primitives.
"""

from migrator.fortran_migrator import (
    _build_split_mask,
    _count_open_parens,
    _find_inline_bang,
    _iter_outside_strings,
    _strip_strings_and_comments,
)
from migrator.pipeline import _strip_inline_bang


def test_iter_outside_strings_skips_doubled_apostrophe():
    text = "a'b''c'd"
    yielded = ''.join(ch for _, ch in _iter_outside_strings(text))
    # Only the chars truly outside the literal survive.
    assert yielded == 'ad'


def test_iter_outside_strings_skips_doubled_double_quote():
    text = 'a"b""c"d'
    yielded = ''.join(ch for _, ch in _iter_outside_strings(text))
    assert yielded == 'ad'


def test_find_inline_bang_ignores_bang_inside_double_quoted_string():
    line = '      X = "a"" !nope" + Y  ! real'
    pos = _find_inline_bang(line)
    assert line[pos:].lstrip().startswith('! real')


def test_count_open_parens_doubled_quote_with_inner_paren():
    # Doubled-quote escape inside a single-quoted string. Without the
    # fix, _count_open_parens flipped state on the first `'`, treated
    # `'(' '` as code, and returned -1 instead of 0.
    line = "      WRITE(*,*) 'O''(unbalanced)'"
    assert _count_open_parens(line) == 0


def test_count_open_parens_doubled_quote_with_inner_paren_double_delim():
    line = '      WRITE(*,*) "a""(b)"'
    assert _count_open_parens(line) == 0


def test_strip_strings_and_comments_handles_doubled_quote_escape():
    # The bang lives inside the string — the cleaner must drop the
    # whole literal and not treat the inner '!' as a comment marker.
    line = "      X = 'foo''!bar'  Y"
    cleaned = _strip_strings_and_comments(line)
    assert 'foo' not in cleaned
    assert 'bar' not in cleaned
    assert 'Y' in cleaned


def test_pipeline_strip_inline_bang_preserves_double_quoted_bang():
    line = 'CALL FOO("a"" !inside") ! tail'
    out = _strip_inline_bang(line)
    # The trailing real comment is gone; the in-string `!` survives.
    assert 'tail' not in out
    assert '!inside' in out


def test_build_split_mask_marks_doubled_double_quote_unsafe():
    body = 'A "x""y" B'
    mask = _build_split_mask(body)
    # Positions of the entire double-quoted run must be unsafe.
    start = body.index('"')
    end = body.rindex('"')
    for i in range(start, end + 1):
        assert mask[i] is False, f'pos {i} should be in-string'
    # Positions outside (the trailing space and 'B') stay safe.
    assert mask[end + 1] is True
    assert mask[end + 2] is True
