"""Regression test for the inline-comment / string-literal masker in
``replace_known_constants``.

The masker has to find the inline ``!`` comment marker while staying
inside string literals — including Fortran's doubled-quote escape
(``''`` for a literal apostrophe). A previous version of the loop
used ``for i, ch in enumerate(line): ... continue`` to skip the
doubled-quote escape, which left the inner index pointing at the
*second* quote and prematurely closed the string. A subsequent ``!``
that was actually inside the string was then mis-identified as the
comment marker, truncating the line.
"""

from migrator.fortran_migrator import replace_known_constants
from migrator.target_mode import load_target


def _mf():
    return load_target('multifloats')


def test_doubled_quote_inside_string_keeps_inline_bang_protected():
    target = _mf()
    # The string literal contains a doubled apostrophe ('') followed
    # by a '!'. The masker must treat that '!' as still being inside
    # the string, not as the start of a Fortran comment. The known
    # constant ZERO appears outside the string and must still be
    # rewritten.
    line = "      X = 'foo'' !ZERO' // ZERO"
    out = replace_known_constants(line, target)
    # Comment-detection should not have truncated the trailing
    # `// ZERO` (which lives in the code portion, not a comment).
    assert '//' in out
    # The trailing ZERO outside the string is rewritten.
    assert out.rstrip().endswith('DD_ZERO')
    # The string literal interior is preserved verbatim — the
    # in-string ZERO must NOT be rewritten.
    assert "'foo'' !ZERO'" in out


def test_inline_comment_after_clean_string_still_truncates():
    target = _mf()
    # No doubled-quote shenanigans: the masker must still strip the
    # inline comment so substitutions inside the comment don't fire.
    line = "      X = ZERO   ! ZERO is the additive identity"
    out = replace_known_constants(line, target)
    # The bare ZERO before the comment is rewritten.
    assert 'X = DD_ZERO' in out
    # The ZERO inside the inline comment is NOT rewritten.
    assert '! ZERO is the additive identity' in out


def test_double_quoted_doubled_quote_inside_string_keeps_inline_bang_protected():
    """Same masker invariant as the ``''`` test, but using Fortran's
    other string delimiter — ``"`` — and its doubled-quote escape
    ``""``. The fix in commit 43c6b7c6 covered this path but no test
    exercised it."""
    target = _mf()
    line = '      X = "foo"" !ZERO" // ZERO'
    out = replace_known_constants(line, target)
    assert '//' in out
    assert out.rstrip().endswith('DD_ZERO')
    assert '"foo"" !ZERO"' in out


def test_kind_based_target_is_a_no_op():
    """Sanity: the function gates on ``is_kind_based`` and returns the
    line unchanged for kind10/kind16."""
    target = load_target('kind16')
    line = "      X = ZERO"
    assert replace_known_constants(line, target) == line
