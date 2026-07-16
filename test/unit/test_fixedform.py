"""Unit tests for fixed-form logical-line segmentation + re-wrapping."""
from migrator.fortran.fixedform import (
    _segment_fixed_form_statements,
    reformat_fixed_line,
)


def _seg(src: str):
    return _segment_fixed_form_statements(src.splitlines(keepends=True))


def test_single_code_line():
    segs = _seg("      X = A + B\n")
    assert len(segs) == 1
    kind, lines, terms, joined = segs[0]
    assert kind == 'code'
    assert joined == "      X = A + B"


def test_blank_comment_pp_kinds():
    segs = _seg("\nC a comment\n#if FOO\n      Y = 1\n")
    assert [s[0] for s in segs] == ['blank', 'comment', 'pp', 'code']


def test_token_boundary_join():
    # A statement split at a token boundary keeps a blank between the pieces.
    joined = _seg("      CALL FOO(A,\n     &B, C)\n")[0][3]
    assert "FOO(A, B, C)".replace(' ', '') in joined.replace(' ', '')


def test_string_split_across_continuation_not_doubled():
    # A character literal split across a continuation must join *tight*: the
    # continuation's column-7+ text follows the head with no inserted blank,
    # or the literal gains a spurious space. Regression: the segmenter used
    # to unconditionally insert a token-boundary blank, doubling the space in
    # WRITE messages such as "...parallel analysis..." -> "...parallel  analysis...".
    src = (
        "      WRITE(6,*) '(\"too small for the parallel\n"
        "     & analysis. Reverting\")'\n"
    )
    joined = _seg(src)[0][3]
    assert "parallel  analysis" not in joined
    assert "parallel analysis" in joined


def test_string_split_no_leading_space_tight():
    # Break falls mid-word inside the literal (no leading blank on the
    # continuation): the two halves must fuse with no gap.
    src = (
        "      MSG = 'DMUMPS_FAC_ASM\n"
        "     &_NIV2'\n"
    )
    joined = _seg(src)[0][3]
    assert "DMUMPS_FAC_ASM_NIV2" in joined
    assert "DMUMPS_FAC_ASM _NIV2" not in joined


def test_reformat_short_line_unchanged():
    line = "      X = A + B"
    assert reformat_fixed_line(line) == line


def test_reformat_does_not_split_string_literal():
    # A long single string with no safe break point outside it must not be
    # split mid-literal; the reformatter falls back to the whole line.
    body = "      MSG = '" + "x" * 90 + "'"
    out = reformat_fixed_line(body)
    # Every emitted physical line, rejoined at col 7+, must reproduce the
    # original string content unchanged.
    rejoined = _seg(out + '\n')[0][3] if '\n' in out else out
    assert rejoined.replace(' ', '') == body.replace(' ', '')


def test_reformat_leaves_preprocessor():
    long_pp = "#define FOO " + "b" * 90
    assert reformat_fixed_line(long_pp) == long_pp
