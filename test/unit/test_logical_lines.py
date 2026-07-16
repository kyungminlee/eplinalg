"""Unit tests for free-form logical-line join/wrap primitives."""
from migrator.fortran.logical_lines import (
    segment_free_form,
    reformat_free_line,
    _ends_in_string,
    _code_and_cont,
)


def _seg(src: str):
    return segment_free_form(src.splitlines(keepends=True))


def test_single_code_line_passthrough():
    segs = _seg("      X = A + B\n")
    assert len(segs) == 1
    kind, lines, terms, joined = segs[0]
    assert kind == 'code'
    assert joined == "      X = A + B"
    assert lines == ["      X = A + B"]


def test_blank_comment_pp_kinds():
    segs = _seg("\n! a comment\n#if FOO\n   Y = 1\n")
    kinds = [s[0] for s in segs]
    assert kinds == ['blank', 'comment', 'pp', 'code']


def test_token_boundary_join_gets_space():
    src = "      X = A + &\n          B * C\n"
    kind, lines, terms, joined = _seg(src)[0]
    assert kind == 'code'
    assert joined == "      X = A + B * C"
    assert len(lines) == 2


def test_split_token_leading_amp_joins_tight():
    src = "      CALL FOO(MPI_DOUBLE_PRE&\n     &CISION, X)\n"
    joined = _seg(src)[0][3]
    assert "MPI_DOUBLE_PRECISION" in joined
    assert "PRE CISION" not in joined


def test_three_line_continuation():
    src = (
        "      CALL MPI_REDUCE(COL, WRK(1+M), N, &\n"
        "          MPI_DOUBLE_PRECISION, &\n"
        "          MPI_MAX, 0, COMM, IERR)\n"
    )
    joined = _seg(src)[0][3]
    assert joined.count('&') == 0
    assert "MPI_REDUCE(COL, WRK(1+M), N, MPI_DOUBLE_PRECISION, MPI_MAX, 0, COMM, IERR)" in joined


def test_inline_comment_after_amp_dropped_from_join():
    src = "      X = A + &  ! keep going\n          B\n"
    joined = _seg(src)[0][3]
    assert '!' not in joined
    assert joined == "      X = A + B"


def test_comment_line_between_continuations_absorbed_not_joined():
    src = (
        "      X = A + &\n"
        "! interjected comment\n"
        "          B\n"
    )
    kind, lines, terms, joined = _seg(src)[0]
    assert kind == 'code'
    assert len(lines) == 3
    assert joined == "      X = A + B"


def test_string_split_across_continuation_tight():
    src = "      MSG = 'HELLO &\n     &WORLD'\n"
    joined = _seg(src)[0][3]
    assert "HELLO WORLD" in joined


def test_ends_in_string():
    assert _ends_in_string("X = 'abc")
    assert not _ends_in_string("X = 'abc'")
    assert not _ends_in_string("X = A + B")


def test_code_and_cont_detects_amp():
    code, cont, in_str = _code_and_cont("      X = A + &")
    assert cont and not in_str
    assert code == "      X = A +"
    code, cont, in_str = _code_and_cont("      X = A + B")
    assert not cont


def test_reformat_short_line_unchanged():
    line = "      X = A + B"
    assert reformat_free_line(line, width=132) == line


def test_reformat_wraps_long_line():
    body = "      CALL FOO(" + ", ".join(f"ARG{i}" for i in range(40)) + ")"
    assert len(body) > 132
    out = reformat_free_line(body, width=132)
    parts = out.split('\n')
    assert len(parts) > 1
    for p in parts[:-1]:
        assert p.rstrip().endswith('&')
    # No physical line exceeds the width.
    for p in parts:
        assert len(p) <= 132
    # Round-trips back to the same logical content when re-joined.
    rejoined = _seg(out + '\n')[0][3]
    assert rejoined.replace(' ', '') == body.replace(' ', '')


def test_reformat_does_not_split_string_literal():
    body = "      MSG = '" + "x" * 200 + "'"
    out = reformat_free_line(body, width=132)
    # The long string is unbreakable at a safe point → emitted as-is.
    assert out == body


def test_reformat_leaves_comment_and_pp():
    long_comment = "! " + "a" * 200
    assert reformat_free_line(long_comment, width=132) == long_comment
    long_pp = "#define FOO " + "b" * 200
    assert reformat_free_line(long_pp, width=132) == long_pp
