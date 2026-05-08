"""Migration pipeline — orchestrates the full migration for any library.

Usage:
    from migrator.pipeline import run_migration
    run_migration(recipe_path, output_dir, target_mode=16)
"""

import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from .config import RecipeConfig, load_recipe
from .symbol_scanner import scan_symbols
from .prefix_classifier import classify_symbols, build_rename_map
from .fortran_migrator import (
    _find_inline_bang,
    migrate_file,
    migrate_file_to_string,
    target_filename,
)
from .c_migrator import migrate_c_directory, migrate_c_file_to_string, _build_sub_vars, _expand_template
from .target_mode import TargetMode

from tqdm import tqdm


def _apply_extra_renames(rename_map: dict[str, str],
                         config: RecipeConfig,
                         target_mode: TargetMode) -> dict[str, str]:
    """Append recipe-declared ``extra_renames`` to ``rename_map``.

    Targets may use ``{RP}/{CP}/{RPU}/{CPU}`` template substitutions.
    Returns ``rename_map`` mutated in place. See
    :class:`RecipeConfig.extra_renames` for the canonical use case
    (precision-prefixed orphan symbols whose S/C sibling does not
    exist upstream so the classifier cannot pair them).
    """
    if not config.extra_renames:
        return rename_map
    template_vars = _build_sub_vars(target_mode)
    for src_upper, tgt_template in config.extra_renames.items():
        rename_map[src_upper] = _expand_template(tgt_template, template_vars).upper()
    return rename_map


def _strip_fortran_comments(text: str, ext: str) -> str:
    """Strip comments from Fortran source so convergence comparisons
    ignore commentary that differs between S/D or C/Z source halves.

    Fixed-form (.f/.for): a line is a full comment when column 1 is
    C, c, *, or !. Inline ! comments are stripped unless inside a
    string literal. (D/d-lines are conditionally-compiled debug code,
    not comments — both halves contain or omit them symmetrically, so
    leave them in place.)

    Free-form (.f90/.F90/.f95): a line is a full comment when the
    first non-blank char is !. Inline ! comments are stripped the
    same way.

    Trailing whitespace and blank lines that result from stripping are
    removed so normalization is stable.
    """
    is_fixed = ext.lower() in ('.f', '.for')
    out_lines: list[str] = []
    for line in text.split('\n'):
        if is_fixed and line[:1] in ('C', 'c', '*', '!'):
            # Column-1 comment marker in fixed-form.
            continue
        if not is_fixed and line.lstrip().startswith('!'):
            continue
        # Strip inline ! comments (respect simple single/double quotes)
        stripped = _strip_inline_bang(line)
        if stripped.strip():
            out_lines.append(stripped.rstrip())

    # Join fixed-form continuation lines so sibling sources that wrap
    # at different columns still normalize to the same text. In
    # fixed-form, column 6 (non-blank, non-zero) indicates a
    # continuation of the previous line — merge it into its parent
    # after stripping the continuation marker and surrounding
    # whitespace.
    if is_fixed:
        joined: list[str] = []
        for ln in out_lines:
            if len(ln) > 5 and ln[5] not in (' ', '0', '\t') and joined:
                joined[-1] = joined[-1].rstrip() + ' ' + ln[6:].lstrip()
            else:
                joined.append(ln)
        out_lines = joined
    else:
        # Free-form: '&' at end of line continues onto next line.
        joined: list[str] = []
        pending = ''
        for ln in out_lines:
            if pending:
                ln = pending + ' ' + ln.lstrip()
                pending = ''
            stripped = ln.rstrip()
            if stripped.endswith('&'):
                pending = stripped[:-1].rstrip()
            else:
                joined.append(stripped)
        if pending:
            joined.append(pending)
        out_lines = joined
    return '\n'.join(out_lines)


def _strip_inline_bang(line: str) -> str:
    return line[:_find_inline_bang(line)]


def _strip_real_cmplx_casts(text: str) -> str:
    """Strip ``REAL(expr, KIND=N)`` and ``CMPLX(expr, KIND=N)`` wrappers,
    replacing them with just ``expr``. Single-argument ``REAL(expr)``
    (integer/implicit-real promotion) is also stripped, but
    ``CMPLX(expr)`` is NOT — it materially changes the result type
    from real to complex and the two halves genuinely may differ in
    how they construct complex literals."""
    pattern = re.compile(r'\b(?:REAL|CMPLX)\s*\(', re.IGNORECASE)
    out = []
    i = 0
    while i < len(text):
        m = pattern.search(text, i)
        if not m:
            out.append(text[i:])
            break
        out.append(text[i:m.start()])
        # Find matching ')'
        depth = 1
        j = m.end()
        while j < len(text) and depth > 0:
            c = text[j]
            if c == '(':
                depth += 1
            elif c == ')':
                depth -= 1
                if depth == 0:
                    break
            j += 1
        if depth != 0:
            # Unmatched — give up on this match
            out.append(text[m.start():m.end()])
            i = m.end()
            continue
        inner = text[m.end():j]
        # Recurse into the inner expression so nested casts are
        # stripped even when the outer call doesn't match the KIND
        # shape we rewrite.
        inner = _strip_real_cmplx_casts(inner)
        parts = _split_top_level_comma(inner)
        is_real_call = text[m.start():m.end()].upper().lstrip().startswith('REAL')
        # Skip declaration-form ``REAL(KIND=N)`` — that's a type
        # specifier, not a cast.
        is_kind_decl = (len(parts) == 1 and
                        re.fullmatch(r'\s*KIND\s*=\s*\d+\s*', parts[0]))
        # Strip trailing ``KIND=N`` from any multi-arg call. This
        # covers ``CMPLX(re, im, KIND=16)`` as well as
        # ``REAL(x, KIND=16)``.
        is_kinded_multi = (len(parts) >= 2 and
                           re.fullmatch(r'\s*KIND\s*=\s*\d+\s*', parts[-1]))
        if is_kinded_multi:
            # Drop the KIND= tail argument but keep the call.
            call_head = text[m.start():m.end()]
            new_inner = ','.join(p.strip() for p in parts[:-1])
            if len(parts) == 2:
                # 1-data-arg form: reduce to bare identifier/expr.
                inner_expr = parts[0].strip()
                if _has_top_level_operator(inner_expr):
                    out.append('(' + inner_expr + ')')
                else:
                    out.append(inner_expr)
            else:
                # Multi-arg CMPLX(a, b, KIND=N) → CMPLX(a, b)
                out.append(call_head + new_inner + ')')
            i = j + 1
            continue
        if not is_kind_decl and (
                (len(parts) == 1 and is_real_call)):
            # Wrap the stripped expression in parens ONLY when the
            # inner contains a top-level additive operator — the
            # surrounding context's precedence then matches the
            # unwrapped side (``REAL(N-1, KIND=16)*T`` → ``(N-1)*T``
            # rather than ``N-1*T``). Single atoms keep their natural
            # form so ``REAL(X, KIND=16)`` → ``X``.
            #
            # Single-arg REAL(X)/CMPLX(X) is also stripped: it is an
            # integer→real promotion that the Fortran compiler would
            # insert implicitly when X is used in a real context, so
            # ``REAL(MAXITR)*UNFL`` is semantically identical to
            # ``MAXITR*UNFL`` after migration.
            inner_expr = parts[0].strip()
            if _has_top_level_operator(inner_expr):
                out.append('(' + inner_expr + ')')
            else:
                out.append(inner_expr)
        else:
            # Not the shape we want — keep original call, but with
            # inner expression already stripped recursively.
            out.append(text[m.start():m.end()] + inner + ')')
        i = j + 1
    return ''.join(out)


def _has_top_level_operator(s: str) -> bool:
    """True if ``s`` contains a top-level ``+`` or ``-`` outside nested
    parens/brackets (ignoring unary signs and exponent markers).

    We only consider additive operators here, not ``*`` / ``/``:
    when the stripped expression is used as a factor in an outer
    multiplication (the usual case for REAL(int_expr, KIND=16)*x),
    a pure product like ``MAXITR*Q*Q`` is associative with the
    outer ``*`` so parens are unnecessary. Only ``a-b`` or ``a+b``
    would change precedence under multiplication.
    """
    depth = 0
    for i, c in enumerate(s):
        if c in '([':
            depth += 1
        elif c in ')]':
            depth -= 1
        elif depth == 0 and c in '+-' and i > 0:
            # Skip if preceded by 'E' / 'D' (exponent sign) or an
            # operator (meaning it's a unary sign).
            prev = s[i - 1]
            if prev in 'eEdD' and i >= 2 and s[i - 2].isdigit():
                continue
            if prev in '+-*/=(,':
                continue
            return True
    return False


def _split_top_level_comma(s: str) -> list[str]:
    """Split on commas that are not inside nested parens or brackets."""
    out: list[str] = []
    depth = 0
    start = 0
    for i, c in enumerate(s):
        if c in '([':
            depth += 1
        elif c in ')]':
            depth -= 1
        elif c == ',' and depth == 0:
            out.append(s[start:i])
            start = i + 1
    out.append(s[start:])
    return out


def _light_normalize(text: str) -> str:
    """Minimal normalization for post-migration convergence checking.

    Applies ONLY cosmetic rules that are never the migrator's job:

    * uppercase everything (Fortran is case-insensitive)
    * collapse runs of horizontal whitespace to a single space
    * strip whitespace adjacent to punctuation (commas, parens,
      arithmetic/relational operators)
    * merge ``END IF``/``END DO``/... → ``ENDIF``/``ENDDO``/...
    * sort ``INTRINSIC``/``EXTERNAL`` argument lists (cross-precision
      co-family halves routinely list the same names in different
      orders — purely cosmetic upstream drift)
    * drop blank lines and trailing whitespace

    No prefix collapse, no literal-form rewrites, no type-cast
    stripping — convergence here means the migrator got it right and
    the co-family halves already agree on identifiers, literal kinds,
    and casts.
    """
    text = text.upper()
    text = re.sub(
        r'\bEND\s+('
        r'IF|DO|SELECT|WHERE|FORALL|SUBROUTINE|FUNCTION|'
        r'MODULE|INTERFACE|PROGRAM|TYPE|BLOCK|ASSOCIATE'
        r')\b',
        r'END\1', text,
    )
    text = re.sub(r'[ \t]+', ' ', text)
    # Horizontal whitespace only ([ \t], not \s) so that a Fortran
    # fixed-form DATA statement ``DATA Z,T/0.E0_16,2.E0_16/`` never
    # gets glued onto the next line via its trailing ``/``.
    text = re.sub(r'[ \t]*([*+\-/=])[ \t]*', r'\1', text)
    text = re.sub(r'[ \t]*\([ \t]*', '(', text)
    text = re.sub(r'[ \t]*\)', ')', text)
    text = re.sub(r'[ \t]*,[ \t]*', ',', text)

    # Fold ``REAL(<expr>,KIND=N)`` and ``CMPLX(<expr>,KIND=N)`` casts to
    # bare ``<expr>``. Co-family halves disagree cosmetically when one
    # writes the explicit kind cast and the other relies on implicit
    # promotion (e.g. D writes ``REAL(NZ,KIND=16)*SAFMIN`` while S
    # writes ``NZ*SAFMIN``). Both expressions evaluate identically when
    # the surrounding type is REAL(KIND=16) — the cast is a no-op. The
    # fold runs on both halves before the diff so the cosmetic drift
    # disappears. Only top-level ``,KIND=N`` matches; parameter-form
    # ``REAL(KIND=N)`` declarations are left alone (no leading expr).
    def _strip_kind_casts(s: str) -> str:
        out: list[str] = []
        pat = re.compile(r'\b(REAL|CMPLX)\(')
        i = 0
        while i < len(s):
            m = pat.search(s, i)
            if not m:
                out.append(s[i:]); break
            out.append(s[i:m.start()])
            depth = 1
            comma_pos = -1
            j = m.end()
            while j < len(s) and depth > 0:
                c = s[j]
                if c == '(':
                    depth += 1
                elif c == ')':
                    depth -= 1
                    if depth == 0:
                        break
                elif c == ',' and depth == 1 and comma_pos == -1:
                    comma_pos = j
                j += 1
            if depth != 0:
                # Unbalanced — keep original prefix and stop folding.
                out.append(s[m.start():])
                return ''.join(out)
            if comma_pos != -1 and re.match(
                    r'KIND=\d+$', s[comma_pos + 1:j]):
                out.append(s[m.end():comma_pos])  # the inner expression
            else:
                out.append(s[m.start():j + 1])
            i = j + 1
        return ''.join(out)
    text = _strip_kind_casts(text)

    # Strip per-line whitespace *before* sorting declarations so the
    # ``^`` anchor sees the statement keyword at the line start.
    text = '\n'.join(ln.strip() for ln in text.split('\n'))

    # Sort INTRINSIC / EXTERNAL argument lists. The ordering of these
    # declarations varies between co-family halves (e.g. C sources
    # list ``INTRINSIC CONJG,MAX,MIN,REAL`` while Z sources write
    # ``INTRINSIC REAL,CONJG,MAX,MIN``) but is semantically irrelevant.
    def _sort_decl(m):
        return f'{m.group(1)} ' + ','.join(sorted(m.group(2).split(',')))

    text = re.sub(
        r'\b(INTRINSIC|EXTERNAL)\s+([A-Z_][A-Z0-9_]*(?:,[A-Z_][A-Z0-9_]*)*)',
        _sort_decl, text,
    )

    # Sort bare-identifier lists after a simple type-spec. Co-family
    # halves reorder symbols between their ``INTEGER``/``REAL(KIND=N)``
    # declarations of named external functions (``REAL(KIND=16)
    # XLANGB,XLANTB,QLAMCH`` vs ``QLAMCH,XLANGB,XLANTB``). Match only
    # pure name lists — lines with ``=``, ``(``, ``*`` (e.g. array
    # specs ``A(LDA,*)``) are left alone.
    def _prec_neutral_key(name: str) -> str:
        """Sort key that maps S↔D and C↔Z to the same character so
        precision-equivalent locals (CI vs ZI) sort to the same position."""
        return name.replace('S', 'D').replace('C', 'Z')

    def _sort_typed_decl(m):
        type_spec, body = m.group(1), m.group(2)
        names = body.split(',')
        return f'{type_spec} ' + ','.join(
            sorted(names, key=_prec_neutral_key))

    text = re.sub(
        r'^(INTEGER|LOGICAL|REAL\(KIND=\d+\)|COMPLEX\(KIND=\d+\)|TYPE\(\w+\))\s+'
        r'([A-Z_][A-Z0-9_]*(?:,[A-Z_][A-Z0-9_]*)+)$',
        _sort_typed_decl, text, flags=re.MULTILINE,
    )

    # Sort type-declaration statement lines to eliminate ordering drift
    # between co-family halves. Only lines that start with a type
    # keyword are reordered; everything else stays in place.
    lines = text.split('\n')
    _DECL_RE = re.compile(
        r'^\s*(?:INTEGER|LOGICAL|REAL\(|COMPLEX\(|CHARACTER|'
        r'DOUBLE\s+PRECISION|TYPE\s*\()\b'
    )
    decl_indices = [k for k, ln in enumerate(lines) if _DECL_RE.match(ln)]
    if decl_indices:
        decl_contents = sorted(
            (lines[k] for k in decl_indices),
            key=_prec_neutral_key,
        )
        for idx, k in enumerate(decl_indices):
            lines[k] = decl_contents[idx]

    lines = [ln.strip() for ln in lines]
    return '\n'.join(ln for ln in lines if ln)


def _apply_local_renames(text: str, renames: dict[str, str]) -> str:
    """Apply recipe-declared local-variable equivalence folds.

    ``renames`` declares equivalence classes: each key is a local
    identifier that should be treated as identical to its value for
    convergence-report purposes only. Applied **to both halves** of
    every pair (D/Z canonical read from disk and S/C other migrated
    in-memory) before light-normalized comparison, so ScaLAPACK's
    S-half ``CR`` / ``CI`` folds onto Z-half ``ZR`` / ``ZI`` in
    pzlattrs AND a D-half file that already uses ``CR`` still
    converges with its S-half sibling. Neither the on-disk canonical
    nor the recipe source is modified; this is comparison-only.
    Matching is case-insensitive; replacement casing tracks the
    source match.
    """
    if not renames:
        return text
    keys = sorted(renames.keys(), key=len, reverse=True)
    pattern = re.compile(
        r'\b(' + '|'.join(re.escape(k) for k in keys) + r')\b',
        re.IGNORECASE,
    )
    upper_map = {k.upper(): v.upper() for k, v in renames.items()}

    def _sub(m: re.Match) -> str:
        src = m.group(0)
        new = upper_map[src.upper()]
        if src.isupper():
            return new
        if src.islower():
            return new.lower()
        return new

    return pattern.sub(_sub, text)


_PRECISION_PAIRS = frozenset({
    frozenset(('S', 'D')),
    frozenset(('C', 'Z')),
})


def _is_precision_equivalent_token(a: str, b: str) -> bool:
    """True if *a* and *b* differ only at S↔D or C↔Z positions.

    Both tokens must have the same length, and every position where
    they disagree must be an (S,D) or (C,Z) pair (case-insensitive).
    Returns False for identical tokens (no drift to fold).
    """
    if len(a) != len(b) or a == b:
        return False
    for ca, cb in zip(a, b):
        if ca == cb:
            continue
        if frozenset((ca.upper(), cb.upper())) not in _PRECISION_PAIRS:
            return False
    return True


def _is_precision_equivalent_line(line_a: str, line_b: str) -> bool:
    """True if two normalized source lines differ only at precision-prefix
    positions within identifier tokens.

    Tokenizes both lines into words (``\\w+``) and punctuation (``\\S``),
    then checks each token pair.  If any non-precision token differs, or
    if the token counts don't match, returns False.
    """
    toks_a = re.findall(r'\w+|\S', line_a)
    toks_b = re.findall(r'\w+|\S', line_b)
    if len(toks_a) != len(toks_b):
        return False
    has_drift = False
    for ta, tb in zip(toks_a, toks_b):
        if ta == tb:
            continue
        if not _is_precision_equivalent_token(ta, tb):
            return False
        has_drift = True
    return has_drift


def _filter_precision_drift(lines_a: list[str],
                            lines_b: list[str]) -> list[str]:
    """Compute diff between *lines_a* and *lines_b*, suppressing hunks
    that are purely precision-prefix local-variable drift (S↔D, C↔Z).

    Uses :func:`difflib.SequenceMatcher` to align the two sides, then
    for each ``'replace'`` block where every replaced line pair is
    precision-equivalent (per :func:`_is_precision_equivalent_line`),
    the block is silently dropped.  Insert, delete, and mixed blocks
    are kept as real divergences.

    Returns a list of diff lines (prefixed with ``-`` / ``+``) for
    genuine divergences only.
    """
    import difflib
    sm = difflib.SequenceMatcher(None, lines_a, lines_b)
    diff_lines: list[str] = []
    for op, i1, i2, j1, j2 in sm.get_opcodes():
        if op == 'equal':
            continue
        if op == 'replace' and (i2 - i1) == (j2 - j1):
            if all(
                lines_a[i1 + k] == lines_b[j1 + k]  # identical
                or _is_precision_equivalent_line(lines_a[i1 + k],
                                                 lines_b[j1 + k])
                for k in range(i2 - i1)
            ):
                continue  # pure precision drift — suppress
        for line in lines_a[i1:i2]:
            diff_lines.append(f'-{line}')
        for line in lines_b[j1:j2]:
            diff_lines.append(f'+{line}')
    return diff_lines


def _light_normalize_c(text: str) -> str:
    """Minimal C normalization for post-migration convergence checking.

    Mirrors :func:`_light_normalize` for Fortran. Assumes the migrator
    has already rewritten identifiers, types, and MPI constants — the
    only drift tolerated here is whitespace around tokens and blank
    lines. Preserves case (C is case-sensitive) but collapses runs of
    horizontal whitespace and strips whitespace adjacent to common
    punctuation, so column-aligned declarations like ``float
    *A,`` vs ``QREAL          *A,`` converge.
    """
    # Tokenize each line into words and single non-whitespace chars,
    # then rejoin with a single space. This collapses column padding
    # and incidental space drift around punctuation uniformly.
    lines = [' '.join(re.findall(r'\w+|\S', ln)) for ln in text.split('\n')]
    return '\n'.join(ln for ln in lines if ln)


def _canonicalize_for_compare(text: str) -> str:
    """Normalize text to ignore cosmetic precision-specific differences
    that remain after migration.
    """
    # 0. Strip module-based type constructors: e.g. real64x2('1.0D0') -> '1.0D0'
    # Matches identifiers containing a digit (type names like real64x2,
    # cmplx64x2, real256, etc.) followed by a single-argument call.
    # This avoids matching Fortran keywords or subroutine calls.
    text = re.sub(r"\b(\w*\d\w*)\s*\(\s*'([^']+)'\s*\)", r"\2", text, flags=re.IGNORECASE)
    text = re.sub(r"\b(\w*\d\w*)\s*\(([^)]+)\)", r"\2", text, flags=re.IGNORECASE)

    # 1. d/D exponent → E in numeric literals
    text = re.sub(r'(\d)[dD]([+-]?\d)', r'\1E\2', text)
    # 2. Strip _N kind suffixes on literals (suffix after a digit)
    text = re.sub(r'(\d)_\d+\b', r'\1', text)
    # 3. Drop insignificant exponents (E0/E+0/e-0/E00 ...)
    text = re.sub(r'([0-9.])[eE][+-]?0+\b', r'\1', text)
    # 4. Canonicalize mantissa exponent marker to uppercase E
    text = re.sub(r'(\d)e([+-]?\d)', r'\1E\2', text)
    # 4b. Strip REAL(expr, KIND=N) / CMPLX(expr, KIND=N) wrappers.
    text = _strip_real_cmplx_casts(text)
    # 4c. Canonicalize numeric literal forms.
    text = re.sub(r'(\d)\.0*(?![\deEdD_])', r'\1', text)
    # 4d. Strip ``(KIND=N)`` suffix or TYPE(...) specifier.
    text = re.sub(
        r'\bREAL\s*\(\s*KIND\s*=\s*\d+\s*\)',
        r'REAL', text, flags=re.IGNORECASE)
    text = re.sub(
        r'\bCOMPLEX\s*\(\s*KIND\s*=\s*\d+\s*\)',
        r'COMPLEX', text, flags=re.IGNORECASE)
    # Canonicalize TYPE(...) declarations to REAL/COMPLEX.
    # Any TYPE(...) that remains after KIND canonicalization is a module-
    # based type (e.g. TYPE(real64x2), TYPE(real256)) — normalize to
    # base precision for comparison purposes.
    text = re.sub(
        r'\bTYPE\s*\(\s*\w+\s*\)',
        r'REAL', text, flags=re.IGNORECASE)
    
    # 5. Canonicalize prefix-dependent identifiers
    text = re.sub(
        r'(?<![A-Za-z0-9])[sdczSDCZ]+(?=[A-Za-z])', '@', text,
    )
    # 5a. iso_fortran_env kind imports: ``only: real32`` vs
    #     ``only: real64`` survive as an unused import (after the
    #     migrator rewrites ``WP = real32`` → ``WP = 16``). Normalize
    #     the kind-suffix digits so both halves compare equal.
    text = re.sub(r'\breal(32|64|128)\b', 'realK', text, flags=re.IGNORECASE)
    # 5b. Uppercase everything. Fortran is case-insensitive, and S vs
    #     D sources are not consistent about casing keywords
    #     (``then`` vs ``THEN``, ``IF`` vs ``if``). Uppercasing both
    #     halves makes them compare equal.
    text = text.upper()
    # 5c. Merge ``END IF`` → ``ENDIF`` etc. LAPACK writes these block
    #     terminators with or without a space arbitrarily between
    #     precision halves.
    text = re.sub(r'\bEND\s+(IF|DO|SELECT|WHERE|FORALL|SUBROUTINE|FUNCTION|MODULE|INTERFACE|PROGRAM|TYPE|BLOCK|ASSOCIATE)\b',
                  r'END\1', text)
    # 5d. Strip bare ``IMPLICIT NONE`` (or any other IMPLICIT spec).
    #     Some S/D halves have it and others don't; it has no bearing
    #     on the migrated numerics.
    text = re.sub(r'^\s*IMPLICIT\s+\S.*$', '', text, flags=re.MULTILINE)
    # 6. Collapse runs of whitespace to a single space — BLAS sources
    #    hand-align declaration lists with varying spacing between
    #    type keyword and argument list, and indent DO loop bodies
    #    differently across precision halves.
    text = re.sub(r'[ \t]+', ' ', text)
    # 6b. Strip whitespace adjacent to arithmetic / relational
    #     operators. LAPACK's S and D sources format the same
    #     expression with different spacing (``NZ*EPS`` vs
    #     ``NZ* EPS``, ``I + 1`` vs ``I+1``). We squeeze both sides
    #     so they compare equal.
    text = re.sub(r'\s*([*+\-/=])\s*', r'\1', text)
    # 6c. Strip whitespace adjacent to parens. LAPACK formats
    #     argument lists with varying spacing (``F( X )`` vs
    #     ``F(X)``, ``IF (X)`` vs ``IF(X)``).
    text = re.sub(r'\s*\(\s*', '(', text)
    text = re.sub(r'\s*\)', ')', text)
    # 6d. Strip whitespace around commas. LAPACK formats complex
    #     literals and argument lists as ``(1.0, 0.0)`` vs
    #     ``(1.0,0.0)`` inconsistently between S/C and D/Z halves,
    #     and occasionally leaves a stray space before the comma
    #     (``1 ,WORK``).
    text = re.sub(r'\s*,\s*', ',', text)
    # 6e. Collapse doubled parens ``((X))`` → ``(X)``. LAPACK
    #     occasionally wraps a boolean / logical expression in an
    #     extra redundant paren on one half.
    prev = None
    while prev != text:
        prev = text
        text = re.sub(r'\(\(([^()]*)\)\)', r'(\1)', text)
    # 7. Sort items in declaration lists whose order isn't
    #    semantically meaningful. BLAS/LAPACK sources list
    #    intrinsics, externals, and type-declared variables in
    #    arbitrary (often precision-specific) order — e.g., S source
    #    writes "INTRINSIC REAL, CONJG, MAX" while Z source writes
    #    "INTRINSIC CONJG, MAX, MIN, REAL".
    def _sort_decl(match: re.Match) -> str:
        kw = match.group(1)
        items = [s.strip() for s in match.group(2).split(',')]
        return f'{kw} ' + ', '.join(sorted(items))
    text = re.sub(
        r'\b(INTRINSIC|EXTERNAL|INTEGER|LOGICAL|COMPLEX|REAL)\s+'
        r'([A-Za-z_@][A-Za-z0-9_@,\s]*?)(?=$)',
        _sort_decl, text, flags=re.MULTILINE,
    )
    return text


def _strip_c_comments(text: str) -> str:
    """Strip // and /* */ comments plus blank/whitespace-only lines."""
    # Remove /* ... */ block comments (possibly multi-line)
    text = re.sub(r'/\*.*?\*/', '', text, flags=re.DOTALL)
    # Remove // line comments
    text = re.sub(r'//[^\n]*', '', text)
    # Collapse blank/whitespace-only lines and trim trailing spaces
    lines = [ln.rstrip() for ln in text.split('\n')]
    lines = [ln for ln in lines if ln.strip()]
    return '\n'.join(lines)


def run_fortran_migration(config: RecipeConfig, rename_map: dict[str, str],
                          output_dir: Path, target_mode: TargetMode,
                          dry_run: bool = False,
                          classification=None,
                          parser: str | None = None,
                          parser_cmd: str | None = None) -> dict:
    """Run Fortran migration pipeline."""
    # Identify precision-independent symbols
    if classification is None:
        symbols = scan_symbols(config.source_dir, config.language,
                               config.extensions, config.library_path,
                               extra_c_return_types=tuple(config.c_return_types))
        classification = classify_symbols(symbols)
    independent = classification.independent

    migrated_count = 0
    copied_count = 0
    skipped: list[str] = []
    divergences: list[str] = []

    # Collect eligible source files. In addition to files with the
    # recipe's declared extensions, include any file whose stem is listed
    # in copy_files (regardless of extension) — this lets libraries
    # carry shared, type-independent headers (e.g. MUMPS's
    # mumps_tags.h / mumps_headers.h) through without having to widen
    # the extensions list and accidentally migrate them.
    #
    # extra_fortran_dirs contributes additional directories whose files
    # are migrated alongside source_dir. MUMPS uses this for
    # per-arithmetic headers (``dmumps_struc.h`` → ``qmumps_struc.h``)
    # that live in a parallel ``include/`` directory but are Fortran
    # content that must go through the full migration pipeline.
    src_dirs = [config.source_dir] + list(
        getattr(config, 'extra_fortran_dirs', []) or []
    )
    src_files = sorted(
        p for d in src_dirs for p in d.iterdir()
        if p.suffix.lower() in config.extensions
        or p.stem.upper() in config.copy_files
    )
    # extra_migrate_files pulls in specific leaf files from outside
    # source_dir / extra_fortran_dirs — used to target individual
    # helpers that live in a shared directory whose other contents
    # belong to a different library (LAPACK's INSTALL/dlamch.f,
    # PTZBLAS's TOOLS/zzdotc.f). The symbol scanner already sees these
    # via extra_symbol_dirs, so the rename map covers the referenced
    # names; this line only widens the emit loop to include them.
    extra_migrate = list(getattr(config, 'extra_migrate_files', []) or [])
    src_files = sorted(set(src_files).union(
        p for p in extra_migrate if p.is_file()
    ))
    # Swap any upstream source whose filename matches a recipe-level
    # override. The override file is in upstream shape (DOUBLE PRECISION,
    # pd*/dz* naming, etc.) and goes through the normal migration
    # pipeline so it produces correctly-renamed output for every target.
    if config.source_overrides:
        src_files = [
            (config.source_overrides[p.name]
             if p.name in config.source_overrides else p)
            for p in src_files
        ]

    # Convergence buffer: first writer of each output name stores its
    # text; subsequent writers must agree or we record a divergence.
    # D/Z sources are preferred as the canonical text: when a pair
    # (SGEMM, DGEMM) both target QGEMM, DGEMM's migrated body is kept
    # and SGEMM's is only consulted for the equality check. The
    # ``prefer_source`` recipe field flips this for individual stems
    # (e.g. to route around a bug that only exists in the D/Z half).
    prefer = config.prefer_source

    def _canonical_rank(stem: str) -> int:
        # 0 = highest priority (recipe-pinned), 1 = D/Z default,
        # 2 = S/C fallback. Sorting ascending makes the preferred
        # source the first writer for each target.
        u = stem.upper()
        if u in prefer:
            return 0
        if u and u[0] in ('D', 'Z'):
            return 1
        return 2

    src_files.sort(key=lambda p: (_canonical_rank(p.stem), p.name))

    canonical_text: dict[str, str] = {}
    canonical_normalized: dict[str, str] = {}
    canonical_source: dict[str, str] = {}

    # Build module-rename regex pairs once. Applied post-migration to
    # every migrated file (copy_files are deliberately untouched so the
    # verbatim upstream module keeps its original name).
    module_rename_pairs: list[tuple[re.Pattern[str], str]] = []
    for old_mod, new_mod in (config.module_renames or {}).items():
        pat = re.compile(r'(?i)(\bUSE\s+)' + re.escape(old_mod) + r'\b')
        module_rename_pairs.append((pat, r'\g<1>' + new_mod))

    def _apply_module_renames(text: str) -> str:
        for pat, repl in module_rename_pairs:
            text = pat.sub(repl, text)
        return text

    # Partition: copy/skip/dry-run decisions stay in the main process;
    # only the Flang-bound migration is dispatched to workers.
    to_migrate: list[Path] = []
    for src_path in src_files:
        stem_upper = src_path.stem.upper()
        if stem_upper in config.skip_files:
            skipped.append(src_path.name)
            continue
        is_copy = stem_upper in config.copy_files or stem_upper in independent
        if dry_run:
            if is_copy:
                tqdm.write(f'  {src_path.name} (copy)')
            else:
                out_name = target_filename(src_path.name, rename_map, target_mode)
                tqdm.write(f'  {src_path.name} → {out_name}')
            continue
        if is_copy:
            (output_dir / src_path.name).write_text(
                src_path.read_text(errors='replace'))
            copied_count += 1
            continue
        to_migrate.append(src_path)

    if dry_run:
        return {
            'migrated': 0, 'copied': copied_count,
            'skipped': skipped, 'divergences': divergences,
        }

    # Parallel migration. Each worker runs parser + regex substitution,
    # which is the dominant cost for large libraries like LAPACK.
    # We reduce results in canonical-first order so D/Z output is
    # the one written to disk.
    workers = max(1, (os.cpu_count() or 4))
    results: dict[Path, tuple[str, str] | None] = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = {
            ex.submit(migrate_file_to_string, p, rename_map, target_mode,
                      parser, parser_cmd,
                      config.keep_kind_lines.get(p.name)): p
            for p in to_migrate
        }
        for fut in tqdm(as_completed(futures), total=len(futures),
                        desc='  Migrating', unit='file',
                        mininterval=1.0, miniters=10):
            p = futures[fut]
            try:
                results[p] = fut.result()
            except Exception as exc:
                print(f'  warning: migration crashed on {p.name}: '
                      f'{type(exc).__name__}: {exc}', file=sys.stderr)
                results[p] = None

    # Reduce in deterministic D/Z-first order.
    for src_path in to_migrate:
        result = results.get(src_path)
        if result is None:
            skipped.append(src_path.name)
            continue
        out_name, migrated = result
        migrated = _apply_module_renames(migrated)
        normalized = _canonicalize_for_compare(
            _strip_fortran_comments(migrated, src_path.suffix)
        )

        prior = canonical_normalized.get(out_name)
        if prior is None:
            # First (canonical) writer for this target name
            (output_dir / out_name).write_text(migrated)
            canonical_text[out_name] = migrated
            canonical_normalized[out_name] = normalized
            canonical_source[out_name] = src_path.name
            migrated_count += 1
        elif prior == normalized:
            # Convergence: co-family member produced identical code
            # (comments may differ and are ignored).
            pass
        else:
            # Divergence: co-family members disagree. Keep the canonical
            # (D/Z) version on disk and record the mismatch.
            divergences.append(
                f'{src_path.name} vs {canonical_source[out_name]} '
                f'→ {out_name}'
            )

    return {
        'migrated': migrated_count,
        'copied': copied_count,
        'skipped': skipped,
        'divergences': divergences,
    }


def _scan_extra_dirs(extra_dirs: list[Path],
                     extra_c_return_types: tuple[str, ...]) -> set[str]:
    """Scan ``extra_symbol_dirs`` plus one level of subdirectories for
    Fortran and C symbols."""
    out: set[str] = set()
    for extra_dir in extra_dirs:
        if not extra_dir.is_dir():
            continue
        for lang, exts in [('fortran', ['.f', '.f90', '.F90']),
                           ('c', ['.c'])]:
            out |= scan_symbols(
                extra_dir, lang, exts,
                extra_c_return_types=extra_c_return_types)
        for sub in sorted(extra_dir.iterdir()):
            if sub.is_dir():
                for lang, exts in [('fortran', ['.f', '.f90', '.F90']),
                                   ('c', ['.c'])]:
                    out |= scan_symbols(
                        sub, lang, exts,
                        extra_c_return_types=extra_c_return_types)
    return out


def _collect_all_symbols(config: RecipeConfig,
                         project_root: Path | None) -> set[str]:
    """Scan the recipe plus every transitive dependency — including
    each dependency's ``extra_symbol_dirs`` — so the rename map
    covers kernels a dependency exposes indirectly (e.g. LAPACK's
    INSTALL/slamch.f reached via scalapack → lapack)."""
    own_c_types = tuple(config.c_return_types)
    symbols = scan_symbols(
        config.source_dir, config.language,
        config.extensions, config.library_path,
        extra_c_return_types=own_c_types,
    )
    symbols |= _scan_extra_dirs(config.extra_symbol_dirs, own_c_types)

    seen: set[Path] = set()
    queue = list(config.depends)
    while queue:
        dep_path = queue.pop(0)
        key = dep_path.resolve()
        if key in seen:
            continue
        seen.add(key)
        dep_cfg = load_recipe(dep_path, project_root)
        dep_c_types = tuple(dep_cfg.c_return_types)
        symbols |= scan_symbols(
            dep_cfg.source_dir, dep_cfg.language,
            dep_cfg.extensions, dep_cfg.library_path,
            extra_c_return_types=dep_c_types,
        )
        symbols |= _scan_extra_dirs(dep_cfg.extra_symbol_dirs, dep_c_types)
        queue.extend(dep_cfg.depends)
    return symbols


def run_divergence_report(recipe_path: Path, target_mode=None,
                           project_root: Path | None = None,
                           parser: str | None = None,
                           parser_cmd: str | None = None) -> list[dict]:
    """Migrate every co-family source pair in-memory and return the
    normalized diff for each pair whose members disagree.

    Returns a list of ``{'target', 'canonical', 'other', 'diff'}``
    dicts, sorted by target filename. ``diff`` is the list of +/-
    lines (without context) from the unified diff of the two
    canonicalized texts.
    """
    import difflib
    config = load_recipe(recipe_path, project_root)

    symbols = _collect_all_symbols(config, project_root)
    classification = classify_symbols(symbols)
    rename_map = classification.build_rename_map(target_mode)
    rename_map = _apply_extra_renames(rename_map, config, target_mode)

    # Group eligible source files by their target output name.
    src_files = sorted(
        p for p in config.source_dir.iterdir()
        if p.suffix.lower() in config.extensions
    )
    by_target: dict[str, list[Path]] = {}
    for p in src_files:
        stem_u = p.stem.upper()
        if stem_u in config.skip_files or stem_u in config.copy_files:
            continue
        if stem_u in classification.independent:
            continue
        by_target.setdefault(target_filename(p.name, rename_map, target_mode), []).append(p)

    # Migrate every member of every multi-member group in parallel.
    pairs: list[tuple[Path, Path]] = []
    for tn, members in by_target.items():
        if len(members) < 2:
            continue
        members.sort(key=lambda p: (
            0 if p.stem.upper() in config.prefer_source
            else (1 if p.stem and p.stem[0].upper() in ('D', 'Z') else 2),
            p.name,
        ))
        canonical = members[0]
        for other in members[1:]:
            pairs.append((canonical, other))

    all_paths = {p for pair in pairs for p in pair}
    texts: dict[Path, str] = {}
    workers = max(1, (os.cpu_count() or 4))
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = {
            ex.submit(migrate_file_to_string, p, rename_map, target_mode,
                      parser, parser_cmd,
                      config.keep_kind_lines.get(p.name)): p
            for p in all_paths
        }
        for fut in tqdm(as_completed(futures), total=len(futures),
                        desc='  Migrating', unit='file',
                        mininterval=1.0, miniters=10):
            p = futures[fut]
            try:
                res = fut.result()
                if res is None:
                    continue
                _, migrated = res
                texts[p] = migrated
            except Exception as exc:
                print(f'  warning: migration crashed on {p.name}: '
                      f'{type(exc).__name__}: {exc}', file=sys.stderr)

    report: list[dict] = []
    for canonical, other in pairs:
        if canonical not in texts or other not in texts:
            continue
        n_can = _canonicalize_for_compare(
            _strip_fortran_comments(texts[canonical], canonical.suffix))
        n_oth = _canonicalize_for_compare(
            _strip_fortran_comments(texts[other], other.suffix))
        if n_can == n_oth:
            continue
        diff = [
            line for line in difflib.unified_diff(
                n_oth.splitlines(), n_can.splitlines(),
                lineterm='', n=0,
            )
            if line.startswith(('-', '+'))
            and not line.startswith(('---', '+++'))
        ]
        report.append({
            'target': target_filename(canonical.name, rename_map, target_mode),
            'canonical': canonical.name,
            'other': other.name,
            'diff': diff,
        })
    report.sort(key=lambda r: (r['other'], r['canonical']))
    return report


def run_convergence_report(recipe_path: Path, output_dir: Path,
                            target_mode=None,
                            project_root: Path | None = None,
                            parser: str | None = None,
                            parser_cmd: str | None = None) -> list[dict]:
    """Verify that already-migrated files converge with a fresh
    migration of their S/C co-family siblings.

    Reads each pair's canonical (D/Z-derived) target from
    ``output_dir``, re-migrates the S/C member in memory, then
    compares the two with :func:`_light_normalize` — whitespace,
    case, and ``END KEYWORD`` merging only. All semantic
    transforms (identifier renames, literal kinds, casts, typedef
    substitutions) are the migrator's responsibility; this report
    flags anything that slipped through.

    Returns ``[{'target', 'canonical', 'other', 'diff', 'status'}]``
    where ``status`` is ``'diverged'`` for normalized mismatches and
    ``'missing'`` for pairs whose on-disk canonical is absent.
    """
    import difflib
    config = load_recipe(recipe_path, project_root)

    if config.language == 'c':
        return run_c_convergence_report(
            recipe_path, output_dir, target_mode, project_root
        )
    if config.language != 'fortran':
        print(f'  (converge currently only supports Fortran/C recipes; '
              f'{config.library} is {config.language} — skipping)')
        return []

    symbols = _collect_all_symbols(config, project_root)
    classification = classify_symbols(symbols)
    rename_map = classification.build_rename_map(target_mode)
    rename_map = _apply_extra_renames(rename_map, config, target_mode)

    src_files = sorted(
        p for p in config.source_dir.iterdir()
        if p.suffix.lower() in config.extensions
    )
    by_target: dict[str, list[Path]] = {}
    for p in src_files:
        stem_u = p.stem.upper()
        if stem_u in config.skip_files or stem_u in config.copy_files:
            continue
        if stem_u in classification.independent:
            continue
        by_target.setdefault(target_filename(p.name, rename_map, target_mode), []).append(p)

    pairs: list[tuple[Path, Path]] = []
    for tn, members in by_target.items():
        if len(members) < 2:
            continue
        members.sort(key=lambda p: (
            0 if p.stem.upper() in config.prefer_source
            else (1 if p.stem and p.stem[0].upper() in ('D', 'Z') else 2),
            p.name,
        ))
        canonical = members[0]
        for other in members[1:]:
            pairs.append((canonical, other))

    # Only the S/C (other) halves need re-migration; the D/Z
    # canonical is already on disk.
    others_to_migrate = {other for _, other in pairs}
    texts: dict[Path, str] = {}
    workers = max(1, (os.cpu_count() or 4))
    # Honor source_overrides during re-migration so that hand-edited
    # upstream replacements normalize both halves of a co-family pair,
    # not just the canonical D/Z half.
    def _resolve(p: Path) -> Path:
        ov = config.source_overrides.get(p.name)
        return ov if ov is not None else p
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = {
            ex.submit(migrate_file_to_string, _resolve(p), rename_map,
                      target_mode, parser, parser_cmd,
                      config.keep_kind_lines.get(p.name)): p
            for p in others_to_migrate
        }
        for fut in tqdm(as_completed(futures), total=len(futures),
                        desc='  Re-migrating S/C', unit='file',
                        mininterval=1.0, miniters=10):
            p = futures[fut]
            try:
                res = fut.result()
                if res is None:
                    continue
                _, migrated = res
                texts[p] = migrated
            except Exception as exc:
                print(f'  warning: re-migration crashed on {p.name}: '
                      f'{type(exc).__name__}: {exc}', file=sys.stderr)

    report: list[dict] = []
    for canonical, other in pairs:
        target_name = target_filename(canonical.name, rename_map, target_mode)
        on_disk = output_dir / target_name
        if not on_disk.is_file():
            report.append({
                'target': target_name,
                'canonical': canonical.name,
                'other': other.name,
                'diff': [],
                'status': 'missing',
            })
            continue
        if other not in texts:
            continue
        disk_text = on_disk.read_text(errors='replace')
        other_text = texts[other]
        # Apply any remaining recipe-level local renames (for non-precision
        # cases like CONJTOPH→TTOPH that aren't S↔D / C↔Z swaps).
        if config.local_renames:
            disk_text = _apply_local_renames(disk_text, config.local_renames)
            other_text = _apply_local_renames(other_text,
                                              config.local_renames)
        n_can = _light_normalize(
            _strip_fortran_comments(disk_text, canonical.suffix))
        n_oth = _light_normalize(
            _strip_fortran_comments(other_text, other.suffix))
        if n_can == n_oth:
            continue
        # Filter out precision-prefix local-variable drift (S↔D, C↔Z).
        # Any properly-classified symbol has the same migrated name in
        # both halves, so precision-prefix differences in the diff can
        # only come from locals not in the rename map.
        diff = _filter_precision_drift(
            n_oth.splitlines(), n_can.splitlines())
        if not diff:
            continue
        report.append({
            'target': target_name,
            'canonical': canonical.name,
            'other': other.name,
            'diff': diff,
            'status': 'diverged',
        })
    report.sort(key=lambda r: (r['status'], r['other'], r['canonical']))
    return report


def run_c_convergence_report(recipe_path: Path, output_dir: Path,
                              target_mode=None,
                              project_root: Path | None = None) -> list[dict]:
    """Verify that migrated C files on disk converge with a fresh
    migration of their S/C co-family siblings.

    C counterpart of :func:`run_convergence_report`. Each eligible
    source ``.c`` file in ``config.source_dir`` is migrated in memory
    via :func:`migrate_c_file_to_string` to discover its target name;
    files that collide on a target form co-family groups. The D/Z
    member is the canonical (read from disk at ``output_dir``); each
    S/C sibling's in-memory migration result is compared with
    :func:`_light_normalize_c` after stripping C comments.

    Supports both scalapack-style recipes (rename-map driven, e.g.
    PBLAS) and direct-style recipes (prefix driven, e.g. BLACS).
    """
    import difflib
    config = load_recipe(recipe_path, project_root)

    if config.language != 'c':
        print(f'  (run_c_convergence_report called on non-C recipe '
              f'{config.library} — skipping)')
        return []

    rename_map: dict[str, str] | None = None
    classification = None
    # Every C recipe except BLACS uses the rename-map-driven path —
    # see the matching gate in run_migration().
    if config.library != 'blacs':
        symbols = _collect_all_symbols(config, project_root)
        classification = classify_symbols(symbols)
        rename_map = classification.build_rename_map(target_mode)
        rename_map = _apply_extra_renames(rename_map, config, target_mode)

    src_files = sorted(
        p for p in config.source_dir.iterdir()
        if p.suffix.lower() == '.c'
    )

    # Migrate every source in memory and group by target filename. For
    # BLACS the target-computing cost is trivial; for PBLAS we'd save
    # work by skipping D/Z here, but the rename_map regex build is the
    # hot part and unavoidable either way.
    by_target: dict[str, list[tuple[Path, str]]] = {}
    for p in tqdm(src_files, desc='  Re-migrating C sources', unit='file',
                  mininterval=1.0, miniters=10):
        try:
            res = migrate_c_file_to_string(
                p, target_mode,
                rename_map=rename_map,
                classification=classification,
                c_type_aliases=config.c_type_aliases,
            )
        except Exception as exc:
            print(f'  warning: C re-migration crashed on {p.name}: '
                  f'{type(exc).__name__}: {exc}', file=sys.stderr)
            continue
        if res is None:
            continue
        target_name, text = res
        by_target.setdefault(target_name, []).append((p, text))

    def _is_dz(src_path: Path) -> bool:
        stem = src_path.stem
        base = stem[3:] if stem.startswith('BI_') else stem
        return bool(base) and base[0].lower() in ('d', 'z')

    report: list[dict] = []
    for target_name, members in by_target.items():
        if len(members) < 2:
            continue
        members.sort(key=lambda m: (not _is_dz(m[0]), m[0].name))
        canonical_path, _ = members[0]
        on_disk = output_dir / target_name
        if not on_disk.is_file():
            for other_path, _ in members[1:]:
                report.append({
                    'target': target_name,
                    'canonical': canonical_path.name,
                    'other': other_path.name,
                    'diff': [],
                    'status': 'missing',
                })
            continue
        disk_text = on_disk.read_text(errors='replace')
        n_can = _light_normalize_c(_strip_c_comments(disk_text))
        for other_path, other_text in members[1:]:
            n_oth = _light_normalize_c(_strip_c_comments(other_text))
            if n_can == n_oth:
                continue
            diff = [
                line for line in difflib.unified_diff(
                    n_oth.splitlines(), n_can.splitlines(),
                    lineterm='', n=0,
                )
                if line.startswith(('-', '+'))
                and not line.startswith(('---', '+++'))
            ]
            report.append({
                'target': target_name,
                'canonical': canonical_path.name,
                'other': other_path.name,
                'diff': diff,
                'status': 'diverged',
            })
    report.sort(key=lambda r: (r['status'], r['other'], r['canonical']))
    return report


def run_c_migration(config: RecipeConfig, output_dir: Path,
                    target_mode: TargetMode, dry_run: bool = False,
                    classification=None,
                    rename_map: dict[str, str] | None = None) -> dict:
    """Run C migration pipeline (clone-and-substitute).

    When `classification` and `rename_map` are supplied, the generic
    rename-map-driven migration is used (for ScaLAPACK-style libraries
    like PBLAS). Otherwise the BLACS-specific path is used.
    """
    if dry_run:
        print('  (dry-run for C migration not yet implemented)')
        return {'cloned': [], 'template_vars': {}}

    overrides = _resolve_overrides(config, target_mode)

    result = migrate_c_directory(
        config.source_dir, output_dir, target_mode,
        copy_originals=config.copy_all_originals,
        classification=classification,
        rename_map=rename_map,
        c_type_aliases=config.c_type_aliases,
        c_pointer_cast_aliases=config.c_pointer_cast_aliases,
        header_patches=config.header_patches,
        overrides=overrides,
        extra_c_dirs=config.extra_c_dirs,
        skip_files=config.skip_files,
        copy_files=config.copy_files,
        source_overrides=config.source_overrides,
    )
    return result


def _apply_fortran_overrides(config: RecipeConfig,
                             target_mode: TargetMode,
                             output_dir: Path) -> None:
    """Copy target-gated override files into the Fortran output dir.

    Used for hand-written shared modules that cannot be produced by
    migration (e.g. MUMPS's extended-precision reallocator module).
    """
    overrides = _resolve_overrides(config, target_mode)
    for src_path, dst_name in overrides:
        (output_dir / dst_name).write_text(src_path.read_text())
        print(f'  Override:  {dst_name} (from {src_path.name})')


def _resolve_overrides(config: RecipeConfig,
                       target_mode: TargetMode) -> list[tuple[Path, str]]:
    """Select the active target's override entries from the recipe.

    The ``overrides`` recipe field is a dict keyed by target name. For
    the currently active target, return a list of ``(src_path, dst_name)``
    pairs with ``src_path`` resolved against the recipe directory.
    Returns an empty list if the recipe has no overrides for this target.
    """
    target_entry = (config.overrides or {}).get(target_mode.name)
    if not target_entry:
        return []
    src_dir_rel = target_entry.get('src_dir', '')
    files = target_entry.get('files', [])
    if config.recipe_dir is None:
        raise RuntimeError(
            'RecipeConfig.recipe_dir is not set; cannot resolve overrides'
        )
    src_dir = (config.recipe_dir / src_dir_rel).resolve()
    return [(src_dir / fname, fname) for fname in files]


def run_migration(recipe_path: Path, output_dir: Path,
                  target_mode=None, dry_run: bool = False,
                  project_root: Path | None = None,
                  parser: str | None = None,
                  parser_cmd: str | None = None) -> dict:
    """Run the full migration pipeline for a library.

    Args:
        recipe_path: Path to the YAML recipe file.
        output_dir: Where to write migrated files.
        target_mode: 10 or 16.
        dry_run: If True, show what would be done without writing.
        project_root: Project root for resolving relative paths.
        parser: Parse tree backend: ``'flang'``, ``'gfortran'``, or
            ``None`` (regex-only).
        parser_cmd: Explicit path to the compiler binary.

    Returns:
        Summary dict with migration statistics.
    """
    config = load_recipe(recipe_path, project_root)
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f'Library:     {config.library}')
    print(f'Language:    {config.language}')
    print(f'Source:      {config.source_dir}')
    print(f'Target:      {target_mode.name}')
    print(f'Prefix:      {config.prefix_style}')
    print()

    # Scan symbols and classify precision families
    print('Scanning symbols...')
    own_symbols = scan_symbols(
        config.source_dir, config.language,
        config.extensions, config.library_path,
        extra_c_return_types=tuple(config.c_return_types),
    )
    own_count = len(own_symbols)
    symbols = _collect_all_symbols(config, project_root)

    classification = classify_symbols(symbols)
    rename_map = classification.build_rename_map(target_mode)
    rename_map = _apply_extra_renames(rename_map, config, target_mode)

    # NOTE: target_mode.known_constants (ZERO/ONE/...) are handled
    # per-file by strip_known_constants_from_decls +
    # replace_known_constants. We intentionally do NOT add them to the
    # global rename_map: doing so would also rewrite the LHS aliases
    # of LAPACK ``USE LA_CONSTANTS, ONLY: zero=>dzero`` clauses (and
    # rename the lowercase ``zero`` reference everywhere) which breaks
    # the alias binding.


    print(f'  {own_count} own symbols + {len(symbols) - own_count} from dependencies')
    print(f'  {len(classification.families)} precision families')
    print(f'  {len(classification.independent)} independent symbols')
    print(f'  {len(rename_map)} renames computed')

    # Dispatch to language-specific migrator
    print(f'\nMigrating to {target_mode.name}...')

    if config.language == 'fortran':
        result = run_fortran_migration(
            config, rename_map, output_dir, target_mode, dry_run,
            classification=classification,
            parser=parser, parser_cmd=parser_cmd,
        )
        if not dry_run:
            _apply_fortran_overrides(config, target_mode, output_dir)
        if not dry_run:
            print(f'\n  Migrated: {result["migrated"]} files')
            print(f'  Copied:   {result["copied"]} files (precision-independent)')
            if result['skipped']:
                print(f'  Skipped:  {len(result["skipped"])} files')
            if result.get('divergences'):
                print(f'  Divergences: {len(result["divergences"])} '
                      f'(co-family members produced differing output)')
                for d in result['divergences'][:10]:
                    print(f'    {d}')
                if len(result['divergences']) > 10:
                    print(f'    ... ({len(result["divergences"]) - 10} more)')
    elif config.language == 'c':
        # BLACS keeps the legacy hardcoded-pattern migrator (it carries
        # MPI typedef patches, Bdef.h rewrites, MPI_REAL16 check
        # generation, and BLACS-specific Cd*/BI_d* routine patterns
        # that have no analogue in other C libraries). Every other C
        # recipe (PBLAS / ScaLAPACK_C / XBLAS / future libraries) uses
        # the generic rename-map-driven cloner — the prefix classifier
        # discovers slot positions empirically, so the cloner is
        # naming-convention-agnostic.
        if config.library == 'blacs':
            result = run_c_migration(config, output_dir, target_mode, dry_run)
        else:
            result = run_c_migration(
                config, output_dir, target_mode, dry_run,
                classification=classification, rename_map=rename_map,
            )
        if not dry_run:
            print(f'\n  Cloned: {len(result["cloned"])} files')
            if result.get('divergences'):
                print(f'  Divergences: {len(result["divergences"])} '
                      f'(co-family members produced differing output)')
                for d in result['divergences'][:10]:
                    print(f'    {d}')
                if len(result['divergences']) > 10:
                    print(f'    ... ({len(result["divergences"]) - 10} more)')
    else:
        raise ValueError(f'Unsupported language: {config.language}')

    return result
