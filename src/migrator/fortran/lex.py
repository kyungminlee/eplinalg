"""Lexical + scanning primitives for the Fortran migrator (bottom layer).

Form-agnostic string/line tools and the shared token regexes that every
per-line rewriter builds on. Depends on nothing else in the package.
Extracted verbatim from ``fortran_migrator.py``.
"""
import re


# Hoisted patterns reused across hot-path callers (replace_literals,
# replace_generic_conversions, _rewrite_int_of_complex / _rewrite_int_kind_on_real64x2,
# _filter_known_constants_from_decl, _unwrap_redundant_constructors). All
# depend only on constants known at module load time or values that are
# stable for the duration of a run; rebuilding them per line was the
# dominant cost in microbenchmarks of the migrator hot loops.
#
# ``.EQ.`` / ``.AND.`` etc. operator detector used in ``replace_literals``
# to mask operator dots out before the bare-literal pass scans for
# floating-point suffixes.
_FORTRAN_OP_RE = re.compile(
    r'\.\s*(EQ|NE|LT|GT|LE|GE|AND|OR|NOT|TRUE|FALSE|EQV|NEQV)\s*\.',
    re.IGNORECASE,
)


# ``INT(`` / ``NINT(`` detectors for the two _rewrite_int_* helpers.
_INT_CALL_RE = re.compile(r'\bINT\s*\(', re.IGNORECASE)


_NINT_CALL_RE = re.compile(r'\bNINT\s*\(', re.IGNORECASE)


_FP_VALUE_RE = re.compile(
    r'^[+-]?\s*(?:'
    r'\d+\.\d*|\d*\.\d+|\d+'
    r')(?:[DdEe][+-]?\d+)?\s*$'
)


def _is_fp_value(text: str, known_constants: dict) -> bool:
    """Heuristic: does the trimmed PARAMETER/DATA value look like an
    FP literal, a known-constant reference, or an expression that
    contains one?

    Tighter than a substring scan: rejects identifiers like ``ELEMENT``
    that happen to contain ``E``. The composite-expression case splits
    on Fortran operators and checks each operand: if ANY operand is a
    known FP constant or an FP literal, the whole expression counts.
    """
    s = text.strip()
    if not s:
        return False
    if _FP_VALUE_RE.match(s):
        # Pure numeric token: FP iff it has a decimal point or exponent.
        return ('.' in s) or ('D' in s.upper()) or ('E' in s.upper())
    # Tokenize on identifiers + numeric literals; ignore operators.
    tokens = re.findall(r'[A-Za-z_]\w*|\d+\.\d*[DdEe][+-]?\d+|\d*\.\d+[DdEe]?[+-]?\d*|\d+\.\d*|\d*\.\d+', s)
    for tok in tokens:
        if tok.upper() in known_constants:
            return True
        if _FP_VALUE_RE.match(tok) and (('.' in tok) or ('D' in tok.upper()) or ('E' in tok.upper())):
            return True
    return False


def _join_continued_lines(lines: list[str], start: int) -> tuple[str, int]:
    """Join a fixed-form statement starting at ``lines[start]`` with any
    continuation lines (column-6 marker). Returns ``(joined_text, end)``
    where ``end`` is the index of the first line NOT consumed.
    """
    out = lines[start].rstrip('\n')
    j = start + 1
    while j < len(lines):
        nxt = lines[j]
        if is_continuation_line(nxt):
            out += ' ' + nxt[6:].rstrip('\n')
            j += 1
            continue
        break
    return out, j


_PROC_HEADER_RE_SCOPE = re.compile(
    r'^\s*(?:RECURSIVE\s+|PURE\s+|ELEMENTAL\s+)*'
    r'(?:(?:INTEGER|REAL|COMPLEX|LOGICAL|CHARACTER|TYPE\s*\([^)]+\)|DOUBLE\s+PRECISION|DOUBLE\s+COMPLEX)'
    r'(?:\s*\*\s*\d+)?\s+)?'
    r'(?:PROGRAM|SUBROUTINE|FUNCTION|MODULE|BLOCK\s+DATA)\b',
    re.IGNORECASE,
)


_INTERFACE_BEGIN_RE = re.compile(r'^\s*(?:ABSTRACT\s+)?INTERFACE\b', re.IGNORECASE)


_INTERFACE_END_RE = re.compile(r'^\s*END\s*INTERFACE\b', re.IGNORECASE)


def _scope_indices(lines: list[str]) -> list[int]:
    """Return a per-line list mapping line index → procedure scope index.

    Scope index 0 is the first SUBROUTINE/FUNCTION header encountered
    when scanning forward from the top of the file. Lines before any
    header are scope -1 (module/global level). Used by
    ``convert_parameter_stmts`` / ``convert_data_stmts`` to tag each
    converted assignment with the scope it belongs to, so that
    ``insert_use_multifloats`` can insert only the assignments
    belonging to the current scope.

    INTERFACE-block inner SUBROUTINE/FUNCTION declarations are NOT
    counted as new scopes — they declare prototypes for external
    procedures, not local scopes that take runtime assignments.

    Replaces the prior per-call ``_scope_index_at`` helper that
    rescanned ``lines[0..i]`` on every invocation (O(N²) when called
    from a per-statement outer loop). Compute the full vector once
    in O(N) and index into it.
    """
    scopes: list[int] = []
    scope = -1
    in_interface = 0
    for ln in lines:
        if _INTERFACE_BEGIN_RE.match(ln):
            in_interface += 1
        elif _INTERFACE_END_RE.match(ln):
            if in_interface > 0:
                in_interface -= 1
        elif in_interface == 0 and _PROC_HEADER_RE_SCOPE.match(ln):
            scope += 1
        scopes.append(scope)
    return scopes


_STMT_FN_RE = re.compile(r'^[A-Za-z_]\w*\s*\(\s*[A-Za-z_]\w*\s*\)\s*=')


def _looks_like_statement_function(stripped: str, lines: list[str], k: int) -> bool:
    """Heuristic: does line k look like a LAPACK statement function
    definition that the declaration-block walker should step over?

    A statement function has the form ``NAME(SCALAR_ARG) = expression``
    and appears in the LAPACK source between a comment marker like
    ``*     .. Statement Function definitions ..`` and the executable
    statements section. We detect by looking back at recent lines for
    that marker (within ~10 lines).
    """
    if not _STMT_FN_RE.match(stripped):
        return False
    look = max(0, k - 12)
    for kk in range(k - 1, look - 1, -1):
        prev = lines[kk]
        if not prev.strip():
            continue
        if 'Statement Function' in prev:
            return True
        if prev and prev[0] in ('C', 'c', '*', '!'):
            continue
        if prev.lstrip().startswith('!'):
            continue
        # Hit a code line that wasn't a statement function — stop.
        if not _STMT_FN_RE.match(prev.lstrip()):
            return False
    return False


# ``__KEEPKIND_DP__`` is the keep-kind sentinel that stands in for
# ``DOUBLE PRECISION`` (see fortran/keepkind.py). ``specialize_use_module``
# scans for locally-declared names while the sentinel is still in place —
# the restore happens only after migrate_*_form returns — so this pattern
# must recognise the sentinel too, or a keep-kind-protected local (e.g. a
# ``DOUBLE PRECISION :: GAMMA`` that collides with a multifloats generic)
# is missed and wrongly imported via ``USE ..., only:``.
_DECL_LINE_RE = re.compile(
    r'^\s+(?:TYPE\s*\([^)]*\)|INTEGER\b|REAL\b|COMPLEX\b|LOGICAL\b|'
    r'CHARACTER\b|DOUBLE\s+PRECISION\b|DOUBLE\s+COMPLEX\b|__KEEPKIND_DP__)',
    re.IGNORECASE,
)


_IDENT_RE = re.compile(r"\b([A-Za-z_]\w*)\b")


_STRING_RE = re.compile(r"'(?:[^']|'')*'|\"(?:[^\"]|\"\")*\"")


def _strip_strings_and_comments(line: str) -> str:
    """Drop string literals and trailing inline comments."""
    out = _STRING_RE.sub('', line)
    bang = out.find('!')
    if bang >= 0:
        out = out[:bang]
    return out


def _scan_local_declared_names(proc_lines: list[str]) -> set[str]:
    """Collect local variable names from type-declaration statements
    inside a procedure. Used to suppress matching multifloats public
    names so that the local variable can shadow the use-associated
    generic interface (gfortran refuses if the name is in scope).
    """
    names: set[str] = set()
    joined: list[str] = []
    cur = ''
    for raw in proc_lines:
        # Skip pure comment lines (fixed and free form).
        if raw[:1] in ('C', 'c', '*', '!'):
            continue
        if cur and is_continuation_line(raw):
            # Fixed-form continuation.
            cur += ' ' + raw[6:].rstrip('\n')
            continue
        if cur and cur.rstrip().endswith('&'):
            cur = cur.rstrip().rstrip('&') + ' ' + raw.lstrip().lstrip('&').rstrip('\n')
            continue
        if cur:
            joined.append(cur)
        cur = raw.rstrip('\n')
    if cur:
        joined.append(cur)

    for stmt in joined:
        if not _DECL_LINE_RE.match(stmt):
            continue
        # Strip the type prefix (everything up to and including the
        # first ``::`` if present, otherwise up to the type keyword).
        if '::' in stmt:
            tail = stmt.split('::', 1)[1]
        else:
            m = _DECL_LINE_RE.match(stmt)
            tail = stmt[m.end():]
        tail = _strip_strings_and_comments(tail)
        # Drop array specs and KIND parameters in parentheses so we
        # don't pick up bound expressions (``WNRM(MAX(M,N))``) as
        # local-variable names. Iterate until no more nested parens
        # remain.
        while True:
            new_tail = re.sub(r'\([^()]*\)', '', tail)
            if new_tail == tail:
                break
            tail = new_tail
        # Drop ``= initializer`` clauses (PARAMETER initializers may
        # contain references to module names like ``complex128x2``
        # that should not be treated as locally-declared variables).
        # Split at top-level commas, then drop everything from ``=``
        # onward in each item.
        items = []
        cur = ''
        for ch in tail + ',':
            if ch == ',':
                items.append(cur)
                cur = ''
            else:
                cur += ch
        for item in items:
            lhs = item.split('=', 1)[0]
            for m in _IDENT_RE.finditer(lhs):
                names.add(m.group(1).lower())
    return names


def _scan_referenced_identifiers(proc_lines: list[str]) -> set[str]:
    """Lower-cased identifiers referenced anywhere in the procedure
    body, excluding comments and string literals."""
    names: set[str] = set()
    for raw in proc_lines:
        if raw[:1] in ('C', 'c', '*', '!'):
            continue
        cleaned = _strip_strings_and_comments(raw)
        for m in _IDENT_RE.finditer(cleaned):
            names.add(m.group(1).lower())
    return names


# Procedure header for SUBROUTINE/FUNCTION/PROGRAM/MODULE/BLOCK DATA.
# ``MODULE PROCEDURE`` (inside an INTERFACE block) is NOT a procedure
# header and is excluded — injecting USE between ``MODULE PROCEDURE foo``
# and ``END INTERFACE`` is illegal Fortran.
_PROC_HEADER_RE = re.compile(
    r'^(\s{6,}|^\s*)(?:RECURSIVE\s+|PURE\s+|ELEMENTAL\s+)*'
    r'(?:(?:INTEGER|REAL|COMPLEX|LOGICAL|CHARACTER|TYPE\s*\([^)]+\)|DOUBLE\s+PRECISION|DOUBLE\s+COMPLEX)'
    r'(?:\s*\*\s*\d+)?\s+)?'
    r'(?:PROGRAM|SUBROUTINE|FUNCTION|MODULE(?!\s+PROCEDURE\b)|BLOCK\s+DATA)\b',
    re.IGNORECASE,
)


# ``END SUBROUTINE FOO``, ``END FUNCTION``, plain ``END``. Crucially
# NOT ``END IF`` / ``END DO`` / ``END SELECT`` — the keyword whitelist
# must be required when any word follows ``END``, otherwise inner
# control-flow ENDs would falsely terminate body scans.
_END_PROC_RE = re.compile(
    r'^\s*END(?:\s+(?:PROGRAM|SUBROUTINE|FUNCTION|MODULE|BLOCK\s*DATA)'
    r'(?:\s+\w+)?)?\s*(?:!.*)?$',
    re.IGNORECASE,
)


def is_comment_line(line: str) -> bool:
    return bool(line) and line[0] in ('C', 'c', '*', '!')


def _iter_outside_strings(text: str):
    """Yield ``(i, ch)`` for each character in ``text`` that is NOT
    inside a Fortran string literal. Recognizes both ``'`` and ``"``
    delimiters and the Fortran doubled-quote escape (``''`` inside a
    ``'``-string, ``""`` inside a ``"``-string).

    Canonical primitive shared by every scanner that needs to find
    inline ``!`` comments, count parens, or build a string-mask while
    respecting quoted-literal boundaries.
    """
    n = len(text)
    in_string = False
    quote = ''
    i = 0
    while i < n:
        ch = text[i]
        if in_string:
            if ch == quote:
                if i + 1 < n and text[i + 1] == quote:
                    # Doubled-quote escape — both chars stay in-string.
                    i += 2
                    continue
                in_string = False
        else:
            if ch in ("'", '"'):
                in_string = True
                quote = ch
            else:
                yield i, ch
        i += 1


def _ends_in_string(text: str) -> bool:
    """True if a character-literal is still open at the end of ``text``.

    Recognizes both ``'`` and ``"`` delimiters and the Fortran
    doubled-quote escape. Used to decide whether a continuation join sits
    inside a string literal (join tight, no inserted blank) or at a token
    boundary (a blank stands in for the break).
    """
    in_string = False
    quote = ''
    n = len(text)
    i = 0
    while i < n:
        ch = text[i]
        if in_string:
            if ch == quote:
                if i + 1 < n and text[i + 1] == quote:
                    i += 2
                    continue
                in_string = False
        elif ch in ("'", '"'):
            in_string = True
            quote = ch
        i += 1
    return in_string


def _find_inline_bang(text: str) -> int:
    """Return the index of the first inline ``!`` comment marker in
    ``text``, or ``len(text)`` if none. Quote-aware (doubled-quote
    escape included).

    Fast path: when the line has no string-delimiter characters at
    all (the common case for code-bearing source lines), fall through
    to ``str.find`` which runs in C. Only when a quote appears do we
    need the full quote-aware scan.
    """
    if "'" not in text and '"' not in text:
        idx = text.find('!')
        return idx if idx != -1 else len(text)
    for i, ch in _iter_outside_strings(text):
        if ch == '!':
            return i
    return len(text)


def _count_open_parens(line: str) -> int:
    """Net paren delta for ``line`` ignoring quoted strings and inline
    fixed-form ``!`` / ``C`` / ``*`` comments. Used to track when a
    SUBROUTINE/FUNCTION formal-arg list is still open across CPP
    ``#if/#endif`` blocks."""
    if not line:
        return 0
    if line[0] in ('C', 'c', '*'):
        return 0
    depth = 0
    for _, ch in _iter_outside_strings(line):
        if ch == '!':
            break
        if ch == '(':
            depth += 1
        elif ch == ')':
            depth -= 1
    return depth


def is_continuation_line(line: str) -> bool:
    # Fixed-form continuation: cols 1-5 blank (spaces only — NOT tabs), col 6
    # holds a non-blank, non-zero marker. A tab in col 1 is the gfortran/intel
    # extension that means "rest of line is normal source from col 7"; such a
    # line is a fresh statement, not a continuation, even when col 6 (i.e. the
    # 6th character after the tab) happens to be alphabetic. See dlarre2.f in
    # ScaLAPACK 2.2.3 for the wild form: `\t    CALL DLARRC(...)`.
    if not line or line[0] == '\t':
        return False
    return len(line) > 5 and line[0:5] == '     ' and line[5] not in (' ', '0', '')


def _build_split_mask(body: str) -> list[bool]:
    """Boolean mask over ``body`` — True at positions safe to split a
    fixed-form line at (i.e. not inside a string literal). The opening
    quote of a literal is itself marked unsafe, matching legacy
    behavior so the splitter never lands a continuation marker
    immediately before a string."""
    mask = [False] * len(body)
    for i, _ in _iter_outside_strings(body):
        mask[i] = True
    return mask
