"""Floating-point literal + INT/NINT kind rewriting (Cluster B).

Retypes floating-point literals to the target kind (honouring keep-kind
sentinels) and fixes INT()/NINT() calls on migrated complex/real64x2 operands.
Extracted verbatim from ``fortran_migrator.py``.
"""
import re

from ..target_mode import TargetMode
from .lex import _FORTRAN_OP_RE, _INT_CALL_RE, _NINT_CALL_RE
from .keepkind import _KK_SENTINEL, _KK_DBLE_SENTINEL, _KK_DCMPLX_SENTINEL


def replace_literals(line: str, target_mode: TargetMode,
                     source_kind: int | None = None) -> str:
    """Replace floating-point literals with target form.

    KIND mode: 1.0D+0 → 1.0E+0_k, 0.0E+0 → 0.0E+0_k, 0.0 → 0.0E0_k.
    Constructor mode: 1.0D+0 → float64x2('1.0D+0') or float64x2(1.0D0).

    Bare unsuffixed literals (no D/E exponent, no ``_kind``) are also
    promoted so that complex constants like ``(0.0,0.0)`` in C sources
    converge with ``(0.0d0,0.0d0)`` in Z sources after migration.

    Rule (a) gating by ``source_kind`` (mirrors :func:`replace_type_decls`):
    a kind4 source half (S/C) only promotes kind4-shaped literals
    (``1.0E0``, bare ``1.0``); a kind8 source half (D/Z) only promotes
    kind8-shaped literals (``1.0D0``). Cross-kind literals inside a
    half are intentional explicit-kind references and must survive
    migration. ``None`` falls back to the legacy promote-everything
    behavior.
    """
    # Keep-kind preserves the LHS type (DOUBLE PRECISION) or call-site
    # wrapper (dble/dcmplx) on this line. In constructor mode, wrapping
    # an RHS literal as ``real64x2(limbs=[...])`` produces a derived-
    # type value the keep-kind LHS / wrapper cannot absorb (no implicit
    # narrowing exists for derived types). Skip the literal pass on
    # such lines and leave the original DP literal in place — it
    # matches the preserved LHS / wrapper type. The kind_suffix path
    # is unaffected because REAL(KIND=N) → REAL(8) narrowing IS legal.
    if (target_mode.literal_mode == 'constructor'
            and (_KK_SENTINEL in line
                 or _KK_DBLE_SENTINEL in line
                 or _KK_DCMPLX_SENTINEL in line)):
        return line

    # ``INTEGER, PARAMETER :: name = <fp_literal>`` lines exist in
    # upstream MUMPS as a degenerate idiom (e.g. dsol_fwd_aux.F:1093
    # ``INTEGER, PARAMETER :: ZERO = 0.0D0``) — gfortran tolerates the
    # FP literal via integer coercion. In multifloats mode, wrapping
    # the literal to ``real64x2(limbs=[...])`` produces a derived-type
    # value that the INTEGER PARAMETER cannot accept. Leave the
    # original FP literal alone on those lines.
    if (target_mode.literal_mode == 'constructor'
            and re.match(r'^\s*INTEGER\b.*\bPARAMETER\b.*::',
                         line, re.IGNORECASE)):
        return line

    def literal_sub(m):
        mantissa = m.group(1)
        suffix = m.group(2)
        exp_rest = m.group(3)
        if exp_rest.startswith('+'):
            exp_rest = exp_rest[1:]

        # Rule (a): the suffix encodes the literal's source kind
        # (E = kind4, D = kind8). Skip the rewrite if the source half
        # doesn't own that kind.
        suffix_kind = 8 if suffix.upper() == 'D' else 4
        if source_kind is not None and source_kind != suffix_kind:
            return m.group(0)

        if target_mode.literal_mode == 'kind_suffix':
            return f'{mantissa}E{exp_rest}_{target_mode.kind_suffix}'
        else:
            # Named-component (``limbs=...``) structure-constructor
            # form. This is a constant expression and is therefore
            # legal in PARAMETER initializers — the alternative
            # ``float64x2('1.0D0')`` is a function call (the
            # ``mf_from_char`` generic) and would be rejected by the
            # compiler in PARAMETER context.
            return (
                f"{target_mode.real_constructor}"
                f"(limbs=[{mantissa}D{exp_rest}, 0.0_8])"
            )

    parts = re.split(r"('(?:[^']|'')*'|\"(?:[^\"]|\"\")*\")", line)
    for idx in range(0, len(parts), 2):
        seg = parts[idx]
        masked = _FORTRAN_OP_RE.sub(
            lambda m: '\x00' + m.group(1) + '\x00', seg,
        )

        # ``1.23_wp`` style literals are normalized first into a
        # placeholder so that the subsequent passes (D/E exponent and
        # bare-literal substitution) do not see the substituted form
        # and re-wrap its inner ``1.23D0`` text. We restore the
        # placeholders at the very end of the loop.
        wp_placeholders: list[str] = []

        def _wp_sub(m):
            val = m.group(1)
            if 'D' not in val.upper() and 'E' not in val.upper():
                if '.' in val:
                    val += 'D0'
                else:
                    val += '.0D0'
            tok = f"\x01WP{len(wp_placeholders)}\x01"
            wp_placeholders.append(
                f"{target_mode.real_constructor}(limbs=[{val}, 0.0_8])"
            )
            return tok

        if target_mode.literal_mode == 'constructor':
            masked = re.sub(r'(\d+\.\d*|\d*\.\d+)_wp', _wp_sub, masked, flags=re.IGNORECASE)

        # The negative lookbehind also rejects ``[`` (so we don't
        # match the inner ``.0D0`` of a previously-wrapped
        # ``float64x2(limbs=[1.0D0, 0.0_8])`` form) and rejects a
        # leading digit (so we don't match the ``.0D0`` substring of
        # a complete ``1.0D0`` literal that's already been consumed).
        masked = re.sub(
            r'(?<![\[\d])(\d+\.\d*|\d*\.\d+)([DEde])([+-]?\d+)',
            literal_sub, masked,
        )
        # Bare unsuffixed literals (e.g. ``1.0``, ``0.5``) are kind4
        # by Fortran default. Skip when ``source_kind == 8`` so a
        # kind8 source half preserves them as kind4 (rule a).
        if source_kind != 8:
            if target_mode.literal_mode == 'kind_suffix':
                masked = re.sub(
                    r'(?<![.\w])(\d+\.\d*|\d*\.\d+)(?![DdEe\w]|_\d)',
                    rf'\1E0_{target_mode.kind_suffix}', masked,
                )
            else:
                def bare_sub(m):
                    val = m.group(1)
                    return (
                        f"{target_mode.real_constructor}"
                        f"(limbs=[{val}D0, 0.0_8])"
                    )
                # The negative lookbehind on ``[`` skips literals that
                # already live inside a previously-wrapped
                # ``float64x2(limbs=[...])`` form.
                masked = re.sub(
                    r'(?<![.\w\[])(\d+\.\d*|\d*\.\d+)(?![DdEe\w]|_\d)',
                    bare_sub, masked,
                )
            
        if target_mode.literal_mode == 'constructor':
            # Wrap complex constants (float64x2(...), float64x2(...)) in complex128x2(...).
            # An optional unary +/- is allowed before each component to
            # cover cases like ``(-1.0D0, 0.0D0)`` from LAPACK.
            # Bare-name complex literals like ``(MF_ZERO, MF_ZERO)`` are
            # wrapped later in the per-line pipeline by
            # _wrap_bare_complex_literals (after replace_known_constants
            # has produced the MF_* names).
            masked = re.sub(
                r'\(\s*([-+]?\s*' + target_mode.real_constructor + r'\([^)]+\))\s*,\s*([-+]?\s*' + target_mode.real_constructor + r'\([^)]+\))\s*\)',
                rf"{target_mode.complex_constructor}(\1,\2)",
                masked,
                flags=re.IGNORECASE
            )

        # Restore the ``_wp`` literal placeholders.
        if wp_placeholders:
            def _wp_restore(m):
                return wp_placeholders[int(m.group(1))]
            masked = re.sub(r'\x01WP(\d+)\x01', _wp_restore, masked)

        parts[idx] = re.sub(
            r'\x00([A-Za-z]+)\x00', r'.\1.', masked,
        )
    return ''.join(parts)


def _rewrite_int_kind_on_real64x2(
    line: str,
    target_mode: TargetMode,
    real_names: set[str] | None = None,
) -> str:
    """Multifloats only — rewrite ``INT(expr, KIND_INT)`` to
    ``INT(dble(expr), KIND_INT)`` when ``expr`` contains a
    ``real64x2(`` token, or when ``expr`` is a single identifier
    declared as ``TYPE(real64x2)`` in the file (``real_names``).

    Multifloats publishes ``int`` only as ``dd_int(real64x2)`` returning
    default-kind integer — no ``int(real64x2, kind)`` overload exists.
    The migrator's wrap_constructor pass replaces literal real arguments
    with ``real64x2(...)`` calls but leaves the surrounding
    ``INT(..., 8)`` shape intact, leading to "Generic function 'int' is
    not consistent with a specific intrinsic interface" build failures
    in MUMPS files like ``mana_reordertree.F`` / ``mumps_ooc.F``.

    Routing the inner through ``dble`` (multifloats's public generic
    over real64x2 → real(dp), via ``dd_to_dp``) lets the standard
    Fortran INT intrinsic handle the kind selector. Values that flow
    through this pattern (cost counters, matrix size estimates) fit
    comfortably in double precision. ``_build_use_only_clause`` adds
    ``dble`` to the only-list whenever ``INT`` and ``real64x2(``
    co-occur, so the import is in place when this rewrite fires.
    """
    if target_mode.intrinsic_mode != 'wrap_constructor':
        return line
    if not line or line[0] in ('C', 'c', '*', '!'):
        return line
    lower = line.lower()
    # ``int (`` with whitespace before the paren is also a valid call;
    # use the regex below as the authoritative gate, not a substring
    # check that misses the spaced form.
    if not re.search(r'\bint\s*\(', lower):
        return line
    real_names_upper = {n.upper() for n in (real_names or ())}
    if 'real64x2' not in lower and not real_names_upper:
        return line
    pattern = _INT_CALL_RE
    out: list[str] = []
    i = 0
    while i < len(line):
        m = pattern.search(line, i)
        if not m:
            out.append(line[i:])
            break
        out.append(line[i:m.start()])
        paren_open = m.end() - 1
        depth, j = 1, paren_open + 1
        while j < len(line) and depth > 0:
            ch = line[j]
            if ch == '(': depth += 1
            elif ch == ')': depth -= 1
            j += 1
        if depth != 0:
            out.append(line[m.start():])
            break
        inner = line[paren_open + 1:j - 1]
        parts: list[str] = []
        cur, d = '', 0
        for ch in inner:
            if ch == '(': d += 1; cur += ch
            elif ch == ')': d -= 1; cur += ch
            elif ch == ',' and d == 0:
                parts.append(cur); cur = ''
            else:
                cur += ch
        if cur:
            parts.append(cur)
        if (len(parts) == 2
                and re.fullmatch(r'\s*(?:KIND\s*=\s*)?\w+\s*', parts[1])):
            expr = parts[0].strip()
            ident_m = re.fullmatch(r'[A-Za-z_]\w*', expr)
            inner_is_real = (
                'real64x2(' in parts[0].lower()
                or (ident_m is not None
                    and ident_m.group(0).upper() in real_names_upper)
            )
            if inner_is_real:
                head = line[m.start():m.end()]
                kind_arg = parts[1].strip()
                replacement = f'{head}dble({expr}), {kind_arg})'
                out.append(replacement)
            else:
                out.append(line[m.start():j])
        else:
            out.append(line[m.start():j])
        i = j
    return ''.join(out)


def _rewrite_int_of_complex(line: str, complex_names: set[str]) -> str:
    """Wrap the argument of ``INT(...)`` / ``NINT(...)`` with ``MF_REAL``
    when its leading identifier is a known complex variable.

    Multifloats's ``int`` interface only accepts float64x2; calling it
    on complex128x2 fails to dispatch (gfortran does NOT fall back to
    the standard intrinsic for derived-type args, the way it does for
    e.g. integer args). ``MF_REAL`` has overloads for both float64x2
    (identity) and complex128x2 (real-part extraction), so wrapping
    the argument is type-safe.
    """
    if not complex_names or not line:
        return line
    if line[0] in ('C', 'c', '*', '!'):
        return line

    def _process(name_re: re.Pattern, src: str) -> str:
        out = src
        pos = 0
        while True:
            m = name_re.search(out, pos)
            if not m:
                break
            paren_open = m.end() - 1
            depth = 1
            i = paren_open + 1
            while i < len(out) and depth > 0:
                ch = out[i]
                if ch == '(':
                    depth += 1
                elif ch == ')':
                    depth -= 1
                i += 1
            if depth != 0:
                break
            paren_close = i - 1
            inner = out[paren_open + 1:paren_close]
            head = re.match(r'\s*([A-Za-z_]\w*)', inner)
            if head and head.group(1).upper() in complex_names:
                # Use multifloats's overloaded ``real`` generic
                # (cdd_to_dd / dd_to_dd) — it extracts the real
                # component of cmplx64x2 as a real64x2, which the
                # ``int`` / ``nint`` overloads then accept. The earlier
                # ``DD_REAL`` spelling targeted an upstream symbol
                # that no longer exists; ``real`` is in the multifloats
                # generic_names list and gets imported via the USE
                # only-clause whenever a file references it.
                replacement = f'real({inner.strip()})'
                out = out[:paren_open + 1] + replacement + out[paren_close:]
                pos = paren_open + 1 + len(replacement) + 1
            else:
                pos = paren_close + 1
        return out

    line = _process(_INT_CALL_RE, line)
    line = _process(_NINT_CALL_RE, line)
    return line
