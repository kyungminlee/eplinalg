"""Complex-literal wrapping / redundant-constructor unwrapping (Cluster F).

Wraps bare complex literals and strips redundant complex constructors so the
migrated type rewrite stays well-typed. Extracted verbatim from
``fortran_migrator.py``.
"""
import re
import functools

from ..target_mode import TargetMode


def _unwrap_ctor_re(ctor: str) -> re.Pattern:
    """Cache the ``_unwrap_redundant_constructors`` matcher per
    constructor name (``target_mode.real_constructor`` is run-stable)."""
    return re.compile(rf'\b{re.escape(ctor)}\s*\(\s*([A-Za-z_]\w*)\s*\)')


def _unwrap_redundant_constructors(
    line: str, target_mode: TargetMode, real_names: set[str] | None = None,
) -> str:
    """Drop ``float64x2(arg)`` wrappers when ``arg`` is already float64x2.

    The migrator generates ``float64x2(arg)`` for ``DBLE(arg)`` /
    ``REAL(arg)`` calls because the float64x2 interface is the
    universal converter. After the rename pass, ``arg`` may be a bare
    ``MF_*`` constant (which is itself float64x2), or a known
    float64x2 local variable. Wrapping those in ``float64x2(...)``
    fails because multifloats has no identity constructor.
    """
    if target_mode.is_kind_based or not target_mode.real_constructor:
        return line
    ctor = target_mode.real_constructor
    pattern = _unwrap_ctor_re(ctor)

    def _sub(m):
        name = m.group(1)
        upper = name.upper()
        if upper.startswith('MF_'):
            return name
        if real_names and upper in real_names:
            return name
        return m.group(0)

    return pattern.sub(_sub, line)


def _wrap_bare_complex_literals(line: str, target_mode: TargetMode) -> str:
    """Wrap ``(MF_NAME, MF_NAME)`` style complex literals.

    After replace_known_constants has rewritten ``(ZERO, ZERO)`` (etc.)
    to ``(MF_ZERO, MF_ZERO)``, the result is still a Fortran complex
    literal whose components are non-numeric — gfortran rejects this
    in PARAMETER and array-constructor contexts. Wrap with
    ``complex128x2(...)`` so the structure constructor takes over.

    The same fix-up applies to literals whose components have been
    rewritten to ``float64x2('...')`` calls (e.g. when the original
    source used ``(1.0_WP, 0.0_WP)``): the parenthesized pair is no
    longer a recognized Fortran complex literal once the components
    contain function calls, so the outer wrap is necessary.
    """
    if target_mode.is_kind_based or not target_mode.complex_constructor:
        return line
    # Use named-component form (``re=..., im=...``) so the result
    # binds to the complex128x2 *structure constructor* and not to the
    # generic interface — the latter is a function call and is illegal
    # inside PARAMETER initializers (e.g. ``ZONE = complex128x2(re=...,
    # im=...)`` in DGEDMD/UGEDMD). The negative lookbehind on the
    # opening paren prevents re-wrapping an existing
    # ``complex128x2(...)`` call, e.g. one that the migrator already
    # produced from ``DCMPLX(ONE, ZERO)``.
    line = re.sub(
        r'(?<![A-Za-z0-9_])\(\s*([-+]?\s*(?:MF|DD)_[A-Z][A-Z0-9_]*)\s*,\s*'
        r'([-+]?\s*(?:MF|DD)_[A-Z][A-Z0-9_]*)\s*\)',
        rf"{target_mode.complex_constructor}(re=\1, im=\2)",
        line,
    )
    real_ctor = re.escape(target_mode.real_constructor or '')
    line = re.sub(
        rf"(?<![A-Za-z0-9_])\(\s*([-+]?\s*{real_ctor}\([^()]*\))\s*,"
        rf"\s*([-+]?\s*{real_ctor}\([^()]*\))\s*\)",
        rf"{target_mode.complex_constructor}(re=\1, im=\2)",
        line,
    )
    return line
