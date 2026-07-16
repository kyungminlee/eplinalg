"""Fortran type-specific intrinsic function mappings.

Complete mapping of every type-specific intrinsic name to its generic
equivalent, derived from GCC gfortran's intrinsic registry
(gcc/fortran/intrinsic.cc).  These are part of the Fortran language
standard (or common extensions) and apply universally to any Fortran
numerical library.

Reference: https://gcc.gnu.org/onlinedocs/gfortran/Intrinsic-Procedures.html
"""

# =====================================================================
# Intrinsic function call replacements
# =====================================================================
#
# Key:   old name (type-specific)
# Value: (new_name, needs_kind_arg)
#
#   needs_kind_arg=False → simple rename: DCONJG(x) → CONJG(x)
#   needs_kind_arg=True  → add KIND argument: DBLE(x) → REAL(x, KIND=k)
#
# Organisation: grouped by generic function, then by input type.

INTRINSIC_MAP: dict[str, tuple[str, bool]] = {

    # -----------------------------------------------------------------
    # ABS — absolute value
    # -----------------------------------------------------------------
    'IABS':   ('ABS',   False),   # INTEGER → INTEGER
    'DABS':   ('ABS',   False),   # DOUBLE  → DOUBLE
    'CABS':   ('ABS',   False),   # COMPLEX → REAL
    'ZABS':   ('ABS',   False),   # DOUBLE COMPLEX → DOUBLE
    'CDABS':  ('ABS',   False),   # DOUBLE COMPLEX (GNU alias)

    # -----------------------------------------------------------------
    # AIMAG — imaginary part of complex
    # -----------------------------------------------------------------
    'DIMAG':  ('AIMAG', False),   # DOUBLE COMPLEX → DOUBLE

    # -----------------------------------------------------------------
    # AINT — truncation to whole number (real result)
    # -----------------------------------------------------------------
    'DINT':   ('AINT',  False),   # DOUBLE → DOUBLE

    # -----------------------------------------------------------------
    # ANINT — nearest whole number (real result)
    # -----------------------------------------------------------------
    'DNINT':  ('ANINT', False),   # DOUBLE → DOUBLE

    # -----------------------------------------------------------------
    # ACOS — arc cosine
    # -----------------------------------------------------------------
    'DACOS':  ('ACOS',  False),   # DOUBLE → DOUBLE

    # -----------------------------------------------------------------
    # ACOSH — inverse hyperbolic cosine  (Fortran 2008)
    # -----------------------------------------------------------------
    'DACOSH': ('ACOSH', False),   # DOUBLE → DOUBLE

    # -----------------------------------------------------------------
    # ASIN — arc sine
    # -----------------------------------------------------------------
    'DASIN':  ('ASIN',  False),   # DOUBLE → DOUBLE

    # -----------------------------------------------------------------
    # ASINH — inverse hyperbolic sine  (Fortran 2008)
    # -----------------------------------------------------------------
    'DASINH': ('ASINH', False),   # DOUBLE → DOUBLE

    # -----------------------------------------------------------------
    # ATAN — arc tangent
    # -----------------------------------------------------------------
    'DATAN':  ('ATAN',  False),   # DOUBLE → DOUBLE

    # -----------------------------------------------------------------
    # ATAN2 — arc tangent of y/x (two arguments)
    # -----------------------------------------------------------------
    'DATAN2': ('ATAN2', False),   # DOUBLE → DOUBLE

    # -----------------------------------------------------------------
    # ATANH — inverse hyperbolic tangent  (Fortran 2008)
    # -----------------------------------------------------------------
    'DATANH': ('ATANH', False),   # DOUBLE → DOUBLE

    # -----------------------------------------------------------------
    # CONJG — complex conjugate
    # -----------------------------------------------------------------
    'DCONJG': ('CONJG', False),   # DOUBLE COMPLEX → DOUBLE COMPLEX

    # -----------------------------------------------------------------
    # COS — cosine
    # -----------------------------------------------------------------
    'DCOS':   ('COS',   False),   # DOUBLE → DOUBLE
    'CCOS':   ('COS',   False),   # COMPLEX → COMPLEX
    'ZCOS':   ('COS',   False),   # DOUBLE COMPLEX → DOUBLE COMPLEX
    'CDCOS':  ('COS',   False),   # DOUBLE COMPLEX (GNU alias)

    # -----------------------------------------------------------------
    # COSH — hyperbolic cosine
    # -----------------------------------------------------------------
    'DCOSH':  ('COSH',  False),   # DOUBLE → DOUBLE

    # -----------------------------------------------------------------
    # DIM — positive difference: max(x-y, 0)
    # -----------------------------------------------------------------
    'IDIM':   ('DIM',   False),   # INTEGER → INTEGER
    'DDIM':   ('DIM',   False),   # DOUBLE → DOUBLE

    # -----------------------------------------------------------------
    # ERF — error function  (Fortran 2008)
    # -----------------------------------------------------------------
    'DERF':   ('ERF',   False),   # DOUBLE → DOUBLE

    # -----------------------------------------------------------------
    # ERFC — complementary error function  (Fortran 2008)
    # -----------------------------------------------------------------
    'DERFC':  ('ERFC',  False),   # DOUBLE → DOUBLE

    # -----------------------------------------------------------------
    # EXP — exponential
    # -----------------------------------------------------------------
    'DEXP':   ('EXP',   False),   # DOUBLE → DOUBLE
    'CEXP':   ('EXP',   False),   # COMPLEX → COMPLEX
    'ZEXP':   ('EXP',   False),   # DOUBLE COMPLEX → DOUBLE COMPLEX
    'CDEXP':  ('EXP',   False),   # DOUBLE COMPLEX (GNU alias)

    # -----------------------------------------------------------------
    # GAMMA — gamma function  (Fortran 2008)
    # -----------------------------------------------------------------
    'DGAMMA': ('GAMMA', False),   # DOUBLE → DOUBLE

    # -----------------------------------------------------------------
    # INT — conversion to integer (truncation)
    # -----------------------------------------------------------------
    'IFIX':   ('INT',   False),   # REAL → INTEGER
    'IDINT':  ('INT',   False),   # DOUBLE → INTEGER

    # -----------------------------------------------------------------
    # LOG — natural logarithm
    # -----------------------------------------------------------------
    'ALOG':   ('LOG',   False),   # REAL → REAL
    'DLOG':   ('LOG',   False),   # DOUBLE → DOUBLE
    'CLOG':   ('LOG',   False),   # COMPLEX → COMPLEX
    'ZLOG':   ('LOG',   False),   # DOUBLE COMPLEX → DOUBLE COMPLEX
    'CDLOG':  ('LOG',   False),   # DOUBLE COMPLEX (GNU alias)

    # -----------------------------------------------------------------
    # LOG10 — base-10 logarithm
    # -----------------------------------------------------------------
    'ALOG10': ('LOG10', False),   # REAL → REAL
    'DLOG10': ('LOG10', False),   # DOUBLE → DOUBLE

    # -----------------------------------------------------------------
    # LOG_GAMMA / LGAMMA — log of absolute value of gamma function
    # -----------------------------------------------------------------
    'ALGAMA': ('LOG_GAMMA', False),  # REAL (GNU)
    'DLGAMA': ('LOG_GAMMA', False),  # DOUBLE (GNU)

    # -----------------------------------------------------------------
    # MAX — maximum value
    # -----------------------------------------------------------------
    'MAX0':   ('MAX',   False),   # INTEGER, INTEGER → INTEGER
    'MAX1':   ('MAX',   False),   # REAL, REAL → INTEGER
    'AMAX0':  ('MAX',   False),   # INTEGER, INTEGER → REAL
    'AMAX1':  ('MAX',   False),   # REAL, REAL → REAL
    'DMAX1':  ('MAX',   False),   # DOUBLE, DOUBLE → DOUBLE

    # -----------------------------------------------------------------
    # MIN — minimum value
    # -----------------------------------------------------------------
    'MIN0':   ('MIN',   False),   # INTEGER, INTEGER → INTEGER
    'MIN1':   ('MIN',   False),   # REAL, REAL → INTEGER
    'AMIN0':  ('MIN',   False),   # INTEGER, INTEGER → REAL
    'AMIN1':  ('MIN',   False),   # REAL, REAL → REAL
    'DMIN1':  ('MIN',   False),   # DOUBLE, DOUBLE → DOUBLE

    # -----------------------------------------------------------------
    # MOD — remainder
    # -----------------------------------------------------------------
    'AMOD':   ('MOD',   False),   # REAL → REAL
    'DMOD':   ('MOD',   False),   # DOUBLE → DOUBLE

    # -----------------------------------------------------------------
    # NINT — nearest integer
    # -----------------------------------------------------------------
    'IDNINT': ('NINT',  False),   # DOUBLE → INTEGER

    # -----------------------------------------------------------------
    # SIGN — transfer of sign: |a| * sign(b)
    # -----------------------------------------------------------------
    'ISIGN':  ('SIGN',  False),   # INTEGER → INTEGER
    'DSIGN':  ('SIGN',  False),   # DOUBLE → DOUBLE

    # -----------------------------------------------------------------
    # SIN — sine
    # -----------------------------------------------------------------
    'DSIN':   ('SIN',   False),   # DOUBLE → DOUBLE
    'CSIN':   ('SIN',   False),   # COMPLEX → COMPLEX
    'ZSIN':   ('SIN',   False),   # DOUBLE COMPLEX → DOUBLE COMPLEX
    'CDSIN':  ('SIN',   False),   # DOUBLE COMPLEX (GNU alias)

    # -----------------------------------------------------------------
    # SINH — hyperbolic sine
    # -----------------------------------------------------------------
    'DSINH':  ('SINH',  False),   # DOUBLE → DOUBLE

    # -----------------------------------------------------------------
    # SQRT — square root
    # -----------------------------------------------------------------
    'DSQRT':  ('SQRT',  False),   # DOUBLE → DOUBLE
    'CSQRT':  ('SQRT',  False),   # COMPLEX → COMPLEX
    'ZSQRT':  ('SQRT',  False),   # DOUBLE COMPLEX → DOUBLE COMPLEX
    'CDSQRT': ('SQRT',  False),   # DOUBLE COMPLEX (GNU alias)

    # -----------------------------------------------------------------
    # TAN — tangent
    # -----------------------------------------------------------------
    'DTAN':   ('TAN',   False),   # DOUBLE → DOUBLE

    # -----------------------------------------------------------------
    # TANH — hyperbolic tangent
    # -----------------------------------------------------------------
    'DTANH':  ('TANH',  False),   # DOUBLE → DOUBLE

    # =================================================================
    # Type conversion intrinsics — need KIND argument to preserve the
    # target precision after migration.
    # =================================================================

    # -----------------------------------------------------------------
    # REAL — conversion to real  (needs KIND to select target precision)
    # -----------------------------------------------------------------
    'DBLE':   ('REAL',  True),    # any → DOUBLE (specific to double)
    'SNGL':   ('REAL',  True),    # DOUBLE → REAL (specific to single)
    'FLOAT':  ('REAL',  True),    # INTEGER → REAL (specific to single)
    'DFLOAT': ('REAL',  True),    # INTEGER → DOUBLE (GNU extension)
    'DREAL':  ('REAL',  True),    # DOUBLE COMPLEX → DOUBLE (GNU ext.)

    # -----------------------------------------------------------------
    # CMPLX — conversion to complex  (needs KIND to select target prec.)
    # -----------------------------------------------------------------
    'DCMPLX': ('CMPLX', True),   # any → DOUBLE COMPLEX (GNU extension)
    # Bare CMPLX needs an explicit KIND argument too. The original
    # precision-specific source used default-kind complex because the
    # file was single-precision; post-migration the result must carry
    # the target precision. The ``needs_kind`` path skips any call
    # that already has a top-level ``KIND=`` argument, so calls left
    # over from the ``DCMPLX``→``CMPLX`` rewrite above are not
    # double-annotated.
    'CMPLX':  ('CMPLX', True),
}


# =====================================================================
# INTRINSIC declaration statement replacements
# =====================================================================
#
# When a name appears in an  INTRINSIC  declaration, the type-specific
# name should be replaced by its generic equivalent.  This map mirrors
# INTRINSIC_MAP but uses plain str→str (no needs_kind flag); the
# ``k != v[0]`` filter drops the bare ``CMPLX → CMPLX`` self-rename,
# which only appears in call-site rewriting (where it gates KIND
# annotation), never in INTRINSIC declarations.

INTRINSIC_DECL_MAP: dict[str, str] = {
    k: v[0] for k, v in INTRINSIC_MAP.items() if k != v[0]
}
