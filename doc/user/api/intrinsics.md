# Fortran Generic Intrinsics Involving REAL(KIND=16) or COMPLEX(KIND=16)

All standard Fortran intrinsic procedures (F77 through F2023) that accept or
return `REAL(16)` or `COMPLEX(16)` values, listed by their **generic
(type-agnostic) names only**.  Type-specific names such as `QABS`, `CQABS`,
`DABS`, `CDABS`, etc. are excluded.

This document is a Fortran-language **reference** — it answers "what
generic intrinsics work on quad types out of the box". For the
migrator's actual rewrite table — which type-specific names
(`DABS` → `ABS`, `ZSQRT` → `SQRT`, `DBLE` → `REAL(..., KIND=16)`,
etc.) get translated and whether each call site needs an explicit
`KIND=…` argument — see `codegen/migrator/intrinsics.py`. The two
complement each other: this file describes the surface; the engine
file is the migration rule table.

The same generic-intrinsic surface works for the migrator's other
two targets — `kind10` (`REAL(KIND=10)` / `COMPLEX(KIND=10)`) and
`multifloats` (`TYPE(real64x2)` / `TYPE(cmplx64x2)`) — with these
caveats:

- **kind10**: the table below applies verbatim, replacing
  `KIND=16`/`R16`/`C16` with `KIND=10`/`R10`/`C10`. gfortran (and
  ifx for kind16, but not kind10 — see release-notes for compiler
  caveats) provide the full standard surface natively.
- **multifloats**: the `multifloats` Fortran module provides
  user-defined generic-name overloads for `ABS`, `SQRT`, `EXP`,
  `LOG`, the trigonometric and inverse trig functions, plus
  comparison and arithmetic operators (see
  `codegen/targets/multifloats.yaml`'s `generic_names` and
  `operator_generics` lists). The migrator's intrinsic-rewriter
  emits these in `wrap_constructor` form
  (`ABS(x)` → `real64x2(ABS(REAL(x)))`); the wrapper module makes
  the surface transparent to caller code.

Notation:

- `R16` = `REAL(KIND=16)`
- `C16` = `COMPLEX(KIND=16)`
- `INT` = default `INTEGER`
- `LOG` = default `LOGICAL`
- Standard column: 77 / 90 / 95 / 03 / 08 / 18 / 23 = first Fortran standard
  where the intrinsic (or the specific overload) appeared.

---

## 1. Elemental Mathematical Functions

| Name | Signature (KIND=16 overloads) | Description | Std |
|------|-------------------------------|-------------|-----|
| `ABS` | `ABS(A: R16) -> R16` | Absolute value | 77 |
| | `ABS(A: C16) -> R16` | Absolute value (modulus) of complex | 77 |
| `AINT` | `AINT(A: R16 [,KIND]) -> R16` | Truncation to whole number | 77 |
| `ANINT` | `ANINT(A: R16 [,KIND]) -> R16` | Nearest whole number | 77 |
| `DIM` | `DIM(X: R16, Y: R16) -> R16` | Positive difference max(X-Y, 0) | 77 |
| `MAX` | `MAX(A1: R16, A2: R16 [,...]) -> R16` | Maximum value | 77 |
| `MIN` | `MIN(A1: R16, A2: R16 [,...]) -> R16` | Minimum value | 77 |
| `MOD` | `MOD(A: R16, P: R16) -> R16` | Remainder: A - INT(A/P)*P | 77 |
| `MODULO` | `MODULO(A: R16, P: R16) -> R16` | Modulo (floor division remainder) | 90 |
| `SIGN` | `SIGN(A: R16, B: R16) -> R16` | Transfer of sign: \|A\| * sign(B) | 77 |
| `SQRT` | `SQRT(X: R16) -> R16` | Square root | 77 |
| | `SQRT(X: C16) -> C16` | Complex square root | 77 |
| `EXP` | `EXP(X: R16) -> R16` | Exponential e^X | 77 |
| | `EXP(X: C16) -> C16` | Complex exponential | 77 |
| `LOG` | `LOG(X: R16) -> R16` | Natural logarithm | 77 |
| | `LOG(X: C16) -> C16` | Complex natural logarithm | 77 |
| `LOG10` | `LOG10(X: R16) -> R16` | Common (base-10) logarithm | 77 |
| `HYPOT` | `HYPOT(X: R16, Y: R16) -> R16` | Euclidean distance sqrt(X^2+Y^2) | 08 |

## 2. Trigonometric Functions

| Name | Signature (KIND=16 overloads) | Description | Std |
|------|-------------------------------|-------------|-----|
| `SIN` | `SIN(X: R16) -> R16` | Sine (radians) | 77 |
| | `SIN(X: C16) -> C16` | Complex sine | 77 |
| `COS` | `COS(X: R16) -> R16` | Cosine (radians) | 77 |
| | `COS(X: C16) -> C16` | Complex cosine | 77 |
| `TAN` | `TAN(X: R16) -> R16` | Tangent (radians) | 77 |
| | `TAN(X: C16) -> C16` | Complex tangent | 08 |
| `ASIN` | `ASIN(X: R16) -> R16` | Arc sine | 77 |
| | `ASIN(X: C16) -> C16` | Complex arc sine | 08 |
| `ACOS` | `ACOS(X: R16) -> R16` | Arc cosine | 77 |
| | `ACOS(X: C16) -> C16` | Complex arc cosine | 08 |
| `ATAN` | `ATAN(X: R16) -> R16` | Arc tangent | 77 |
| | `ATAN(Y: R16, X: R16) -> R16` | Arc tangent of Y/X (two-argument form) | 08 |
| | `ATAN(X: C16) -> C16` | Complex arc tangent | 08 |
| `ATAN2` | `ATAN2(Y: R16, X: R16) -> R16` | Arc tangent of Y/X | 77 |

### 2b. Degree-Based Trigonometric Functions (Fortran 2023)

| Name | Signature (KIND=16 overloads) | Description | Std |
|------|-------------------------------|-------------|-----|
| `SIND` | `SIND(X: R16) -> R16` | Sine (degrees) | 23 |
| `COSD` | `COSD(X: R16) -> R16` | Cosine (degrees) | 23 |
| `TAND` | `TAND(X: R16) -> R16` | Tangent (degrees) | 23 |
| `ASIND` | `ASIND(X: R16) -> R16` | Arc sine, result in degrees | 23 |
| `ACOSD` | `ACOSD(X: R16) -> R16` | Arc cosine, result in degrees | 23 |
| `ATAND` | `ATAND(X: R16) -> R16` | Arc tangent, result in degrees | 23 |
| `ATAN2D` | `ATAN2D(Y: R16, X: R16) -> R16` | Arc tangent of Y/X, result in degrees | 23 |

### 2c. Half-Revolution (Pi-Based) Trigonometric Functions (Fortran 2023)

| Name | Signature (KIND=16 overloads) | Description | Std |
|------|-------------------------------|-------------|-----|
| `SINPI` | `SINPI(X: R16) -> R16` | sin(X * pi) | 23 |
| `COSPI` | `COSPI(X: R16) -> R16` | cos(X * pi) | 23 |
| `TANPI` | `TANPI(X: R16) -> R16` | tan(X * pi) | 23 |
| `ASINPI` | `ASINPI(X: R16) -> R16` | asin(X) / pi | 23 |
| `ACOSPI` | `ACOSPI(X: R16) -> R16` | acos(X) / pi | 23 |
| `ATANPI` | `ATANPI(X: R16) -> R16` | atan(X) / pi | 23 |
| `ATAN2PI` | `ATAN2PI(Y: R16, X: R16) -> R16` | atan2(Y, X) / pi | 23 |

## 3. Hyperbolic Functions

| Name | Signature (KIND=16 overloads) | Description | Std |
|------|-------------------------------|-------------|-----|
| `SINH` | `SINH(X: R16) -> R16` | Hyperbolic sine | 77 |
| | `SINH(X: C16) -> C16` | Complex hyperbolic sine | 08 |
| `COSH` | `COSH(X: R16) -> R16` | Hyperbolic cosine | 77 |
| | `COSH(X: C16) -> C16` | Complex hyperbolic cosine | 08 |
| `TANH` | `TANH(X: R16) -> R16` | Hyperbolic tangent | 77 |
| | `TANH(X: C16) -> C16` | Complex hyperbolic tangent | 08 |
| `ASINH` | `ASINH(X: R16) -> R16` | Inverse hyperbolic sine | 08 |
| | `ASINH(X: C16) -> C16` | Complex inverse hyperbolic sine | 08 |
| `ACOSH` | `ACOSH(X: R16) -> R16` | Inverse hyperbolic cosine | 08 |
| | `ACOSH(X: C16) -> C16` | Complex inverse hyperbolic cosine | 08 |
| `ATANH` | `ATANH(X: R16) -> R16` | Inverse hyperbolic tangent | 08 |
| | `ATANH(X: C16) -> C16` | Complex inverse hyperbolic tangent | 08 |

## 4. Special Mathematical Functions

| Name | Signature (KIND=16 overloads) | Description | Std |
|------|-------------------------------|-------------|-----|
| `ERF` | `ERF(X: R16) -> R16` | Error function | 08 |
| `ERFC` | `ERFC(X: R16) -> R16` | Complementary error function 1-ERF(X) | 08 |
| `ERFC_SCALED` | `ERFC_SCALED(X: R16) -> R16` | Scaled complementary error function e^(X^2)*ERFC(X) | 08 |
| `GAMMA` | `GAMMA(X: R16) -> R16` | Gamma function | 08 |
| `LOG_GAMMA` | `LOG_GAMMA(X: R16) -> R16` | Logarithm of absolute value of gamma function | 08 |
| `BESSEL_J0` | `BESSEL_J0(X: R16) -> R16` | Bessel function of first kind, order 0 | 08 |
| `BESSEL_J1` | `BESSEL_J1(X: R16) -> R16` | Bessel function of first kind, order 1 | 08 |
| `BESSEL_JN` | `BESSEL_JN(N: INT, X: R16) -> R16` | Bessel function of first kind, order N (elemental) | 08 |
| | `BESSEL_JN(N1: INT, N2: INT, X: R16) -> R16(:)` | Orders N1 through N2 (transformational) | 08 |
| `BESSEL_Y0` | `BESSEL_Y0(X: R16) -> R16` | Bessel function of second kind, order 0 | 08 |
| `BESSEL_Y1` | `BESSEL_Y1(X: R16) -> R16` | Bessel function of second kind, order 1 | 08 |
| `BESSEL_YN` | `BESSEL_YN(N: INT, X: R16) -> R16` | Bessel function of second kind, order N (elemental) | 08 |
| | `BESSEL_YN(N1: INT, N2: INT, X: R16) -> R16(:)` | Orders N1 through N2 (transformational) | 08 |

## 5. Type Conversion Functions

| Name | Signature (KIND=16 overloads) | Description | Std |
|------|-------------------------------|-------------|-----|
| `REAL` | `REAL(A: INT, KIND=16) -> R16` | Convert integer to R16 | 77 |
| | `REAL(A: R16) -> R16` | Real part (identity for real) | 77 |
| | `REAL(A: C16) -> R16` | Real part of complex | 77 |
| | `REAL(A: any, KIND=16) -> R16` | Convert any numeric to R16 | 90 |
| `CMPLX` | `CMPLX(X: any [,Y: any], KIND=16) -> C16` | Construct COMPLEX(16) from numeric args | 77 |
| `INT` | `INT(A: R16 [,KIND]) -> INT` | Convert to integer (truncation) | 77 |
| | `INT(A: C16 [,KIND]) -> INT` | Convert complex to integer | 77 |
| `NINT` | `NINT(A: R16 [,KIND]) -> INT` | Nearest integer | 77 |
| `CEILING` | `CEILING(A: R16 [,KIND]) -> INT` | Smallest integer >= A | 90 |
| `FLOOR` | `FLOOR(A: R16 [,KIND]) -> INT` | Largest integer <= A | 90 |
| `TRANSFER` | `TRANSFER(SOURCE: any, MOLD: R16 [,SIZE]) -> R16` | Bitwise reinterpretation to R16 | 90 |
| | `TRANSFER(SOURCE: any, MOLD: C16 [,SIZE]) -> C16` | Bitwise reinterpretation to C16 | 90 |
| | `TRANSFER(SOURCE: R16, MOLD: any [,SIZE]) -> type(MOLD)` | Bitwise reinterpretation from R16 | 90 |
| | `TRANSFER(SOURCE: C16, MOLD: any [,SIZE]) -> type(MOLD)` | Bitwise reinterpretation from C16 | 90 |
| `OUT_OF_RANGE` | `OUT_OF_RANGE(X: R16, MOLD: any [,ROUND]) -> LOG` | Test if X is out of range for MOLD type | 18 |
| | `OUT_OF_RANGE(X: any, MOLD: R16 [,ROUND]) -> LOG` | Test if X is out of range for R16 | 18 |

## 6. Complex Number Functions

| Name | Signature (KIND=16 overloads) | Description | Std |
|------|-------------------------------|-------------|-----|
| `CONJG` | `CONJG(Z: C16) -> C16` | Complex conjugate | 77 |
| `AIMAG` | `AIMAG(Z: C16) -> R16` | Imaginary part of complex | 77 |

## 7. Numeric Inquiry Functions

These functions return properties of the numeric model for REAL(16).
The argument value is not used; only its type and kind matter.

| Name | Signature (KIND=16 overloads) | Description | Std |
|------|-------------------------------|-------------|-----|
| `DIGITS` | `DIGITS(X: R16) -> INT` | Number of significant binary digits | 90 |
| `EPSILON` | `EPSILON(X: R16) -> R16` | Smallest X such that 1+X /= 1 | 90 |
| `HUGE` | `HUGE(X: R16) -> R16` | Largest representable number | 90 |
| `MAXEXPONENT` | `MAXEXPONENT(X: R16) -> INT` | Maximum model exponent | 90 |
| `MINEXPONENT` | `MINEXPONENT(X: R16) -> INT` | Minimum model exponent | 90 |
| `PRECISION` | `PRECISION(X: R16) -> INT` | Decimal precision | 90 |
| | `PRECISION(X: C16) -> INT` | Decimal precision of complex kind | 90 |
| `RADIX` | `RADIX(X: R16) -> INT` | Base of the numeric model (always 2) | 90 |
| `RANGE` | `RANGE(X: R16) -> INT` | Decimal exponent range | 90 |
| | `RANGE(X: C16) -> INT` | Decimal exponent range of complex kind | 90 |
| `TINY` | `TINY(X: R16) -> R16` | Smallest positive model number | 90 |
| `KIND` | `KIND(X: R16) -> INT` | Returns kind type parameter (16) | 90 |
| | `KIND(X: C16) -> INT` | Returns kind type parameter (16) | 90 |
| `STORAGE_SIZE` | `STORAGE_SIZE(A: R16 [,KIND]) -> INT` | Storage size in bits | 08 |
| | `STORAGE_SIZE(A: C16 [,KIND]) -> INT` | Storage size in bits | 08 |

## 8. Numeric Model Manipulation Functions

These functions manipulate the floating-point representation.

| Name | Signature (KIND=16 overloads) | Description | Std |
|------|-------------------------------|-------------|-----|
| `EXPONENT` | `EXPONENT(X: R16) -> INT` | Exponent part of floating-point number | 90 |
| `FRACTION` | `FRACTION(X: R16) -> R16` | Fractional part of floating-point model | 90 |
| `NEAREST` | `NEAREST(X: R16, S: R16) -> R16` | Nearest different representable number | 90 |
| `RRSPACING` | `RRSPACING(X: R16) -> R16` | Reciprocal of relative spacing | 90 |
| `SPACING` | `SPACING(X: R16) -> R16` | Absolute spacing of model numbers near X | 90 |
| `SCALE` | `SCALE(X: R16, I: INT) -> R16` | X * RADIX^I | 90 |
| `SET_EXPONENT` | `SET_EXPONENT(X: R16, I: INT) -> R16` | Set exponent of X to I | 90 |

## 9. Array Reduction Functions

These operate on arrays and produce scalar or reduced-rank results.

| Name | Signature (KIND=16 overloads) | Description | Std |
|------|-------------------------------|-------------|-----|
| `SUM` | `SUM(ARRAY: R16 [,DIM] [,MASK]) -> R16` | Sum of array elements | 90 |
| | `SUM(ARRAY: C16 [,DIM] [,MASK]) -> C16` | Sum of complex array elements | 90 |
| `PRODUCT` | `PRODUCT(ARRAY: R16 [,DIM] [,MASK]) -> R16` | Product of array elements | 90 |
| | `PRODUCT(ARRAY: C16 [,DIM] [,MASK]) -> C16` | Product of complex array elements | 90 |
| `DOT_PRODUCT` | `DOT_PRODUCT(A: R16, B: R16) -> R16` | Dot product SUM(A*B) | 90 |
| | `DOT_PRODUCT(A: C16, B: C16) -> C16` | Complex dot product SUM(CONJG(A)*B) | 90 |
| `MAXVAL` | `MAXVAL(ARRAY: R16 [,DIM] [,MASK]) -> R16` | Maximum value in array | 90 |
| `MINVAL` | `MINVAL(ARRAY: R16 [,DIM] [,MASK]) -> R16` | Minimum value in array | 90 |
| `NORM2` | `NORM2(X: R16 [,DIM]) -> R16` | L2 norm (Euclidean norm) | 08 |
| `REDUCE` | `REDUCE(ARRAY: R16, OP [,DIM] [,MASK] [,IDENTITY] [,ORDERED]) -> R16` | General reduction | 18 |
| | `REDUCE(ARRAY: C16, OP [,DIM] [,MASK] [,IDENTITY] [,ORDERED]) -> C16` | General reduction (complex) | 18 |

## 10. Array Location Functions

These accept REAL(16) arrays and return INTEGER index arrays.

| Name | Signature (KIND=16 overloads) | Description | Std |
|------|-------------------------------|-------------|-----|
| `MAXLOC` | `MAXLOC(ARRAY: R16 [,DIM] [,MASK] [,KIND] [,BACK]) -> INT(:)` | Location of maximum value | 90 |
| `MINLOC` | `MINLOC(ARRAY: R16 [,DIM] [,MASK] [,KIND] [,BACK]) -> INT(:)` | Location of minimum value | 90 |
| `FINDLOC` | `FINDLOC(ARRAY: R16, VALUE: R16 [,DIM] [,MASK] [,KIND] [,BACK]) -> INT(:)` | Location of a value | 08 |
| | `FINDLOC(ARRAY: C16, VALUE: C16 [,DIM] [,MASK] [,KIND] [,BACK]) -> INT(:)` | Location of a complex value | 08 |

## 11. Array Transformation Functions

These are type-polymorphic but operate on REAL(16) / COMPLEX(16) arrays.

| Name | Signature (KIND=16 overloads) | Description | Std |
|------|-------------------------------|-------------|-----|
| `MATMUL` | `MATMUL(A: R16, B: R16) -> R16` | Matrix multiplication | 90 |
| | `MATMUL(A: C16, B: C16) -> C16` | Complex matrix multiplication | 90 |
| `TRANSPOSE` | `TRANSPOSE(MATRIX: R16) -> R16` | Matrix transpose | 90 |
| | `TRANSPOSE(MATRIX: C16) -> C16` | Complex matrix transpose | 90 |
| `CSHIFT` | `CSHIFT(ARRAY: R16/C16, SHIFT [,DIM]) -> same` | Circular shift | 90 |
| `EOSHIFT` | `EOSHIFT(ARRAY: R16/C16, SHIFT [,BOUNDARY] [,DIM]) -> same` | End-off shift | 90 |
| `MERGE` | `MERGE(TSOURCE: R16/C16, FSOURCE: R16/C16, MASK) -> same` | Elementwise merge by mask | 90 |
| `PACK` | `PACK(ARRAY: R16/C16, MASK [,VECTOR]) -> R16/C16(:)` | Pack into rank-1 array | 90 |
| `UNPACK` | `UNPACK(VECTOR: R16/C16, MASK, FIELD) -> R16/C16` | Unpack rank-1 into array | 90 |
| `RESHAPE` | `RESHAPE(SOURCE: R16/C16, SHAPE [,PAD] [,ORDER]) -> R16/C16` | Reshape array | 90 |
| `SPREAD` | `SPREAD(SOURCE: R16/C16, DIM, NCOPIES) -> R16/C16` | Replicate along a dimension | 90 |

## 12. Coarray Collective Subroutines (Fortran 2018)

| Name | Signature (KIND=16 overloads) | Description | Std |
|------|-------------------------------|-------------|-----|
| `CO_SUM` | `CO_SUM(A: R16 [,RESULT_IMAGE] [,STAT] [,ERRMSG])` | Sum across images | 18 |
| | `CO_SUM(A: C16 [,RESULT_IMAGE] [,STAT] [,ERRMSG])` | Complex sum across images | 18 |
| `CO_MAX` | `CO_MAX(A: R16 [,RESULT_IMAGE] [,STAT] [,ERRMSG])` | Maximum across images | 18 |
| `CO_MIN` | `CO_MIN(A: R16 [,RESULT_IMAGE] [,STAT] [,ERRMSG])` | Minimum across images | 18 |
| `CO_BROADCAST` | `CO_BROADCAST(A: R16/C16, SOURCE_IMAGE [,STAT] [,ERRMSG])` | Broadcast from one image | 18 |
| `CO_REDUCE` | `CO_REDUCE(A: R16/C16, OP [,RESULT_IMAGE] [,STAT] [,ERRMSG])` | General reduction across images | 18 |

## 13. Intrinsic Subroutines

| Name | Signature (KIND=16 overloads) | Description | Std |
|------|-------------------------------|-------------|-----|
| `RANDOM_NUMBER` | `CALL RANDOM_NUMBER(HARVEST: R16)` | Pseudo-random number(s) in [0, 1) | 90 |
| `CPU_TIME` | `CALL CPU_TIME(TIME: R16)` | Processor time in seconds | 95 |

## Appendix A. IEEE_ARITHMETIC Module Procedures (Fortran 2003+)

These are **module procedures** from `USE IEEE_ARITHMETIC`, not intrinsic
procedures.  They are included here for completeness because they are part of
the Fortran standard library and operate on `REAL(16)` when
`IEEE_SUPPORT_DATATYPE(X)` is `.TRUE.` for KIND=16.

### A1. IEEE Elemental Functions

| Name | Signature (KIND=16 overloads) | Description | Std |
|------|-------------------------------|-------------|-----|
| `IEEE_COPY_SIGN` | `IEEE_COPY_SIGN(X: R16, Y: R16) -> R16` | X with sign of Y | 03 |
| `IEEE_FMA` | `IEEE_FMA(A: R16, B: R16, C: R16) -> R16` | Fused multiply-add (A*B)+C | 18 |
| `IEEE_INT` | `IEEE_INT(A: R16, ROUND [,KIND]) -> INT` | IEEE-conformant conversion to integer | 18 |
| `IEEE_LOGB` | `IEEE_LOGB(X: R16) -> R16` | Unbiased exponent as real | 03 |
| `IEEE_MAX_NUM` | `IEEE_MAX_NUM(X: R16, Y: R16) -> R16` | Maximum, NaN-aware | 18 |
| `IEEE_MAX_NUM_MAG` | `IEEE_MAX_NUM_MAG(X: R16, Y: R16) -> R16` | Max by magnitude, NaN-aware | 18 |
| `IEEE_MIN_NUM` | `IEEE_MIN_NUM(X: R16, Y: R16) -> R16` | Minimum, NaN-aware | 18 |
| `IEEE_MIN_NUM_MAG` | `IEEE_MIN_NUM_MAG(X: R16, Y: R16) -> R16` | Min by magnitude, NaN-aware | 18 |
| `IEEE_NEXT_AFTER` | `IEEE_NEXT_AFTER(X: R16, Y: R16) -> R16` | Next representable number from X toward Y | 03 |
| `IEEE_NEXT_DOWN` | `IEEE_NEXT_DOWN(X: R16) -> R16` | Next representable number less than X | 18 |
| `IEEE_NEXT_UP` | `IEEE_NEXT_UP(X: R16) -> R16` | Next representable number greater than X | 18 |
| `IEEE_REAL` | `IEEE_REAL(A: R16 [,KIND]) -> REAL(KIND)` | IEEE-conformant real conversion | 18 |
| `IEEE_REM` | `IEEE_REM(X: R16, Y: R16) -> R16` | IEEE remainder: X - Y*NINT(X/Y) | 03 |
| `IEEE_RINT` | `IEEE_RINT(X: R16) -> R16` | Round to integer per current rounding mode | 03 |
| `IEEE_SCALB` | `IEEE_SCALB(X: R16, I: INT) -> R16` | X * 2^I | 03 |
| `IEEE_VALUE` | `IEEE_VALUE(X: R16, CLASS) -> R16` | Generate IEEE special value (NaN, Inf, etc.) | 03 |

### A2. IEEE Classification / Query Functions

| Name | Signature (KIND=16 overloads) | Description | Std |
|------|-------------------------------|-------------|-----|
| `IEEE_CLASS` | `IEEE_CLASS(X: R16) -> IEEE_CLASS_TYPE` | Classify value (NaN, Inf, normal, etc.) | 03 |
| `IEEE_IS_FINITE` | `IEEE_IS_FINITE(X: R16) -> LOG` | Is value finite? | 03 |
| `IEEE_IS_NAN` | `IEEE_IS_NAN(X: R16) -> LOG` | Is value NaN? | 03 |
| `IEEE_IS_NEGATIVE` | `IEEE_IS_NEGATIVE(X: R16) -> LOG` | Is value negative? | 03 |
| `IEEE_IS_NORMAL` | `IEEE_IS_NORMAL(X: R16) -> LOG` | Is value normal (not subnormal/NaN/Inf)? | 03 |
| `IEEE_SIGNBIT` | `IEEE_SIGNBIT(X: R16) -> LOG` | Is sign bit set? | 18 |
| `IEEE_UNORDERED` | `IEEE_UNORDERED(X: R16, Y: R16) -> LOG` | Is either value NaN? | 03 |

### A3. IEEE Inquiry Functions

| Name | Signature (KIND=16 overloads) | Description | Std |
|------|-------------------------------|-------------|-----|
| `IEEE_SUPPORT_DATATYPE` | `IEEE_SUPPORT_DATATYPE(X: R16) -> LOG` | Is IEEE arithmetic supported for R16? | 03 |
| `IEEE_SUPPORT_DENORMAL` | `IEEE_SUPPORT_DENORMAL(X: R16) -> LOG` | Are denormals supported? | 03 |
| `IEEE_SUPPORT_DIVIDE` | `IEEE_SUPPORT_DIVIDE(X: R16) -> LOG` | Is IEEE divide supported? | 03 |
| `IEEE_SUPPORT_INF` | `IEEE_SUPPORT_INF(X: R16) -> LOG` | Are infinities supported? | 03 |
| `IEEE_SUPPORT_NAN` | `IEEE_SUPPORT_NAN(X: R16) -> LOG` | Are NaNs supported? | 03 |
| `IEEE_SUPPORT_SQRT` | `IEEE_SUPPORT_SQRT(X: R16) -> LOG` | Is IEEE sqrt supported? | 03 |
| `IEEE_SUPPORT_STANDARD` | `IEEE_SUPPORT_STANDARD(X: R16) -> LOG` | Full IEEE standard supported? | 03 |
| `IEEE_SUPPORT_SUBNORMAL` | `IEEE_SUPPORT_SUBNORMAL(X: R16) -> LOG` | Are subnormals supported? (alias) | 18 |

---

## Summary Count

| Category | Count |
|----------|-------|
| Elemental mathematical | 17 generic names |
| Trigonometric (radians) | 7 generic names |
| Trigonometric (degrees, F2023) | 7 generic names |
| Trigonometric (pi-based, F2023) | 7 generic names |
| Hyperbolic | 6 generic names |
| Special mathematical | 10 generic names |
| Type conversion | 8 generic names |
| Complex number | 2 generic names |
| Numeric inquiry | 12 generic names |
| Numeric model manipulation | 7 generic names |
| Array reduction | 6 generic names |
| Array location | 3 generic names |
| Array transformation | 8 generic names |
| Coarray collective | 5 generic names |
| Intrinsic subroutines | 2 generic names |
| **Total (intrinsic procedures)** | **107 generic names** |
| | |
| IEEE elemental functions (module) | 16 generic names |
| IEEE classification/query (module) | 7 generic names |
| IEEE inquiry (module) | 8 generic names |
| **Total (incl. IEEE module)** | **138 generic names** |

---

## Notes

1. **KIND=16 availability**: Whether `REAL(KIND=16)` is IEEE binary128 quad
   precision or some other 128-bit format is processor-dependent.  GFortran on
   x86-64 uses IEEE binary128; some platforms may use double-double or
   extended-80 padded to 128 bits.  Use `SELECTED_REAL_KIND(33)` to portably
   request ~33 decimal digits.

2. **DPROD** is excluded because it is a type-specific intrinsic (the "D"
   denotes double precision).  Its signature is fixed:
   `DPROD(X: REAL(default), Y: REAL(default)) -> REAL(dp)`.

3. **F2008 complex extensions**: Fortran 2008 extended `TAN`, `ASIN`, `ACOS`,
   `ATAN`, `SINH`, `COSH`, `TANH`, `ASINH`, `ACOSH`, `ATANH` to accept
   `COMPLEX` arguments.  Prior standards only supported `REAL`.

4. **Fortran 2023 trigonometric additions**: Degree-based functions (`SIND`,
   `COSD`, `TAND`, `ASIND`, `ACOSD`, `ATAND`, `ATAN2D`) and pi-based
   functions (`SINPI`, `COSPI`, `TANPI`, `ASINPI`, `ACOSPI`, `ATANPI`,
   `ATAN2PI`) were added in Fortran 2023.  They accept `REAL(16)` but are
   not yet widely implemented.  Some compilers already support degree-based
   functions as extensions.

5. **`SELECTED_REAL_KIND`** is not listed because its arguments and return
   value are `INTEGER`, not `REAL(16)` or `COMPLEX(16)`.

6. **IEEE_ARITHMETIC module procedures** (e.g., `IEEE_VALUE`,
   `IEEE_IS_FINITE`, `IEEE_CLASS`) are module procedures from
   `USE IEEE_ARITHMETIC`, not intrinsic procedures.  They are listed in
   Appendix A for completeness.

7. **`RANDOM_NUMBER`** and **`CPU_TIME`** are intrinsic subroutines (not
   functions) whose `INTENT(OUT)` argument can be `REAL(16)`.

---

## Cross-reference: migrator's rewrite table

The migrator maps **type-specific intrinsic names** (which appear in
upstream Netlib code as `DABS`, `DCONJG`, `DSQRT`, `CDABS`, `IDINT`,
etc.) to **generic names** (`ABS`, `CONJG`, `SQRT`, …) that work on
the target precision automatically. The mapping table is data-only,
in `codegen/migrator/intrinsics.py:INTRINSIC_MAP`. As of 2026-04-30 it
has ~140 entries. Each entry is keyed on the upstream type-specific
name and stores `(generic_name, needs_kind_arg)`:

```python
'DABS':   ('ABS',   False),   # DOUBLE  → DOUBLE
'CABS':   ('ABS',   False),   # COMPLEX → REAL
'ZABS':   ('ABS',   False),   # DOUBLE COMPLEX → DOUBLE
'DSQRT':  ('SQRT',  False),
'DCONJG': ('CONJG', False),
'DBLE':   ('REAL',  True),    # DBLE(x) → REAL(x, KIND=…)
'IDINT':  ('INT',   False),
...
```

`needs_kind_arg = True` means the migrator inserts a `KIND=…` arg at
the call site (e.g. `DBLE(x)` → `REAL(x, KIND=16)` for kind16,
`REAL(x, KIND=10)` for kind10, or just `REAL(x)` wrapped via the
multifloats module's overload). `False` means the generic intrinsic
already dispatches correctly on the argument type.

The table is curated against gfortran's intrinsic registry plus a
handful of GNU/IBM aliases (`CDABS`, `IDFIX`, `DREAL`, …). When a new
upstream library introduces an intrinsic name not yet in the table,
the migrator will leave the call unchanged — fix is to add the entry
(one line) and re-run staging. Test with
`test/integration/lapack/TODO.md`-tracked phases or the convergence report.

This file (`intrinsics.md`) is a Fortran-language reference; the
engine table (`intrinsics.py`) is the migration ruleset. Both are
maintained by hand and may drift if either is edited without
checking the other — when in doubt, the engine table is authoritative
for what the migrator *does*; this doc is authoritative for what the
language *provides*.
