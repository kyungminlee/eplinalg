# Developer Guide

See [`../README.md`](../README.md) for the full documentation index.

## Project Goal

**eplinalg** converts BLAS, LAPACK, BLACS, ScaLAPACK, and similar
numerical Fortran libraries from standard precision (`REAL`/`DOUBLE PRECISION`,
`COMPLEX`/`DOUBLE COMPLEX`) to extended precision (`KIND=10` for 80-bit,
`KIND=16` for 128-bit) or multiword floating-point (`float64x2` via the
multifloats target). This involves:

1. Rewriting type declarations
2. Renaming subroutine/function names (both definitions and call sites)
3. Updating literal constants
4. Converting type-specific intrinsic calls
5. Updating `EXTERNAL` declarations
6. Renaming output file names

The tool must preserve all preprocessor directives, comments, and formatting.

---

## Approach: Hybrid (Parser + Source-Level Rewriting)

### Why Not Pure String Manipulation?

Fixed-form Fortran has complex rules (columns 1-5 labels, column 6 continuation,
columns 7-72 statements). A type declaration can span multiple continuation
lines, appear inside preprocessor-guarded blocks, or use non-standard syntax
like `REAL*8`, `COMPLEX*16`, `DOUBLE PRECISION`. Simple regex is brittle and
error-prone — it risks false positives in comments, string literals, or
identifiers that happen to contain type keywords.

### Why Not Pure Flang Unparse?

Flang's `flang-new -fdebug-unparse` can round-trip source code, but it
normalizes formatting, rewrites continuation lines, and may alter whitespace
in ways that make the output differ from the original. For libraries like
LAPACK that have very precise fixed-form formatting and extensive comment
blocks, this is unacceptable. Preprocessing directives can also be lost or
rearranged.

### The Hybrid Strategy

1. **Parse with a Fortran compiler** (Flang or GFortran, selected via `--parser`)
   to produce a parse tree with source provenance. Without `--parser`, the tool
   falls back to regex-only scanning.
2. **Walk the parse tree** to identify all locations that need transformation:
   type declarations, function/subroutine names, `EXTERNAL` statements, literal
   constants, intrinsic function calls, `CALL` statements, and string literals
   containing routine names (e.g., `CALL XERBLA('DGEMM ', INFO)`).
3. **Apply transformations as source-level regex replacements**, preserving
   everything else (comments, preprocessor directives, whitespace, column
   alignment) exactly as-is.
4. **Handle column-width constraints** for fixed-form: when a replacement string
   is longer than the original, adjust spacing or split into continuation lines
   as needed while staying within columns 7-72.

This gives us the accuracy of a parser-guided approach with the preservation
guarantees of text-level patching.

---

## What Gets Transformed

### 1. Type Declarations

Source precision types and their target mappings:

| Source Type            | KIND=10 Target     | KIND=16 Target     |
|------------------------|--------------------|--------------------|
| `REAL`                 | `REAL(KIND=10)`    | `REAL(KIND=16)`    |
| `REAL*4`               | `REAL(KIND=10)`    | `REAL(KIND=16)`    |
| `REAL*8`               | `REAL(KIND=10)`    | `REAL(KIND=16)`    |
| `DOUBLE PRECISION`     | `REAL(KIND=10)`    | `REAL(KIND=16)`    |
| `REAL(KIND=WP)`        | `REAL(KIND=10)`    | `REAL(KIND=16)`    |
| `COMPLEX`              | `COMPLEX(KIND=10)` | `COMPLEX(KIND=16)` |
| `COMPLEX*8`            | `COMPLEX(KIND=10)` | `COMPLEX(KIND=16)` |
| `COMPLEX*16`           | `COMPLEX(KIND=10)` | `COMPLEX(KIND=16)` |
| `DOUBLE COMPLEX`       | `COMPLEX(KIND=10)` | `COMPLEX(KIND=16)` |
| `COMPLEX(KIND=WP)`     | `COMPLEX(KIND=10)` | `COMPLEX(KIND=16)` |

Note: the source files may use *either* single-precision (`REAL`, `COMPLEX`)
or double-precision (`DOUBLE PRECISION`, `DOUBLE COMPLEX`, `REAL*8`,
`COMPLEX*16`) depending on the routine variant (`s*`/`c*` vs `d*`/`z*`).
All floating-point types get mapped to the target kind.

For the multifloats target, `REAL` → `TYPE(float64x2)` and
`COMPLEX` → `TYPE(complex64x2)`, with constructor wrapping instead of
KIND suffixes.

### 2. Subroutine/Function Name Prefixes

Standard BLAS/LAPACK naming convention uses a single-character prefix:

| Source Prefix | Meaning           | KIND=10 | KIND=16 | Multifloats |
|---------------|-------------------|---------|---------|-------------|
| `S`           | Single real       | `E`     | `Q`     | `DD`        |
| `D`           | Double real       | `E`     | `Q`     | `DD`        |
| `C`           | Single complex    | `Y`     | `X`     | `ZZ`        |
| `Z`           | Double complex    | `Y`     | `X`     | `ZZ`        |

Names must be updated in:
- `SUBROUTINE` / `FUNCTION` / `ENTRY` definition lines
- `CALL` statements
- `EXTERNAL` declarations
- Function references in expressions (e.g., `TEMP = DDOT(N, X, 1, Y, 1)`)
- String literals passed to `XERBLA` (e.g., `'DGEMM '` → `'QGEMM '`)

For ScaLAPACK, the prefix is `P` + precision character (e.g., `PDGESV` →
`PQGESV`). The `P` prefix is preserved; only the precision character changes.

### 3. Literal Constants

| Source Literal  | KIND=10 Replacement    | KIND=16 Replacement    |
|-----------------|------------------------|------------------------|
| `1.0E+0`        | `1.0E+0_10`            | `1.0E+0_16`            |
| `1.0D+0`        | `1.0E+0_10`            | `1.0E+0_16`            |
| `0.0D0`         | `0.0E0_10`             | `0.0E0_16`             |

The `D` exponent letter is specific to `DOUBLE PRECISION` and must be replaced
with `E` when using explicit `KIND`. For the multifloats target, literals are
wrapped in constructors instead: `1.0D+0` → `float64x2(1.0D+0)`.

### 4. Intrinsic Functions

Type-specific intrinsics must be converted to their generic or
extended-precision equivalents:

| Source Intrinsic | Generic Equivalent | Notes |
|------------------|-------------------|-------|
| `DBLE(x)`        | `REAL(x, KIND=k)` | Type conversion |
| `REAL(x)`        | `REAL(x, KIND=k)` | Add KIND argument |
| `DCMPLX(x,y)`   | `CMPLX(x, y, KIND=k)` | |
| `DCONJG(z)`      | `CONJG(z)` | Generic works for all kinds |
| `DIMAG(z)`       | `AIMAG(z)` | Generic works for all kinds |
| `DABS(x)`        | `ABS(x)` | Generic |
| `DSQRT(x)`       | `SQRT(x)` | Generic |
| `DEXP(x)`        | `EXP(x)` | Generic |
| `DLOG(x)`        | `LOG(x)` | Generic |
| `DSIN(x)` / `DCOS(x)` | `SIN(x)` / `COS(x)` | Generic |
| `DSIGN(x,y)`     | `SIGN(x,y)` | Generic |
| `DMAX1(x,y)`     | `MAX(x,y)` | Generic |
| `DMIN1(x,y)`     | `MIN(x,y)` | Generic |
| `DNINT(x)`       | `ANINT(x)` | Generic |
| `IDNINT(x)`      | `NINT(x)` | Generic |

See [intrinsics.md](intrinsics.md) for the complete catalog.

### 5. Machine-Parameter Routines

`DLAMCH` / `SLAMCH` provide machine epsilon, overflow/underflow thresholds, etc.
These routines must be regenerated for the target precision — they cannot simply
be renamed because the constants they return are precision-dependent. The
migrated version (`QLAMCH` / `ELAMCH` / `XLAMCH` / `YLAMCH`) needs to compute
and return the correct machine parameters for the target `KIND`.

### 6. File Names

Output files are renamed to reflect the new prefix:
- `dgemm.f` → `qgemm.f`
- `zgesv.f` → `xgesv.f`
- `pdgesv.f` → `pqgesv.f`

---

## Symbol Database

The tool uses a **symbol scanner** (`src/migrator/symbol_scanner.py`) to discover
routine names from source files or compiled libraries, then a **prefix
classifier** (`src/migrator/prefix_classifier.py`) to classify each by its precision
prefix (S/D/C/Z) and build a rename map. Only symbols with confirmed precision
variants are renamed, avoiding false positives on identifiers that happen to
start with `d`/`s`/`c`/`z`.

For source-only distributions, `scan_source` is the default method. For
libraries with pre-built binaries, the `nm_library` method uses `nm` to extract
symbols from `.a`/`.so`/`.dylib` files.

---

## Target Libraries

### BLAS (Basic Linear Algebra Subprograms)
- **Source format**: Fixed-form (`.f`)
- **Naming**: Single-character prefix (`DGEMM`, `SGEMV`, `ZHER2K`)
- **Location**: `external/lapack-3.12.1/BLAS/SRC/`
- **Notes**: Pure Fortran, no preprocessor usage.

### LAPACK (Linear Algebra PACKage)
- **Source format**: Fixed-form (`.f`)
- **Naming**: Single-character prefix (`DGESV`, `ZHEEV`, `SGEQRF`)
- **Location**: `external/lapack-3.12.1/SRC/`
- **Notes**: Large library (~2000 source files). Pure Fortran, heavy use of
  `EXTERNAL` declarations and cross-routine calls.

### BLACS (Basic Linear Algebra Communication Subprograms)
- **Source format**: C (`.c`) with headers
- **Location**: `external/scalapack-2.2.3/BLACS/SRC/`
- **Notes**: Implemented entirely in C. Migration via template-based cloning.
  The `Bdef.h` header is patched to add extended-precision typedefs and MPI
  datatypes. For multifloats, `complex128x2_t` uses `.r`/`.i` members matching
  the BLACS `DCOMPLEX`/`SCOMPLEX` struct convention, so `Rabs`/`Cabs` macros
  work without overrides.

### ScaLAPACK (Scalable LAPACK)
- **Source format**: Fixed-form Fortran (`.f`) and some C (`.c`)
- **Naming**: `P` + precision prefix (`PDGESV`, `PZHEEV`)
- **Location**: `external/scalapack-2.2.3/SRC/`
- **Notes**: Follows same patterns as LAPACK with distributed-memory descriptors.
  Calls BLACS and PBLAS routines.

### PBLAS (Parallel BLAS)
- **Source format**: Mixed C and Fortran
- **Location**: `external/scalapack-2.2.3/PBLAS/SRC/`
- **Notes**: C wrapper layer around BLAS with MPI. Uses preprocessor macros
  extensively. Contains Fortran kernel subdirectories (PBBLAS, PTZBLAS).

---

## Fortran Format Considerations

### Fixed-Form (`.f`, `.for`)
- **Columns 1-5**: Statement label (numeric)
- **Column 6**: Continuation character (any non-blank, non-zero; typically `+`, `$`, `&`, or `*`)
- **Columns 7-72**: Statement body
- **Column 73+**: Ignored (historically sequence numbers)
- **Comments**: `C`, `c`, `*`, or `!` in column 1

When a replacement makes a line exceed column 72, the tool must:
1. Try to reduce whitespace within the statement
2. If still too long, split into continuation lines
3. Preserve the continuation character style used in the file

### Free-Form (`.f90`, `.f95`, `.F90`)
- No column restrictions (but typically 132 characters max)
- `&` at end of line for continuation; `&` at start of next line optional
- `!` for inline comments
- Case insensitive (like all Fortran)

### Preprocessor Directives
- Lines starting with `#` (after optional whitespace in free-form)
- Must be preserved verbatim — never modified
- Code inside `#ifdef` blocks should still be transformed

---

## Edge Cases and Pitfalls

1. **XERBLA string arguments**: Routine names appear as string literals in
   `CALL XERBLA('DGEMM ', INFO)`. These must be renamed too.

2. **Mixed-precision routines**: Some routines call across precisions (e.g.,
   `DSDOT` calls both single and double). These need careful analysis — the
   routine itself may not have a direct quad equivalent.

3. **Machine constants**: `DLAMCH` returns precision-specific constants.
   The migrated version must recompute these for the target precision.

4. **SAVE'd variables and DATA statements**: Type-specific constants in `DATA`
   statements must be updated.

5. **EQUIVALENCE statements**: These create type aliasing that may break if
   type sizes change. Flag these for manual review.

6. **Fortran `INCLUDE`**: Some files use `INCLUDE 'filename'` which is not a
   preprocessor directive but a Fortran statement. Included files may contain
   type declarations that need migration.

7. **Comment-only name references**: Names in comments (especially the header
   documentation blocks) should ideally be updated for consistency but are not
   semantically critical.

8. **`PARAMETER` constants**: e.g., `PARAMETER (ONE=1.0D+0, ZERO=0.0D+0)` —
   both the literal format and the variable type must be consistent. For
   multifloats, these can be replaced with module-provided constants
   (`DD_ONE`, `DD_ZERO`) via the `known_constants` recipe mechanism.

9. **C source files (BLACS, PBLAS)**: Migration via template-clone with
   mechanical text substitution on C types (`double`→`__float128` or
   `float64x2_t`), MPI datatypes, and function-name prefixes. No clang parser
   is needed. Regex boundaries use `[^a-zA-Z_0-9]` to prevent re-matching
   inside replacement types (e.g., `float` inside `float64x2_t`).

10. **Non-prefixed helper routines**: Some LAPACK internal routines don't
    follow the prefix convention (e.g., `LSAME`, `XERBLA`, `ILAENV`).
    These are type-independent and should *not* be renamed.

---

## Testing Strategy

### 1. Round-Trip Test
Parse → identify → apply null transformation → verify output is byte-identical
to input.

### 2. Known-Good Conversions
Manually create expected output for a small set of files (`dgemm.f` → `qgemm.f`,
`zgesv.f` → `xgesv.f`) and verify the tool produces identical results.

### 3. Dual-Origin Convergence Test

Every BLAS/LAPACK routine exists in both single- and double-precision variants.
When both are migrated to the same target, they should produce identical output.
This is implemented as the `converge` CLI command (see [usage.md](usage.md)).

**Definition — divergence.** Take the two precision variants of the same
routine (the `s`/`c` single-precision source and the `d`/`z` double-precision
source). Run type migration on both into the same target type (e.g. `q` /
`KIND=16`). Compare the two resulting files. Any remaining difference — in
declarations, literals, intrinsics, parameters, or anywhere else — is a
*divergence*. Zero divergence means the migration is source-agnostic: the
target code does not carry a memory of which precision it started from.

This matters because (a) a divergence-free migration can be driven from
either side alone, and (b) divergences usually flag a place where a purely
syntactic rule is papering over a semantic asymmetry between the sources
(e.g. `DOUBLE PRECISION` meaning "working precision" in a `d*` file but
"timing / MPI interface" in an `s*` file — the same token, two intents).

Run `python -m migrator diverge <recipe>` for current per-target divergence counts.

### 3b. Divergence Patterns by Target

The dual-origin convergence test (§3) reveals different divergence counts
depending on the target mode.  Two target families are supported:

| Target | Type Declarations | Literal Mode | Intrinsic Mode |
|--------|-------------------|--------------|----------------|
| `kind16` | `REAL(KIND=16)` | `_16` suffix | `REAL(x, KIND=16)` |
| `multifloats` | `TYPE(float64x2)` | `float64x2(...)` constructor | Constructor wrapping |

**Why multifloats has fewer divergences than KIND-based targets:**

The multifloats target eliminates certain S/D asymmetries that KIND targets
preserve.  Two mechanisms drive this:

1. **Named-constant replacement.**  The `known_constants` map in the target
   YAML replaces local `PARAMETER` constants (`ONE`, `ZERO`, `HALF`, etc.)
   with module-provided values (`DD_ONE`, `DD_ZERO`, `DD_HALF`).  Many
   LAPACK D/Z routines define these constants locally while their S/C
   counterparts do not.  In KIND mode the D-half's `PARAMETER (ONE=1.0E0_16)`
   declaration survives into the output, creating a divergence.  In
   multifloats mode both halves use `DD_ONE` from `USE MULTIFLOATS`, so
   the local declaration is removed and the pair converges.

   *Example:* `dgehd2.f` defines `DOUBLE PRECISION ONE; PARAMETER(ONE=1.0D+0)`.
   `sgehd2.f` has no such constant.
   - KIND=16: D-half output has `REAL(KIND=16) ONE; PARAMETER(ONE=1.0E0_16)` → diverges (+2).
   - Multifloats: `ONE` is replaced by `DD_ONE` in both halves → converges.

2. **Constructor wrapping normalizes type conversions.**  In KIND mode,
   `REAL(x, KIND=16)` preserves the original intrinsic call structure, so
   differences in whether the S-half or D-half uses an explicit `REAL()`
   show up as divergences.  In multifloats mode, both `REAL(x)` and bare
   usage get wrapped uniformly in `float64x2(x)`, collapsing the
   distinction.

### 4. Compilation Test
Compile the migrated library with a Fortran compiler that supports `KIND=16`
(e.g., GFortran with `-freal-8-real-16`) and verify it compiles without errors.

### 5. Symbol Completeness
After migration, extract symbols from the migrated library and verify every
`D*` symbol has a corresponding `Q*` symbol.

### 6. Numerical Test
For a subset of routines, run the LAPACK test suite with the migrated library
and verify results match expected precision.
