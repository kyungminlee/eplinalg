# BLACS Type Migration Recipe

This document describes how to migrate the BLACS
(`extern/scalapack-2.2.3/BLACS/SRC/`) from standard precision
(`float`/`double`, `SCOMPLEX`/`DCOMPLEX`) to extended precision
(`long double` / `__float128` for KIND=10/16).

## Why BLACS C Migration Doesn't Need Clang

The Fortran migrator uses Flang's parser because fixed-form Fortran has
complex lexical rules (column constraints, multi-word type keywords like
`DOUBLE PRECISION`, continuation lines). C has none of these problems:

- **Types are single tokens:** `double`, `float`, `SCOMPLEX`, `DCOMPLEX`
- **No column restrictions:** standard free-form C
- **Files are template clones:** `dgebr2d_.c` and `sgebr2d_.c` are
  structurally identical — they differ only in type names and MPI datatype
  constants
- **The header (`Bdef.h`) is the hub:** typedefs, complex struct definitions,
  and Fortran naming `#define`s are centralized

A clang AST parse would be overkill. The migration is best done as:
1. **Modify `Bdef.h`** to add extended-precision type definitions
2. **Clone+sed** the type-specific `.c` files with mechanical substitution
3. **Add naming entries** to the Fortran-calling-convention `#define` blocks

## Source Overview

- **Location:** `extern/scalapack-2.2.3/BLACS/SRC/`
- **Files:** ~172 C source files, 2 header files
- **Language:** C89/C99, no C++
- **Dependencies:** MPI

## Type System

### Current Types (defined in `Bdef.h`)

```c
typedef struct {float r, i;}  SCOMPLEX;   // Fortran COMPLEX
typedef struct {double r, i;} DCOMPLEX;   // Fortran COMPLEX*16
```

Scalar types used directly: `float` (single), `double` (double), `Int` (integer).

### Target Types for KIND=16

```c
// Option A: long double (80-bit on x86, 128-bit on some platforms)
typedef long double QREAL;
typedef struct {long double r, i;} QCOMPLEX;

// Option B: __float128 (GCC/ICC 128-bit, requires libquadmath)
typedef __float128 QREAL;
typedef struct {__float128 r, i;} QCOMPLEX;

// Option C: _Float128 (C23 standard, GCC ≥ 7)
typedef _Float128 QREAL;
typedef struct {_Float128 r, i;} QCOMPLEX;
```

### Target Types for KIND=10

```c
typedef long double EREAL;                   // 80-bit on x86
typedef struct {long double r, i;} ECOMPLEX;
```

### MPI Datatype Mapping

| Source          | KIND=16 Target           | KIND=10 Target      |
|-----------------|--------------------------|---------------------|
| `MPI_FLOAT`     | `MPI_C_LONG_DOUBLE` *    | `MPI_C_LONG_DOUBLE` |
| `MPI_DOUBLE`    | `MPI_C_LONG_DOUBLE` *    | `MPI_C_LONG_DOUBLE` |
| `MPI_COMPLEX`   | custom MPI type **       | custom MPI type **  |
| `MPI_DOUBLE_COMPLEX` | custom MPI type ** | custom MPI type **  |

\* For `__float128`, there is no standard MPI datatype — a custom
`MPI_Datatype` must be created via `MPI_Type_contiguous(sizeof(QREAL), MPI_BYTE, &MPI_QREAL)`.

\** Complex types always need custom MPI types since they are C structs:
`MPI_Type_contiguous(2, MPI_QREAL, &MPI_QCOMPLEX)`.

## File Classification

### Header Files (modify in place)

| File        | Changes Needed                                           |
|-------------|----------------------------------------------------------|
| `Bdef.h`    | Add `QCOMPLEX`/`QREAL` typedefs, new copy macros, new function prototypes, new naming `#define` entries |
| `Bconfig.h` | May need new data type constants (`QUAD`, `COMPLEX32`)   |

### Type-Specific User-Facing Routines (clone d→q, z→x)

Each operation has 5 variants: `i`, `s`, `d`, `c`, `z`.
We add `q` (real extended) and `x` (complex extended).

**General matrix send/recv/broadcast (24 files → +8 new):**

| Operation | d-variant       | → q-variant      | z-variant       | → x-variant      |
|-----------|-----------------|------------------|-----------------|------------------|
| send      | `dgesd2d_.c`    | `qgesd2d_.c`     | `zgesd2d_.c`    | `xgesd2d_.c`     |
| recv      | `dgerv2d_.c`    | `qgerv2d_.c`     | `zgerv2d_.c`    | `xgerv2d_.c`     |
| bcast-send| `dgebs2d_.c`    | `qgebs2d_.c`     | `zgebs2d_.c`    | `xgebs2d_.c`     |
| bcast-recv| `dgebr2d_.c`    | `qgebr2d_.c`     | `zgebr2d_.c`    | `xgebr2d_.c`     |

**Trapezoid matrix send/recv/broadcast (24 files → +8 new):**

| Operation | d-variant       | → q-variant      | z-variant       | → x-variant      |
|-----------|-----------------|------------------|-----------------|------------------|
| send      | `dtrsd2d_.c`    | `qtrsd2d_.c`     | `ztrsd2d_.c`    | `xtrsd2d_.c`     |
| recv      | `dtrrv2d_.c`    | `qtrrv2d_.c`     | `ztrrv2d_.c`    | `xtrrv2d_.c`     |
| bcast-send| `dtrbs2d_.c`    | `qtrbs2d_.c`     | `ztrbs2d_.c`    | `xtrbs2d_.c`     |
| bcast-recv| `dtrbr2d_.c`    | `qtrbr2d_.c`     | `ztrbr2d_.c`    | `xtrbr2d_.c`     |

**Combine operations (12 files → +6 new):**

| Operation | d-variant       | → q-variant      | z-variant       | → x-variant      |
|-----------|-----------------|------------------|-----------------|------------------|
| sum       | `dgsum2d_.c`    | `qgsum2d_.c`     | `zgsum2d_.c`    | `xgsum2d_.c`     |
| max       | `dgamx2d_.c`    | `qgamx2d_.c`     | `zgamx2d_.c`    | `xgamx2d_.c`     |
| min       | `dgamn2d_.c`    | `qgamn2d_.c`     | `zgamn2d_.c`    | `xgamn2d_.c`     |

### Type-Specific Internal Helpers (clone d→q, z→x)

**Vector-vector operations (clone for q and x):**

| d-variant          | → q-variant          | Purpose                    |
|--------------------|----------------------|----------------------------|
| `BI_dvvsum.c`      | `BI_qvvsum.c`        | Vector sum                 |
| `BI_dvvamn.c`      | `BI_qvvamn.c`        | Vector absolute min        |
| `BI_dvvamn2.c`     | `BI_qvvamn2.c`       | Vector abs min (variant)   |
| `BI_dvvamx.c`      | `BI_qvvamx.c`        | Vector absolute max        |
| `BI_dvvamx2.c`     | `BI_qvvamx2.c`       | Vector abs max (variant)   |

| z-variant          | → x-variant          | Purpose                    |
|--------------------|----------------------|----------------------------|
| `BI_zvvsum.c`      | `BI_xvvsum.c`        | Complex vector sum         |
| `BI_zvvamn.c`      | `BI_xvvamn.c`        | Complex vector abs min     |
| `BI_zvvamn2.c`     | `BI_xvvamn2.c`       | Complex vector abs min (v) |
| `BI_zvvamx.c`      | `BI_xvvamx.c`        | Complex vector abs max     |
| `BI_zvvamx2.c`     | `BI_xvvamx2.c`       | Complex vector abs max (v) |

**MPI operation wrappers:**

| d-variant          | → q-variant          |
|--------------------|----------------------|
| `BI_dMPI_amn.c`    | `BI_qMPI_amn.c`      |
| `BI_dMPI_amn2.c`   | `BI_qMPI_amn2.c`     |
| `BI_dMPI_amx.c`    | `BI_qMPI_amx.c`      |
| `BI_dMPI_amx2.c`   | `BI_qMPI_amx2.c`     |
| `BI_dMPI_sum.c`    | `BI_qMPI_sum.c`      |

(Same pattern for z→x variants.)

**Matrix copy operations:**

| d-variant          | → q-variant          |
|--------------------|----------------------|
| `BI_dmvcopy.c`     | `BI_qmvcopy.c`       |
| `BI_dvmcopy.c`     | `BI_qvmcopy.c`       |

(Complex variants use macros that delegate to real copy with 2× width.)

### Precision-Independent Files (no change)

All `BI_*.c` infrastructure files that use `char *` / `void *` buffers
and `MPI_Datatype` parameters:
- `BI_ArgCheck.c`, `BI_Arecv.c`, `BI_Asend.c`, `BI_BeComb.c`
- `BI_BlacsAbort.c`, `BI_BlacsErr.c`, `BI_BlacsWarn.c`
- `BI_BuffIsFree.c`, `BI_ContxtNum.c`, `BI_GetBuff.c`
- `BI_GetMpiGeType.c`, `BI_GetMpiTrType.c`, `BI_GlobalVars.c`
- `BI_HypBR.c`, `BI_HypBS.c`, `BI_IdringBR.c`, `BI_IdringBS.c`
- `BI_MpathBR.c`, `BI_MpathBS.c`, `BI_Pack.c`, `BI_Unpack.c`
- `BI_RingBR.c`, `BI_RingBS.c`, `BI_SringBR.c`, `BI_SringBS.c`
- `BI_Ssend.c`, `BI_TransDist.c`, `BI_TreeBR.c`, `BI_TreeBS.c`
- `BI_TreeComb.c`, `BI_UpdateBuffs.c`
- All `blacs_*.c` (grid management, barrier, etc.)
- Integer variants: `BI_i*.c`, `i*2d_.c`

These work for any type through MPI_Datatype abstraction — no changes needed.

## Transformation Details

### 1. `Bdef.h` Modifications

**Add type definitions:**
```c
// After existing SCOMPLEX/DCOMPLEX:
#ifdef HAVE_FLOAT128
typedef __float128 QREAL;
typedef struct {__float128 r, i;} QCOMPLEX;
#else
typedef long double QREAL;
typedef struct {long double r, i;} QCOMPLEX;
#endif
```

**Add data type constants (in `Bconfig.h` or `Bdef.h`):**
```c
#define QUAD      8     // extended precision real
#define COMPLEX32 9     // extended precision complex
```

**Add copy macros:**
```c
#define BI_qmvcopy(m, n, A, lda, buff)  /* new, uses QREAL */
#define BI_qvmcopy(m, n, A, lda, buff)
#define BI_xmvcopy(m, n, A, lda, buff) \
        BI_qmvcopy(2*(m), (n), (QREAL *)(A), 2*(lda), (QREAL *)(buff))
#define BI_xvmcopy(m, n, A, lda, buff) \
        BI_qvmcopy(2*(m), (n), (QREAL *)(A), 2*(lda), (QREAL *)(buff))
```

**Add function prototypes:**
```c
void BI_qmvcopy(Int m, Int n, QREAL *A, Int lda, QREAL *buff);
void BI_qvmcopy(Int m, Int n, QREAL *A, Int lda, QREAL *buff);
```

**Add naming `#define` entries (in each naming convention block):**
```c
// ADD_ convention:
#define qgesd2d_   qgesd2d_
#define qgerv2d_   qgerv2d_
// ... all q* and x* variants
```

### 2. Type-Specific File Cloning

The substitution rules for cloning `d*` → `q*` files:

| Pattern in d-file         | Replacement in q-file       |
|---------------------------|-----------------------------|
| `double *A`               | `QREAL *A`                  |
| `double *`                | `QREAL *`                   |
| `double diff`             | `QREAL diff`                |
| `(double*)`               | `(QREAL*)`                  |
| `sizeof(double)`          | `sizeof(QREAL)`             |
| `MPI_DOUBLE`              | `MPI_QREAL` (custom type)   |
| `dgebr2d_`                | `qgebr2d_`                  |
| `Cdgebr2d`                | `Cqgebr2d`                  |
| `BI_dvvsum`               | `BI_qvvsum`                 |
| `BI_dvvamn`               | `BI_qvvamn`                 |
| `BI_dMPI_sum`             | `BI_qMPI_sum`               |
| `BI_dmvcopy`              | `BI_qmvcopy`                |
| `BI_dvmcopy`              | `BI_qvmcopy`                |

For cloning `z*` → `x*` (complex extended):

| Pattern in z-file         | Replacement in x-file       |
|---------------------------|-----------------------------|
| `DCOMPLEX`                | `QCOMPLEX`                  |
| `(double*)`               | `(QREAL*)`                  |
| `sizeof(double)`          | `sizeof(QREAL)`             |
| `MPI_DOUBLE_COMPLEX`      | `MPI_QCOMPLEX` (custom)     |
| `zgebr2d_`                | `xgebr2d_`                  |
| `Czgebr2d`                | `Cxgebr2d`                  |
| `BI_zvvsum`               | `BI_xvvsum`                 |
| `BI_zmvcopy`              | `BI_xmvcopy`                |
| `BI_dvvsum` (in z files)  | `BI_qvvsum`                 |
| `BI_dmvcopy` (in z files) | `BI_qmvcopy`                |

### 3. MPI Custom Datatype Registration

A new initialization function is needed to create MPI datatypes:

```c
static MPI_Datatype MPI_QREAL;
static MPI_Datatype MPI_QCOMPLEX;

void BI_init_quad_mpi_types(void) {
    MPI_Type_contiguous(sizeof(QREAL), MPI_BYTE, &MPI_QREAL);
    MPI_Type_commit(&MPI_QREAL);
    MPI_Type_contiguous(2, MPI_QREAL, &MPI_QCOMPLEX);
    MPI_Type_commit(&MPI_QCOMPLEX);
}
```

This must be called during BLACS initialization (`blacs_gridinit_`).

## Why Not Use Clang Parser

Given the above, the BLACS migration is essentially a **templated clone
operation** — not a semantic code transformation. The reasons a parser
isn't needed:

1. **No ambiguity:** `double` always means `double` in C. Unlike Fortran's
   `REAL` which could be `REAL*4`, `REAL*8`, or unspecified.

2. **No format constraints:** C doesn't have column rules. Lines are
   as long as needed.

3. **Mechanical substitution is safe:** Each type-specific file is a
   near-identical copy of its sibling. The substitutions are a fixed set
   of token replacements that don't require understanding control flow.

4. **Macros handle complexity:** The `BI_xmvcopy` → `BI_qmvcopy` delegation
   for complex types is already handled by macros in `Bdef.h`.

5. **No literals to convert:** C doesn't have `1.0D+0` notation. Floating
   constants like `0.0` are generic (they'll be promoted by the compiler
   when assigned to `QREAL`). For explicit precision, we'd use
   `0.0Q` (GCC extension) or cast `(QREAL)0.0`.

A clang parser would be useful if we needed to:
- Understand which variables flow into which function calls
- Track type propagation through complex expressions
- Handle template metaprogramming

None of that applies here. The BLACS C code is intentionally simple.

## Execution Steps

### Step 1: Modify headers

Add extended-precision types, macros, and naming defines to `Bdef.h`.

### Step 2: Clone type-specific files

```bash
# Real: d → q
for f in d*.c BI_d*.c; do
    target=$(echo "$f" | sed 's/^d/q/' | sed 's/BI_d/BI_q/')
    sed -e 's/double/QREAL/g' \
        -e 's/MPI_DOUBLE/MPI_QREAL/g' \
        -e 's/\bd\([a-z]\)/q\1/g' \
        -e 's/Cd\([a-z]\)/Cq\1/g' \
        -e 's/BI_d/BI_q/g' \
        "$f" > "$target"
done

# Complex: z → x
for f in z*.c BI_z*.c; do
    target=$(echo "$f" | sed 's/^z/x/' | sed 's/BI_z/BI_x/')
    sed -e 's/DCOMPLEX/QCOMPLEX/g' \
        -e 's/double/QREAL/g' \
        -e 's/MPI_DOUBLE_COMPLEX/MPI_QCOMPLEX/g' \
        -e 's/MPI_DOUBLE/MPI_QREAL/g' \
        -e 's/\bz\([a-z]\)/x\1/g' \
        -e 's/Cz\([a-z]\)/Cx\1/g' \
        -e 's/BI_z/BI_x/g' \
        -e 's/BI_d/BI_q/g' \
        "$f" > "$target"
done
```

### Step 3: Add MPI type initialization

Create `BI_InitQuadTypes.c` and call from `blacs_gridinit_`.

### Step 4: Update CMakeLists.txt

Add the new `q*` and `x*` source files to the build.

### Step 5: Compile and test

```bash
gcc -c -I. BI_qvvsum.c BI_xvvsum.c qgebr2d_.c xgebr2d_.c ...
```
