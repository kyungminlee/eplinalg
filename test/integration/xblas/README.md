# tests/xblas — Differential precision tests for migrated XBLAS

XBLAS (Netlib's Extra-Precise BLAS, 1.0.248) is a C library of
BLAS-shaped routines that take an additional `prec` argument
selecting the internal accumulation precision. The migrated archive
`${LIB_PAIR_PREFIX}xblas` (`libqxxblas.a` for kind16, `libeyxblas.a` for
kind10, `libmwxblas.a` for multifloats) carries Q/X/E/Y/M/W clones of
every uniform-precision `_x` routine with type promotion applied —
e.g. `BLAS_dgemv_x(double *, ...)` → `BLAS_qgemv_x(QREAL *, ...)`
where `QREAL` is `__float128`.

Each test drives one entry point, runs the same operation through
the migrated `_x` routine and through a quad reference, and compares
the results in `REAL(KIND=16)`. Tolerance is
`safety · n_ops · target_eps` per the project-wide convention (see
`tests/README.md`).

## Reference strategy

XBLAS `_x` routines compute the *same mathematical operations* as
ordinary BLAS — `_x` only changes how internally the result is
accumulated. So for routines with a direct BLAS analogue we use the
standard `-freal-8-real-16`-promoted Netlib BLAS already built by
`tests/blas/refblas/` (target `refblas_quad`):

| XBLAS routine | Quad reference                      |
|---|---|
| `BLAS_dgemv_x` / `BLAS_zgemv_x` | `dgemv` / `zgemv`         |
| `BLAS_dgbmv_x` / `BLAS_zgbmv_x` | `dgbmv` / `zgbmv`         |
| `BLAS_dsymv_x`                  | `dsymv`                   |
| `BLAS_zhemv_x`                  | `zhemv`                   |
| `BLAS_dgemm_x` / `BLAS_zgemm_x` | `dgemm` / `zgemm`         |

For routines that XBLAS extends or invents (`ddot_x` adds `alpha`/
`beta` pre-scaling, `axpby`/`sum`/`waxpby`/`gemv2`/`gbmv2`/`symv2`/
`hemv2`/`ge_sum_mv` have no BLAS sibling) we ship hand-coded quad
references in `common/ref_quad_xblas.f90`.

The `prec` argument is always passed `blas_prec_extra` (= 214) — the
highest-accuracy mode the library knows. Under the migrator's
working-precision promotion this collapses to working-precision-on-
working-precision (kind16 head/tail = quad+quad), but the API still
requires the parameter and the dispatch path is the same one
production callers will exercise.

## Calling convention

The migrated archive ships C entry points with the `BLAS_<P>...`
names plus Fortran-callable bridges (compiled with
`-DCONFIG_FC_UNDERSCORE`) named `blas_<P>..._`. Each test wrapper in
`common/target_xblas_body.fypp` declares an explicit `bind(c)`
interface to the bridge, takes inputs in `REAL(KIND=ep)` /
`COMPLEX(KIND=ep)`, converts them to the target's working precision
via `target_conv`, calls the bridge with `prec=blas_prec_extra`, and
converts the result back to `ep`.

XBLAS uses integer enums for `trans`/`uplo`/`conj` (not BLAS's
character flags). The wrappers expose named integer constants
(`blas_no_trans`, `blas_upper`, `blas_no_conj`, …) so test programs
read fluently.

## Coverage

**100% user-facing — 43 family heads.**

| Level | Routines | Tests |
|---|---|---|
| level1 | axpby, dot, sum, waxpby                                                      |  8 |
| level2 | gemv, gemv2, gbmv, gbmv2, ge_sum_mv, symv, symv2, sbmv, spmv, hemv, hemv2, hbmv, hpmv, trmv, tpmv, tbsv, trsv | 32 |
| level3 | gemm, symm, hemm                                                              |  5 |
| | | **43** |

Each routine has one test driver per family head (real and complex).
For routines with no BLAS analogue (XBLAS extensions: axpby, sum,
waxpby, the `*2` head/tail forms, `ge_sum_mv`, complex symmetric
matvec) the reference is hand-coded at REAL(KIND=16) in
`common/ref_quad_xblas.f90`.

| Target       | Status         |
|--------------|---------------|
| kind16       | 43/43 ✓ (passing)  |
| kind10       | 43/43 ✓ (passing)  |
| multifloats  | not yet run        |

## Layout

```
../common/
    prec_kinds.f90       – ep / dp kinds   (shared across suites)
    compare.f90          – rel_err_*       (shared across suites)
    prec_report.F90      – JSON reporter   (shared; strict-exact variant here)
    test_data.f90        – random gens     (shared across suites)
common/
    ref_quad_xblas.f90   – quad refs (BLAS interfaces + hand-coded extras)
    target_xblas_body.fypp – per-target wrapper module template
target_kind10/target_xblas.fypp     – sets prefixes e/y for 80-bit reals
target_kind16/target_xblas.fypp     – sets prefixes q/x for IEEE quad
target_multifloats/target_xblas.fypp – sets prefixes m/w for double-double
level{1,2,3}/test_blas_*.f90         – the test programs
```

## Build / run

```bash
cd $REPO/src
uv run python -m migrator stage /tmp/stg --target kind16 --libraries blas xblas
cmake -S /tmp/stg -B /tmp/stg/build -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/stg/build -j8
ctest --test-dir /tmp/stg/build -R 'xblas_' --output-on-failure
```

Per-test JSON reports land in `/tmp/stg/build/precision_reports/`.
