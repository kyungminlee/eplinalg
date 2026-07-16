# tests/xblas — CHANGELOG

Resolved items, reverse-chronological. Open work lives in `TODO.md`.

## 2026-04-30 — Other targets

All three targets (kind10 / kind16 / multifloats) build and pass without
per-target adjustment. 43/43 tests on multifloats, 85/85 case-level PASS.
The wrapper module's fypp parameterization on the prefix-character pair
was sufficient; no header_patch / tolerance bumps required for the
dd-on-double target.

## LAPACK xx-family unblocking (P17xx)

`tests/lapack/TODO.md` §P17xx documented 18 LAPACK extra-precise
iterative-refinement drivers (`gesvxx`, `gbsvxx`, `gerfsx`, ...) blocked
on XBLAS being available. Now landed — 18/18 PASS on all three targets.

Resolution: `${LIB_PREFIX}lapack` PUBLIC-links `${LIB_PREFIX}xblas` so
migrated qla/ela/mla `*_extended` objects resolve via
qxblas/exblas/mxblas, and a quad-precision Fortran bridge
(`tests/lapack/reflapack/refxblas_quad_bridge.f90`) supplies the
upstream-named `blas_dgemv_x_` / `blas_zhemv_x_` / etc. that
reflapack_quad's `dla_*_extended.f` calls — forwarding to the
quad-promoted refblas_quad standard L2 routines.
