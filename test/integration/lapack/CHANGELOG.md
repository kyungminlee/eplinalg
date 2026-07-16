# tests/lapack — CHANGELOG

Resolved items, reverse-chronological. Open work lives in `TODO.md`.

## 2026-05-02 — L23 audit clean

The L23 audit (script below) returns only 2 deferred kernel-internal
routines (`dsb2st_kernels` / `zhb2st_kernels` — see TODO §sb2st
kernels). All other user-facing L6..L23 routines now ship a test
driver.

All 43 missing routines landed (37 in the initial closeout + 6
orbdb2/3/4 + unbdb2/3/4 closing the L15 caveat below):

| Stage | Phase | Routines landed |
|-------|-------|---|
| 1 | L18 | zlansy, zlansb, zlansp (z complex-symmetric); dlansf (RFP) |
| 1 | L19 | dlacon, zlacon (legacy 1-norm est) |
| 1 | L8  | dsytri_3, zsytri_3 |
| 1 | L13 | dtpqrt2, ztpqrt2, dtplqt2, ztplqt2 (and tpmlqt was already in tree) |
| 2 | L9  | dsytri2, zsytri2, dsytrs2, zsytrs2, dsytrs_3, zsytrs_3, dsycon_3, zsycon_3 |
| 3 | L15 | dorbdb1/5/6 + zunbdb1/5/6 + dorm22/zunm22 (initial pass: 8) |
| 3 | L15 | dorbdb2/3/4 + zunbdb2/3/4 (closeout: 6 — see entry below) |
| 4 | L6+L7 | dggsvd3, zggsvd3, dggsvp3, zggsvp3, dorcsd2by1, zuncsd2by1 |
| 5 | L23 | audit clean |

L17 (P17xx XX drivers) and L22 (DMD) were already shipped in earlier
phases. L10/L11/L12/L16/L21 turned out to be already covered when the
audit was actually performed — none of those phases needed new work.

Audit script:

```bash
awk 'NR>=89 && NR<=1107 && /^\| ✓ \|/ {gsub(/^\| ✓ \| /,""); split($0,a,/ \|/); print a[1]}' tests/RESULT.md \
  | while read r; do
      ls tests/lapack/*/test_${r}.f90 >/dev/null 2>&1 || echo "MISSING: $r"
    done
```

## 2026-05-02 — L15 orbdb2/orbdb3/orbdb4 + unbdb2/unbdb3/unbdb4

Bit-exact PASS on kind16, ~18 digits on kind10, ~31 digits on
multifloats for all 6 routines (m=12 and m=16 sweeps). Tests at
`tests/lapack/factorization/test_{d,z}{orbdb,unbdb}{2,3,4}.f90`.

Resolution had three load-bearing pieces (the original ~50% theta
divergence was misdiagnosed as a migration bug):

1. **Test bug (orbdb2/3/4)** — `theta` is declared dimension `Q` per
   the LAPACK API but only `min(P, M-P, M-Q)` entries are filled by
   the routine. The remaining entries are output as undefined. The
   original test compared all `Q` entries → uninit memory diff between
   the two calls → spurious ~50% relative error. Fix: each test now
   slices to the routine's actual write range: `theta(1:p)` /
   `theta(1:m-p)` / `theta(1:m-q)`.
2. **Upstream LAPACK 3.12.1 bug in {s,d,c,z}orbdb3/unbdb3** — the
   early `DROT/SROT/CSROT/ZDROT` call inside the `I=2..M-P` loop
   passes `LDX11` for *both* INCX and INCY:

   ```fortran
   CALL DROT( Q-I+1, X11(I-1,I), LDX11, X21(I,I), LDX11, C, S )
   ```

   The rotated row sits in X21, so INCY must be `LDX21`. With
   `LDX11 /= LDX21` (typical when M-P < P, the actual orbdb3 regime),
   the upstream form strides past the X21 buffer and corrupts heap.
   Manifests as `munmap_chunk(): invalid pointer` on the next free.
   Fix touches three places:
   - `recipes/lapack/source_overrides/{s,d,c,z}orbdb3.f` — patched
     bodies wired through the migrator's normal precision-substitution
     pipeline. Patch changes only the single `LDX11` → `LDX21`
     argument in each variant.
   - `recipes/lapack.yaml` — new `source_overrides:` block.
   - `tests/lapack/reflapack/overrides/{dorbdb3,zunbdb3}.f` — same
     patched bodies for the `reflapack_quad` reference path. Without
     this the upstream bug also corrupts heap on the reference side.
3. **Test drivers** — `tests/lapack/factorization/test_d{orbdb,
   orbdb3,orbdb4}.f90` and `test_z{unbdb2,unbdb3,unbdb4}.f90` (6
   files) sized to satisfy each routine's invariant: P=M/4 for orbdb2
   (P ≤ min(M-P, Q, M-Q)), P=3M/4 for orbdb3
   (M-P ≤ min(P, Q, M-Q)), Q=3M/4 for orbdb4
   (M-Q ≤ min(P, M-P, Q)).

| Target      | dorbdb2/3/4 + zunbdb2/3/4 (m=12,16) |
|-------------|--------------------------------------|
| kind16      | 12/12 PASS (bit-exact)               |
| kind10      | 12/12 PASS (~18 digits)              |
| multifloats | 12/12 PASS (~31 digits)              |

## P22 — Dynamic mode decomposition (4 routines) delivered

  d/z × gedmd, gedmdq

Tests landed in `tests/lapack/eigenvalue/test_{d,z}gedmd{,q}.f90` plus
helpers (`gen_orthogonal_quad` in `test_data.f90`, `sort_eig_lex_z` in
`compare.f90`), 4 ref interfaces in `ref_quad_lapack.f90`, and 4 fypp
wrappers in `target_lapack_body.fypp`. Layer-(a) spectrum + layer-(b)
differential residuals; layer-(c) reconstruction not implemented
(would be a future tightening). All four variants pass bit-exact on
all three targets after the reference-side override at
`tests/lapack/reflapack/overrides/zgedmd.f90`.

### P22 — complex DMD precision divergence (kind16 / multifloats)

Bit-exact PASS on kind16 / kind10 / multifloats after fixing the
**reference** path. The migrated `xgedmd` / `xgedmdq` were correct all
along; the divergence was a latent fp64 leak in the quad-promoted
reference.

Root cause: `external/lapack-3.12.1/SRC/zgedmd.f90:1006`

```fortran
W(i,j) = CMPLX(RWORK(N+i),ZERO,KIND=WP)*W(i,j)   ! WP = real64 = 8
```

Under `-freal-8-real-16`, gfortran 13 aliases `real(KIND=8)` *storage*
to 16 bytes, but the explicit `CMPLX(...,KIND=8)` intrinsic still
rounds through fp64 internally before storing into the (aliased)
16-byte slot. That bleeds ~1 ULP at fp64 into the Rayleigh quotient,
which the eigenvalue solver then exposes. The migrator rewrites
`WP = real64` to `WP = 16` so the migrated form `CMPLX(...,KIND=16)`
hits the truly-quad branch and is correct. `zgedmdq` inherits the
leak by calling `zgedmd` recursively.

Fix: a reference-only override at
`tests/lapack/reflapack/overrides/zgedmd.f90` that uses
`RWORK(N+i)*W(i,j)` (mathematically identical) to avoid the intrinsic.

| Target      | test_dgedmd | test_zgedmd | test_dgedmdq | test_zgedmdq |
|-------------|-------------|-------------|--------------|--------------|
| kind16      | 48/48 PASS  | 48/48 PASS  | 32/32 PASS   | 32/32 PASS   |
| kind10      | 48/48 PASS  | 48/48 PASS  | 32/32 PASS   | 32/32 PASS   |
| multifloats | 48/48 PASS  | 48/48 PASS  | 32/32 PASS   | 32/32 PASS   |

Test design notes (for future tightening):

  - Sizes for `*gedmdq` capped at `(12, 8)`: single-trajectory snapshot
    conditioning `κ(F) ∝ |λ_max/λ_min|^N` collapses fast.
  - Layer-(a) tolerance: `100 * n^2 * eps` for `*gedmd`,
    `10000 * n^2 * eps` for `*gedmdq` (trajectory data needs the bump).
  - Layer-(b) is the differential
    `‖res_got − res_ref‖∞ / max(‖A‖_F, eps)` rather than a
    relative-to-zero ratio.
  - Layer-(c) (reconstruction
    `‖Y − Z·diag(λ)·Z⁺·X‖_F / ‖Y‖_F`) is not implemented in v1.

## P17xx — Extra-precise iterative-refinement family (18 routines)

  d/z × gbrfsx, gbsvxx, gerfsx, gesvxx, porfsx, posvxx, syrfsx, sysvxx
  z × herfsx, hesvxx

18/18 PASS on kind10 / kind16 / multifloats (~17–19 digits on kind10,
~32–34 on kind16, ~30–32 on multifloats). Tests live at
`tests/lapack/linear_solve/test_{d,z}{gesvxx,gerfsx,gbsvxx,gbrfsx,
posvxx,porfsx,sysvxx,syrfsx}.f90` plus `test_zhesvxx.f90` /
`test_zherfsx.f90`. Each driver compares X / RCOND (and RPVGRW for the
SVXX drivers); BERR / ERR_BNDS_NORM / ERR_BNDS_COMP are computed on
both sides but not differenced (algorithm-internal estimates, not
deterministic across precision targets).

Resolution had three load-bearing pieces:

1. **Migrated lapack now PUBLIC-links migrated xblas**
   (`cmake/CMakeLists.txt`). The xx-family inside qlapack/elapack/
   mlapack reaches the migrated `qla_*_extended` /
   `ela_*_extended` / etc., which call the migrated XBLAS L2 entry
   points (`BLAS_qgemv_x_`, `BLAS_egemv_x_`, ...). Adding xblas to
   the lapack archive's interface link makes those symbols resolve
   transitively in every test executable that pulls in lapack.
2. **Quad-precision Fortran XBLAS bridge for reflapack_quad**
   (`tests/lapack/reflapack/refxblas_quad_bridge.f90`). The reference
   side `reflapack_quad` is built from the upstream Netlib LAPACK with
   `-freal-8-real-16`, so its `dla_*_extended.f` files call the
   upstream-named `blas_dgemv_x_` / `blas_zhemv_x_` / etc. Those names
   are NOT in qxblas (q-prefixed) and NOT in exblas/mxblas (different
   working precision). The bridge defines them as Fortran subroutines
   compiled with the same `-freal-8-real-16` ABI, forwarding to the
   quad-promoted refblas_quad routines (DGEMV/ZGEMV/DGBMV/ZGBMV/
   DSYMV/ZHEMV) for everything with a BLAS analogue and hand-coding
   the complex-symmetric matvec (no upstream BLAS analog). The
   Fortran-callable signature OMITS the C interface's leading `order`
   arg (the f2c bridge in
   `external/xblas-1.0.248/src/<routine>/BLAS_<routine>_x-f2c.c`
   hardcodes `blas_colmajor`); getting this wrong manifests as
   `xerbla DGEMV parameter 6 (LDA) had an illegal value` because all
   subsequent args shift by one slot.
3. **Wrapper module entries** in
   `tests/lapack/common/target_lapack_body.fypp` and matching
   interfaces in `tests/lapack/common/ref_quad_lapack.f90`. NPARAMS=0
   in every test driver — disables the PARAMS knob so the upstream
   default (refinement on, threshold 10.0, normwise-only error
   bounds) applies. The 4-call workspace footprint is `4*N` real for
   d-side, `2*N` complex + `2*N` real for z-side, except `zhesvxx`
   which needs `5*N` complex.

Build wiring: `cmake/CMakeLists.txt` makes `${LIB_PREFIX}lapack`
PUBLIC-link `${LIB_PREFIX}xblas` whenever both targets exist. No
opt-in build flag — XBLAS is always built and the bridge is small
(~280 lines of Fortran).

## P6/P7 — zgesvdq runtime crash

PASS on all three targets (kind10/kind16/multifloats), ~31-34 digits.

Root cause: same migrator bug as L5 — `REAL(x, wp)` /
`REAL(x, KIND=wp)` in upstream zgesvdq's body (e.g.
`SQRT(REAL(M, KIND=WP))` for column-norm scaling) was being wrapped to
a 2-arg structure constructor `real64x2(x, wp)` on multifloats, and
the analogous wrap on kind16 sometimes corrupted an `INFO` /
workspace-size argument because the second arg landed in the wrong
slot. Fix landed in `replace_intrinsic_calls` and
`replace_generic_conversions` (drop the kind-spec second argument
when wrapping to the multifloats constructor / when emitting
`REAL(x, KIND=...)` for kind targets).

| Target      | test_zgesvdq      |
|-------------|-------------------|
| kind16      | 2/2 PASS (~34d)   |
| kind10      | 2/2 PASS (~18d)   |
| multifloats | 2/2 PASS (~31d)   |

## L5 — *gesvj returns Infinity on multifloats target only

PASS on all three targets (kind10/kind16/multifloats), ~31 digits on
multifloats.

Root cause: upstream multifloats's Blue's-style scaling constants were
defined inverted-in-magnitude: `DD_SBIG = 1e+50` (a SCALE-UP factor)
and `DD_SSML = 1e-50` (a SCALE-DOWN factor) — opposite to LAPACK's
convention where `sbig << 1` (scale-down for big values) and
`ssml >> 1` (scale-up for small values). Compounding that, the
1e±100 / 1e±50 magnitudes only kept `(ax * sbig)**2` representable for
`ax` up to ~1e+104. dgesvj's protect-from-underflow up-scaling pass
routinely produces column norms up to ~`huge(double)/sqrt(N)`, so
`(ax * sbig)**2` overflowed to `+Inf` on every column.

Fix shape (touches three places):

1. `recipes/lapack/mf_helpers/la_constants_mf.f90` — `mtsml/mtbig/
   mssml/msbig` now compute the LAPACK formula at compile time from
   the leading-limb double's `radix/minexponent/maxexponent/digits`,
   instead of re-exporting upstream multifloats's `DD_*` (which were
   wrong). Fixes the LAPACK side via `USE LA_CONSTANTS_MF` aliases in
   mlassq / wlassq / etc.
2. `recipes/blas/source_overrides/{dnrm2,dznrm2}.f90` (new) — replace
   upstream BLAS's inline `parameter` declarations of
   `tsml/tbig/ssml/sbig` (which the migrator's nuke-rename pass
   substitutes with imports of broken `DD_*`) with renamed locals
   `btsml/btbig/bssml/bsbig` initialized at first call from the LAPACK
   formula (SAVE-vars init'd in an `if (.not. initialized)` block —
   multifloats's `**` operator is user-defined and not allowed in
   PARAMETER initializer expressions). The renamed identifiers also
   dodge the migrator's nuke-rename match.
3. `src/migrator/fortran_migrator.py` — minor migrator fix in
   `replace_intrinsic_calls` and `replace_generic_conversions`:
   `REAL(x, wp)` was being wrapped to `real64x2(x, wp)` (a 2-arg
   structure constructor — fails because `real64x2` has only one
   component, `limbs`). Now drops the kind-spec second argument so the
   result is `real64x2(x)`, which dispatches through the multifloats
   generic interface (e.g. `dd_from_int`).

| Target      | test_dgesvj      | test_zgesvj      |
|-------------|------------------|------------------|
| kind16      | 2/2 PASS         | 2/2 PASS         |
| kind10      | 2/2 PASS         | 2/2 PASS         |
| multifloats | 2/2 PASS (~31d)  | 2/2 PASS (~31d)  |
