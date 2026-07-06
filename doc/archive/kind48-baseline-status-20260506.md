# kind4 / kind8 baseline columns — status

> **Archived snapshot** (moved from `tests/KIND48_BASELINE_STATUS.md`).
> Describes the state as of 2026-05-06 and is not kept up to date;
> current baseline coverage is whatever `migrator stage --target
> kind{4,8}` + ctest reports today.

*Status as of 2026-05-06.* The plumbing is in place and produces sane numbers,
but the baseline column is **not yet promoted** into the published precision
matrix — residual upstream issues (catalogued under "Recommended treatment"
below) keep it parked until the parked items are addressed or explicitly
accepted.

## Current state

The plumbing is fully wired and the kind4 / kind8 baseline target
produces sane numbers from a `migrator stage --target kind8` tree:

- `migrator stage --target kind8 /tmp/x` short-circuits the migration
  step and writes a `target_config.cmake` with empty `LIB_PREFIX` so
  the unified CMake build resolves the standard archives (`blas`,
  `lapack`, `pblas`, `scalapack`) as the libraries under test.
  Same shape for `--target kind4`.
- `cmake -S /tmp/x -B /tmp/x/build && cmake --build /tmp/x/build -j8`
  links cleanly. `ctest` runs the framework subsets that the baseline
  supports (BLAS + the non-XBLAS LAPACK tests, plus PBLAS / ScaLAPACK
  once those are wired analogously).

Smoke results on the kind8 BLAS + LAPACK suite (725 tests):

| Suite  | Pass | Fail | Notes |
|--------|-----:|-----:|-------|
| BLAS   |   75 |    0 | Every Level 1/2/3 routine reports ~16-digit agreement against the quad reference (the IEEE binary64 floor — exactly the expected sanity-check shape). |
| LAPACK |  720 |    5 | After the symbol-rename pass was extended to gfortran module-procedure symbols (`__<mod>_MOD_<proc>`) AND `refquad_alias.py` was taught to recognize typed function declarations (`logical function disnan(...)`), the SIGSEGV cluster and the `disnan` ABI mismatch cleared. The 5 residual failures split into 2 traceable to a missing baseline application of `recipes/lapack/source_overrides/{d,z}orbdb3.f` (parked in `tests/lapack/TODO.md`) and 3 genuine precision findings. |

`kind16` regression check stays clean — `test_dasum`, `test_dgemm`,
`test_dsyevr`, and `test_dgbsvxx` all pass after the rename and
exclude-list changes.

## What was needed (and now exists)

1. **Per-target shims.** `tests/{blas,lapack,pblas,scalapack}/target_kind{4,8}/`
   fypp shims setting `prefix_real='s'/'d'`, `prefix_complex='c'/'z'`,
   etc. Per-target `tests/common/target_kind{4,8}/target_conv.f90`
   conversion modules (trivial real(4)/(8)↔real(16) intrinsic casts).
2. **Baseline staging path.** `migrator stage --target kind4|kind8`
   skips migration entirely and writes a target_config.cmake with
   `LIB_PREFIX=""`. `cmake/CMakeLists.txt` gates every
   `add_migrated_*` and the install loop on a `BASELINE_BUILD` flag
   so the std archives stand alone. Each `tests/<lib>/CMakeLists.txt`
   has a baseline branch that links the bare std archive instead of
   `${LIB_PREFIX}<lib>`.
3. **Reference symbol rename.** `librefblas_quad.a` and
   `libreflapack_quad.a` get every Fortran-mangled symbol they
   define renamed to `*_quad_` via a POST_BUILD `objcopy
   --redefine-syms` step, so the un-renamed std archives can coexist
   on the same link line. The rename is driven by
   `scripts/refquad_rename_archive.sh`:
   - Step 1 (scan) — each archive's POST_BUILD writes a
     `.quad-symbols` snapshot of its T+W defined Fortran symbols
     *before* renaming. A sibling rename invoked later reads that
     snapshot rather than nm-ing the (already renamed) archive.
   - Step 2 (rename) — the redefine list is the union of own +
     sibling snapshots, minus the T+W defs of any `-x` exclude
     archive. The exclude list pulls in the migrated `*_common`
     archives (`lapack_common`, `mumps_common`, …) so their
     precision-independent helpers (iparam2stage_, la_constants_ep,
     mumps_abort_) stay un-renamed and resolve normally in the
     migrated link cycle.
   - `ref_quad_*.f90` modules are rewritten by
     `scripts/refquad_alias.py`: each routine gets the original
     interface declaration suffixed `_quad` plus a thin module-
     procedure wrapper under the original name. Test files keep
     their existing `use ref_quad_blas, only: dasum` calls — gfortran
     emits those as module-procedure references
     (`__ref_quad_blas_MOD_dasum`), no longer colliding with the
     standard archive's `dasum_`.
4. **XBLAS baseline stubs.** Std lapack's iterative-refinement
   helpers (dla_*_extended.f, zla_*_extended.f) reference XBLAS
   entry points (blas_dgemv_x_, …). Those .o files get pulled into
   baseline tests transitively via `target_lapack`'s public xx /
   rfsx wrappers, even though those tests themselves are filtered
   out of the baseline ctest set. `tests/lapack/reflapack/xblas_baseline_stubs.c`
   provides aborting C stubs for the entry points; the archive is
   linked with `--whole-archive` so ld extracts every stub
   regardless of demand timing.
5. **Migrated cycle additions.** `tests/mumps/CMakeLists.txt` adds
   std `blas` and `lapack` to its archive cycle so precision-
   independent helpers (lsame_, ilaenv_, dlamch_) resolve from a
   non-renamed source. `tests/pblas/CMakeLists.txt` adds std
   `scalapack` to the kind4 / kind8 cycle for the zzdotc_ /
   zzdotu_ TOOLS routines that PTZBLAS calls.

## Known residual issues

### LAPACK SIGSEGV cluster (kind8) — RESOLVED

The earlier cluster of CSD / 2-stage / DMD SIGSEGVs (dorbdb1-4,
dorcsd2by1, dgedmd, dgges, zheevr_2stage, …) traced to module-proc
symbol collisions: la_xisnan.F90 and la_constants.f90 are compiled
into both std lapack (native) and reflapack_quad (quad-promoted),
producing identical `__la_xisnan_MOD_disnan` / `__la_constants_MOD_dp`
symbols with mismatched ABIs. Whichever .o ld picked first determined
the receiver's argument width, breaking the other path's callers.

`scripts/refquad_rename_archive.sh` was extended to match
`__<mod>_MOD_<proc>` alongside the trailing-underscore Fortran-mangled
pattern. After the rename:
- std lapack continues to expose `__la_xisnan_MOD_disnan` (native) for
  its own native callers (`dlassq_`, `dorbdb5_`, …).
- reflapack_quad exposes `__la_xisnan_MOD_disnan_quad` (quad) for its
  own quad-promoted callers.
- Both paths resolve consistently; the cluster is gone.

### Residual failures — categorized

`disnan` resolved (2026-05-06): `refquad_alias.py`'s `SIG_START_RE`
was untyped-function-only (`\s+(function|subroutine)\s+`), so the
typed declaration `logical function disnan(din)` was skipped during
the alias rewrite — `disnan` was never renamed to `disnan_quad`, no
module-procedure wrapper was emitted, and `use ref_quad_lapack, only:
disnan` resolved to std lapack's native `disnan_` (real-8 ABI) called
with a quad-precision argument. Cases 5 and 7 happened to expose the
ABI mismatch as wrong booleans; the original "NaN-payload artifact"
explanation was incorrect. Fix: extended the regex to accept an
optional return-type prefix and taught `render_wrapper` to preserve
the typed-function form (no `result()` clause, return via assignment
to the function name). Regenerated `tests/lapack/common/ref_quad_lapack.f90`.

`zunbdb3` and `dorbdb3` re-categorized (2026-05-06): both trace to
the same upstream LAPACK 3.12.1 typo in `{d,z}orbdb3.f` —
`{D,Z}ROT(..., X21(I,I), LDX11, ...)` should pass `LDX21` as the
`INCY` for `X21`. When `LDX11 /= LDX21` the rotation strides past
`X21`'s buffer and corrupts the next heap chunk. The fix exists in
`recipes/lapack/source_overrides/{d,z}orbdb3.f` and is applied for
every migrated target, but the kind4 / kind8 baseline path skips
migration and stages upstream `_reflapack_src/` verbatim. With the
override overlaid, both `dorbdb3` and `zunbdb3` pass cleanly.
Parked under `tests/lapack/TODO.md` until the baseline column is
promoted to the published precision matrix.

| # | Test | Class | Severity |
|---|------|-------|----------|
| 1 | `lapack_test_zunbdb3` | upstream LAPACK orbdb3 typo (heap overrun) | parked — migrator/staging gap |
| 2 | `lapack_test_dorbdb3` | upstream LAPACK orbdb3 typo (silent corruption presents as wrong digits) | parked — migrator/staging gap |
| 3 | `lapack_test_dorcsd`  | precision genuinely insufficient (orthogonal-completion ambiguity at this shape) | finding |
| 4 | `lapack_test_dgedmdq` | precision insufficient on `JOBR=R` residual-refinement path | finding |
| 5 | `lapack_test_zgedmdq` | precision insufficient on `JOBR=R` residual-refinement path | finding |

#### 1, 2 — `zunbdb3` / `dorbdb3` (upstream LAPACK orbdb3 typo)

```
zunbdb3:  SIGSEGV in arena_for_chunk → __libc_free
            __target_lapack_MOD_target_zunbdb3 (auto-deallocate)
dorbdb3 [m=12]  max_rel_err = 4.86E-02   digits = 1.31  FAIL
dorbdb3 [m=16]  max_rel_err = 1.03E-01   digits = 0.99  FAIL
```

Same root cause for both. Upstream `{d,z}orbdb3.f` has a typo in the
`I .GT. 1` branch's `{D,Z}ROT` call: passes `LDX11` as `INCY` for
`X21`, where `INCY` must be `LDX21`. When `LDX11 /= LDX21` (always
true in the orbdb3 regime since `M-P < P`), the rotation strides past
`X21`'s allocated extent and writes into whatever heap chunk follows.

Heap-layout chance decides the manifestation:
- `zunbdb3` happens to corrupt the chunk metadata of an array freed
  during the wrapper's auto-deallocate, surfacing as SIGSEGV in
  `__libc_free`.
- `dorbdb3` corrupts data words of an adjacent buffer (typically
  `theta` or `phi` slack tail), and the test sees `theta_g` / `phi_g`
  with stomped values, producing the `~1` digit agreement.

The migrator already carries the one-character fix in
`recipes/lapack/source_overrides/{d,z}orbdb3.f` (see
`recipes/lapack.yaml`'s `source_overrides` block) and applies it for
every migrated target. The kind4 / kind8 baseline path skips
migration and stages upstream `_reflapack_src/` verbatim, so the std
`lapack` archive linked under test still has the typo. Confirmed
both failures clear when the override is overlaid into
`_reflapack_src/`.

Parked in `tests/lapack/TODO.md`. Until promoted, suppress these two
from the published baseline column.

#### 3 — `dorcsd` (full real CSD assembly)

```
dorcsd [m=12]  max_rel_err = 2.92    digits = -0.47  FAIL
dorcsd [m=16]  max_rel_err = 8.70E+1 digits = -1.94  FAIL
```

`dorcsd` is the user-facing CSD driver that calls `dorbdb*` underneath
and then assembles the orthogonal factors. Independent of the
`dorbdb3` typo: re-tested with `recipes/lapack/source_overrides/`
overlaid into the baseline `_reflapack_src/` and `dorcsd` still
fails with the same negative digit counts. The kind8 result and the
quad reference end up picking different canonical representatives
among the equivalent orthogonal completions for these problem shapes.

Strongest signal in the table that **CSD at double precision is
unreliable for general shapes**. For routines like `dorcsd` the
cell should probably show `(unstable)` rather than a raw digit count.

#### 4, 5 — `dgedmdq` / `zgedmdq` (residual-refined DMD)

```
dgedmdq [m=12,n=8,JOBR=R,layer=b] max_rel_err = 4.03E-03  digits = 2.40  FAIL
zgedmdq [m=12,n=8,JOBR=R,layer=b] max_rel_err = 4.93E-02  digits = 1.31  FAIL
```

Every other case combination passes (`14-15` digits — normal kind8
floor). Failures are isolated to `JOBR=R` (Rayleigh-Ritz residual
refinement) at `m=12, n=8, layer=b` (the second test layer with a
specific input distribution).

The residual-refinement step subtracts two nearly-equal quantities
to compute Ritz residuals; in the quad reference those quantities
are computed with ~33 digits of headroom and the cancellation
cleanly exposes the small residual. In kind8 native, the cancellation
loses most of the digits and the residual ends up dominated by
roundoff. The quad-vs-native discrepancy is real but expected — it's
exactly the precision phenomenon DMD-Q's residual mode is *designed*
to surface.

For the kind4/kind8 columns this should probably be reported per
parameter combination rather than as a single cell; users picking a
double-precision DMD should know `JOBR=R` carries this caveat.

### Recommended treatment in `tests/RESULT.md`

When the kind4 / kind8 columns are turned on for the public report:

- `zunbdb3`, `dorbdb3` — suppress until the staging path overlays
  `recipes/lapack/source_overrides/` on top of `_reflapack_src/`.
  See `tests/lapack/TODO.md`.
- `dorcsd` — show the digit count honestly. The negative values
  are the headline finding (CSD orthogonal-completion ambiguity at
  double precision).
- `dgedmdq`, `zgedmdq` — consider reporting per-`JOBR` or annotating
  the cell with `(JOBR=R: 2.4)`; the aggregate worst-case is
  unrepresentative.

### Other framework coverage (PBLAS / ScaLAPACK)

The PBLAS / ScaLAPACK kind8 cycles are wired (`scalapack` added
to the PBLAS cycle, std blas/lapack added to the mumps cycle), but
not yet smoke-tested. The same module-procedure rename now in place
should cover them; expect runtime results in line with the LAPACK
suite (≥99% pass for the routines that don't compound CSD-style
near-degenerate sub-blocks).
