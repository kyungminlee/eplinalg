# Differential precision tests for migrated ScaLAPACK

This directory extends the `tests/pblas/` scaffold to the **distributed**
extended-precision ScaLAPACK libraries (`${LIB_PAIR_PREFIX}scalapack` —
`qxscalapack`, `eyscalapack`, `mwscalapack`). Each test runs a ScaLAPACK
operation over a BLACS process grid, gathers the distributed result
back to rank 0, and compares against serial quad-precision Netlib
LAPACK applied to the full global matrix.

The structure mirrors `tests/pblas/` closely — the grid/distribution
helpers, per-target wrapper contract, JSON report format, and cyclic
link-group wiring are all inherited. This README only covers what is
*different* from the PBLAS suite.


## Testing model

Identical to PBLAS: every rank materializes the global problem at
REAL(KIND=16) from a shared seed, copies its block-cyclic slab into
a local array, calls `target_pdXXX` (which hands through to the
migrated `pq*` / `px*` entry point), then rank 0 gathers the result
and compares against serial `reflapack_quad` applied to the global
matrix.

Problems chosen to exercise block-cyclic edge cases at small sizes
without blowing out runtime — most drivers iterate
`n ∈ {32, 64, 96}` (or rectangular variants) with `mb = nb = 8`.
Per-driver shapes vary; see each ``test_*.f90`` for specifics.


## Routines covered

| Group          | Routines                                              |
| -------------- | ----------------------------------------------------- |
| linear_solve   | pdgesv, pdgetrf, pdgetrs, pdpotrf, pdpotrs            |
| factorization  | pdgeqrf                                               |
| eigenvalue     | pdsyev (eigenvalues only — `JOBZ='N'`)                |
| svd            | pdgesvd (singular values only)                        |
| auxiliary      | pdlange, pdlacpy, pdlaset                             |

Complex wrappers (`pxgesv`, `pxgeqrf`, `pxheev`) are exposed by the
target module for future expansion; no test drivers are provided yet.


## Reference strategy

Uses `reflapack_quad` (from `tests/lapack/reflapack/`) for the
LAPACK-level comparisons. When `_reflapack_src/` is not staged,
`reflapack_quad` degrades to the system `-llapack` (double precision),
which caps the meaningful precision at ~1e-16 regardless of target —
the tests still run, but the digit count reported will plateau at 15–16.

`refblas_quad` is also linked because the ScaLAPACK↔LAPACK↔BLAS cycle
has to resolve together under `LINK_GROUP(RESCAN)`.


## Dependencies

ScaLAPACK sits on top of the full stack, so every layer below must be
present. The CMake guards skip the suite (with a `STATUS` message) if
any of these is missing:

- `HAVE_REAL16`
- `MPI_Fortran_FOUND`
- `refblas_quad` (from `tests/blas/`)
- `reflapack_quad` (from `tests/lapack/`)
- `${LIB_PREFIX}scalapack`, `${LIB_PREFIX}pblas`, `${LIB_PREFIX}blacs`,
  `${LIB_PREFIX}lapack`
- A per-target wrapper directory `target_${TARGET_NAME}/` with
  `target_scalapack.f90` (`target_kind10/`, `target_kind16/`, and
  `target_multifloats/` are all populated)


## Layout

```
tests/scalapack/
├── CMakeLists.txt                    # guards, shim, wrapper, LINK_GROUP cycle
├── common/                           # copied verbatim from tests/pblas/common
│   ├── prec_kinds.f90
│   ├── compare.f90
│   ├── prec_report.f90               # module pblas_prec_report (rank-0 JSON)
│   ├── test_data.f90
│   ├── ref_quad_blas.f90
│   ├── ref_quad_lapack.f90           # copied from tests/lapack/common
│   ├── pblas_grid.f90                # BLACS init + numroc + descinit
│   └── pblas_distrib.f90             # block-cyclic scatter/gather
├── target_kind10/
│   └── target_scalapack.f90          # passthrough to pe*/py* routines
├── target_kind16/
│   └── target_scalapack.f90          # passthrough to pq*/px* routines
├── target_multifloats/
│   └── target_scalapack.f90          # passthrough to pt*/pv* routines
├── linear_solve/  test_pdgesv.f90, test_pdgetrf.f90, test_pdgetrs.f90,
│                  test_pdpotrf.f90, test_pdpotrs.f90
├── factorization/ test_pdgeqrf.f90
├── eigenvalue/    test_pdsyev.f90
├── svd/           test_pdgesvd.f90
└── auxiliary/     test_pdlange.f90, test_pdlacpy.f90, test_pdlaset.f90
```

A dedicated `Fortran_MODULE_DIRECTORY`
(`${PROJECT_BINARY_DIR}/fmod/${FORTRAN_MOD_COMPAT_TAG}/scalapack_tests`)
prevents `.mod` collisions with `tests/pblas/` — the helper modules
reuse the same names (`pblas_grid`, `pblas_distrib`, `compare`, …)
because renaming them would defeat the point of reusing the sources.


## Caveats and bring-up notes

The original draft of this suite shipped with three local workarounds
(`scalapack_test_blacs_c_api`, `scalapack_test_pilaenv`,
`scalapack_test_shim`, plus `target_kind16/{lamch_shim,ptzblas_stubs}.f90`
and `common/pilaenv.f`). All have been retired; the migrated stack now
provides the underlying symbols directly. A few non-obvious things
surfaced during that work and are worth keeping written down:

1. **REDIST/SRC needs file-stem-based clone naming.** Files like
   `pdgemr.c` / `pdgemr2.c` / `pdtrmr2.c` expose their Fortran-callable
   entry points only behind ``#define fortran_mr2dnew pdgemr2d_``
   macros, so the symbol scanner never sees `pdgemr2d_` and
   `c_migrator` can't find a rename-map entry for the file stem. The
   migrator now falls back to `_redist_clone_stem`, which substitutes
   the precision letter at position 1 of the stem (`pdgemr2` →
   `pqgemr2` at kind16) using `target_mode.prefix_map`. Without this,
   only `pdgemr.c` clones (matching the first `#define`'s `pdgemr2do_`
   target) and the actual `Cpdgemr2d` implementation in `pdgemr2.c` is
   silently dropped — pdgesvd then heap-corrupts because the migrated
   call ends up freeing kind16 buffers through the original
   double-precision allocator.

2. **`dge*` helper renames are recipe-driven.** The same `#define`
   pattern (`#define freememory dgefreememory`) hides
   `dgescanD0` / `dgedispmat` / `dgesetmemory` / `dgefreememory` /
   `dgescan_intervals` / `Cdgelacpy` from the symbol scanner.
   `codegen/recipes/scalapack_c.yaml`'s `c_type_aliases` block spells out the
   per-family substitution; do not rely on the rename map alone.

3. **BLACS `BUFFALIGN=8` must become 16 for kind16 / multifloats.**
   `BI_qvmcopy` segfaults on the first reduction-buffer write (any
   `qgsum2d_` path — `pdlange`, `pdsyev`, `pdgesvd`, `pdgeqrf` all
   trip it) because `__float128` requires 16-byte alignment. Patched
   directly in `extern/scalapack-2.2.3/BLACS/SRC/Bdef.h`; no-op for
   double callers.

4. **`pdsyev` / `pdgesvd` need a workspace query.** The naive heuristics
   the test draft used (`8*N + 2*max(locm,locn)` and friends)
   undersize `LWORK` and trigger ``parameter number 14 had an illegal
   value`` / heap corruption respectively. Both drivers now do an
   `LWORK = -1` query, read `WORK(1)`, allocate, then call again.

5. **`-DAdd_` is mandatory for REDIST.** The Fortran-callable entry
   points in REDIST/SRC use the upstream `Add_` / `UpCase` /
   `NoChange` switch to pick a name-mangling convention; without an
   explicit `target_compile_definitions(... PRIVATE Add_)` the symbols
   come out as `pqgemr2d` (no trailing underscore) and gfortran can't
   resolve them. Wired in `cmake/CMakeLists.txt` on the C OBJECT libs
   that fold into the `${LIB_PREFIX}scalapack` archive.

6. **Multifloats scalapack C clones are in.** REDIST/SRC originally shipped
   K&R-style function definitions that didn't survive C-as-C++
   compilation. They've all been ANSI-fied in place
   (`extern/scalapack-2.2.3/REDIST/SRC/p*gemr*.c`, `p*trmr*.c`),
   and the `target_multifloats/target_scalapack.f90` wrapper is now
   populated. The library now builds and tests pass under the
   multifloats target — see `tests/REPORT.md` for the per-routine
   digits-of-agreement.

7. **`sizeof(double)` in C migrator promotion.** The
   `_apply_c_type_subs` cast-protection regex `\(\s*double\s*\)`
   originally over-matched `sizeof(double)`, so REDIST allocations like
   `mr2d_malloc(blocksize * sizeof(double))` stayed at 8-byte sizing
   while the buffer was cast to a 16-byte `float64x2*`. PDGEMR2D's
   recvbuff overran 2× and crashed inside `pdsyev` / `pdgesvd`
   teardown. Fixed by adding a negative lookbehind for `sizeof` /
   `alignof` in the protect regex (commit `ef3bf81`).

8. **`DDDOT` symbol collision drove the prefix switch away from
   `DD`/`ZZ`; multifloats now uses single-letter `M`/`W`.**
   ScaLAPACK ships its own `DDDOT` subroutine in `TOOLS/dddot.f`.
   With the original `DD`/`ZZ` multifloats prefixes, BLAS `DDOT`
   renamed to `DDDOT` and collided with the orphaned ScaLAPACK
   wrapper — the linker silently picked the SUBROUTINE form,
   corrupting `pddpotf2`'s diagonal during panel Cholesky. The
   migrator iterated through a transient `T`/`V` step before
   settling on `M`/`W` as the final single-letter pair (no upstream
   BLAS / LAPACK / ScaLAPACK symbol begins with either letter),
   eliminating the collision class entirely. The migrator now also
   emits a stderr warning when an orphaned symbol's name shadows
   another symbol's rename target (`prefix_classifier.build_rename_map`).


## Running

Same as PBLAS — configure a staged target (`uv run python -m migrator
stage <dir> --target kind16`), build, and run with `mpirun`:

```
cmake -S <stage> -B <stage>/build
cmake --build <stage>/build --target test_pdpotrf -j
ctest --test-dir <stage>/build -R '^scalapack_test_pdpotrf$' --output-on-failure
```

Default `SCALAPACK_TEST_NPROC=4` (2×2 grid). Override at configure time
with `-DSCALAPACK_TEST_NPROC=<n>`.
