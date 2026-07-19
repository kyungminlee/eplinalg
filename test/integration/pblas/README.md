# Differential precision tests for migrated PBLAS

This directory extends the `tests/blas/` scaffold to the **distributed**
extended-precision PBLAS libraries (`${LIB_PAIR_PREFIX}pblas` —
`qxpblas`, `eypblas`, `mwpblas`). Each test runs a PBLAS operation over
a BLACS process grid, gathers the distributed result back to rank 0,
and compares against serial quad-precision BLAS applied to the full
global matrix.

The structure mirrors `tests/blas/` closely — the core testing model
(reference vs. target, tolerance formula, JSON report, per-target
wrapper contract) is unchanged. Everything below describes only the
things that are *different* because the routines under test are
distributed.


## Testing model

Each test runs in the following order (on every MPI rank):

1. `grid_init()` — BLACS + MPI initialization, picks a 2D grid shape
   as close to square as possible for the given process count
   (2×2 for NP=4, 1×1 for NP=1, 1×N for prime NP).
2. For each problem shape:
   - Every rank independently materializes the **full global**
     vector/matrix at REAL(KIND=16) using a deterministic seed — so
     ranks don't need to scatter, they each read their own
     block-cyclic piece from their own mirrored copy.
   - The per-target wrapper (`target_pdgemm` etc.) converts the local
     arrays to the target's precision, calls the migrated PBLAS
     routine (`pqgemm` / `pegemm` / `pddgemm`), and converts the
     in/out buffers back to quad.
   - Rank 0 gathers the distributed result from every process via
     point-to-point `MPI_BYTE` transfers (MPI_REAL16 is not portable).
   - Rank 0 computes the reference by running the serial
     **quad-precision Netlib BLAS** (`refblas_quad`, built by
     `tests/blas/refblas/`) on the full global matrix.
   - Rank 0 computes `max_rel_err` and calls `report_case`; non-zero
     ranks skip it.
3. `report_finalize()` is collective — rank 0 broadcasts the pass/fail
   flag so every rank exits with the same code and ctest observes a
   uniform failure/success across the entire `mpirun` invocation.
4. `grid_exit()` tears down BLACS + MPI.

Tolerance still follows `safety_factor * n_ops * target_eps`, but the
safety factor is doubled over the BLAS suite (16 → 32 for most cases,
64 for triangular solves) to absorb the small amount of extra rounding
that comes from SUMMA-style summation orders in PBLAS kernels.


### Reference strategy

PBLAS is a C library, so the `-freal-8-real-16` promotion that builds
`refblas_quad` from Netlib Fortran sources does not apply. Rather than
porting 122 C files to `__float128` plus MPI types, this test suite
uses the **serial** quad-precision BLAS as a differential reference:
the migrated distributed routine at target precision should agree with
serial BLAS at quad precision up to rounding and SUMMA summation order.

Consequence: `tests/blas/` must also be staged for PBLAS tests to
build — the CMakeLists gates on `TARGET refblas_quad` and skips
gracefully otherwise.


## Directory layout

```
tests/pblas/
├── README.md                     ← you are here
├── CMakeLists.txt                ← gated build, registers ctest entries
├── ../common/                    ← shared harness modules (compiled per
│   │                               suite into pblas' own moddir)
│   ├── prec_kinds.f90            ← ep=16, dp=8 (shared with tests/blas)
│   ├── compare.f90               ← rel-error helpers (shared)
│   ├── test_data.f90             ← seeded RNG (shared — same seed → same
│   │                               data on every rank, no scatter needed)
│   ├── prec_report.F90           ← rank-0-only JSON writer (built with
│   │                               PREC_REPORT_MPI here)
│   ├── ref_quad_blas.f90         ← serial quad-BLAS interfaces
│   ├── pblas_grid.F90            ← BLACS init/exit + grid shape picker
│   │                               + inline numroc / descinit / g2l
│   └── pblas_distrib.fypp        ← block-cyclic scatter (mirror → local)
│                                   + gather (local → rank-0 global);
│                                   fypp: real + complex from one body
├── common/target_pblas_body.fypp ← shared wrapper template body
├── target_kind10/target_pblas.f90
├── target_kind16/target_pblas.f90
├── target_multifloats/target_pblas.f90
├── level1/                       ← 18 Level 1 PBLAS tests
│   ├── test_pddot.f90            ← scalar out
│   ├── test_pdnrm2.f90
│   ├── test_pdasum.f90
│   ├── test_pdamax.f90           ← (amax,indx) tuple out
│   ├── test_pdscal.f90           ← in-place vector
│   ├── test_pdaxpy.f90
│   ├── test_pdcopy.f90           ← out-of-place vector
│   ├── test_pdswap.f90
│   ├── test_pdzasum.f90          ← mixed real-from-complex
│   ├── test_pdznrm2.f90          ← mixed real-from-complex
│   ├── test_pzdotc.f90           ← complex
│   ├── test_pzdotu.f90
│   ├── test_pzaxpy.f90
│   ├── test_pzcopy.f90
│   ├── test_pzswap.f90
│   ├── test_pzscal.f90
│   ├── test_pzdscal.f90          ← mixed complex-with-real-scalar
│   └── test_pzamax.f90
├── level2/                       ← 15 Level 2 PBLAS tests
│   ├── test_pdgemv.f90
│   ├── test_pdger.f90            ← rank-1 update
│   ├── test_pdsymv.f90
│   ├── test_pdsyr.f90            ← symmetric rank-1
│   ├── test_pdsyr2.f90           ← symmetric rank-2
│   ├── test_pdtrmv.f90
│   ├── test_pdtrsv.f90           ← triangular solve (diag-dominant)
│   ├── test_pzgemv.f90
│   ├── test_pzhemv.f90           ← Hermitian mat-vec
│   ├── test_pzgerc.f90           ← conjugated rank-1
│   ├── test_pzgeru.f90           ← unconjugated rank-1
│   ├── test_pzher.f90            ← Hermitian rank-1 (alpha real)
│   ├── test_pzher2.f90           ← Hermitian rank-2 (alpha complex)
│   ├── test_pztrmv.f90
│   └── test_pztrsv.f90           ← complex tri solve (diag-dominant)
└── level3/                       ← 15 Level 3 PBLAS tests
    ├── test_pdgemm.f90
    ├── test_pdsymm.f90
    ├── test_pdsyrk.f90
    ├── test_pdsyr2k.f90
    ├── test_pdtrmm.f90
    ├── test_pdtrsm.f90
    ├── test_pzgemm.f90
    ├── test_pzhemm.f90
    ├── test_pzsymm.f90
    ├── test_pzsyrk.f90
    ├── test_pzsyr2k.f90
    ├── test_pzherk.f90
    ├── test_pzher2k.f90
    ├── test_pztrmm.f90
    └── test_pztrsm.f90
```


## Per-target wrapper contract

Every `target_<target>/target_pblas.f90` exposes the same module
`target_pblas` with:

```fortran
module target_pblas
    public :: target_name, target_eps
    public :: target_pddot, target_pdgemm, …  ! one per routine
end module
```

- `target_name`: used in JSON filename (`<routine>.<target_name>.json`).
- `target_eps`: target's machine epsilon for tolerance sizing.
- Each wrapper takes REAL(KIND=16) / COMPLEX(KIND=16) local
  distributed arrays + unchanged descriptor/index arguments, casts
  element-wise to the target native type, calls the migrated PBLAS
  routine, and casts the in/out buffer back to quad. Descriptors are
  distribution info — they pass through unchanged regardless of target
  precision.

The migrated name scheme (scalapack-style single-character slot
substitution):

| original   | kind10     | kind16     | multifloats   |
|------------|------------|------------|---------------|
| `pdgemm`   | `pegemm`   | `pqgemm`   | `pddgemm`     |
| `pddot`    | `pedot`    | `pqdot`    | `pdddot`      |
| `pzgemm`   | `pygemm`   | `pxgemm`   | `pzzgemm`     |
| `pzherk`   | `pyherk`   | `pxherk`   | `pzzherk`     |
| `pdzasum`  | `peyasum`  | `pqxasum`  | `ptvasum`     |
| `pdznrm2`  | `peynrm2`  | `pqxnrm2`  | `ptvnrm2`     |
| `pzdscal`  | `pyescal`  | `pxqscal`  | `pvtscal`     |

The mixed-precision rows (last three) substitute both the real and
complex slots independently. `target_pblas_body.fypp` derives those
prefixes from `prefix_real` / `prefix_complex` automatically — the
per-target shim files only need to set the two single-letter prefixes.


## Running

```bash
python -m migrator stage /tmp/staging-kind16 --target kind16 \
       --libraries blas blacs ptzblas pbblas pblas
cmake -S /tmp/staging-kind16 -B /tmp/staging-kind16/build \
      -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/staging-kind16/build -j8
ctest --test-dir /tmp/staging-kind16/build -R '^pblas_' --output-on-failure
python scripts/precision_report.py /tmp/staging-kind16/build/precision_reports/
```

The process count is fixed at configure time via `-DPBLAS_TEST_NPROC=<N>`
(default 4 → 2×2 grid). `-DPBLAS_TEST_NPROC=1` runs every test in a
degenerate 1×1 grid, useful as a sanity check.

Staging **must** include `blas`, `blacs`, `ptzblas`, `pbblas`, and
`pblas`. `cmake/CMakeLists.txt` wires `${LIB_PREFIX}pblas` to depend
only on `blas` + `blacs`, but the migrated pblas C code calls into
PTZBLAS (e.g. `qvasum_`) and PBBLAS kernels at runtime — the tests
link those two libraries explicitly to satisfy the references. Without
the `blas` library staged, `refblas_quad` isn't built and the PBLAS
tests skip.


## Reading the report

The JSON output uses the same filename convention as the BLAS suite
(`<routine>.<target>.json` in `<build>/precision_reports/`) and
`scripts/precision_report.py` scans that directory without any
PBLAS-specific changes — PBLAS rows appear in the same summary table
as BLAS rows.

| routine | kind10 | kind16 | multifloats |
|---------|--------|--------|-------------|
| ddot    | 17.59  | exact  | 30.68       |
| pddot   | ~17.2  | exact  | ~29.8       |
| pdgemm  | ~17.0  | exact  | ~29.5       |

For kind16 PBLAS rows, "exact" is the expected outcome when the
tolerance comfortably absorbs SUMMA summation-order differences —
reference and target are both at REAL(KIND=16), and the migrated
distributed algorithm differs from serial BLAS only in associativity.


## Adding a new test

For a PBLAS routine `pxroutine` already wrapped in
`common/target_pblas_body.fypp` and declared in
`../common/ref_quad_blas.f90`:

1. Drop a new program into the appropriate `level{1,2,3}/` directory
   (one of the existing `test_pd*.f90` files makes a suitable
   starting template).
2. `CMakeLists.txt` picks up `test_*.f90` via
   `file(GLOB CONFIGURE_DEPENDS …)`, so no manual registration is
   needed.

For a routine **not** yet declared:
1. Add an explicit interface to `../common/ref_quad_blas.f90` for the
   corresponding serial reference routine.
2. Add an abstract interface and a concrete `target_p<routine>`
   wrapper inside `common/target_pblas_body.fypp` (one fypp file
   shared across all targets) — signature is uniformly
   REAL(KIND=16) / COMPLEX(KIND=16) on the outside; the body uses
   `q2t_r` / `q2t_c` to convert into the target's native type before
   the migrated PBLAS call. For mixed-precision routines (`pdz*` or
   `pzd*`) use `${RC}$` / `${CR}$` instead of `${R}$` / `${C}$`.
3. Export the new wrapper from the `public ::` list at the top of
   the same fypp file.


## Open issues

See `TODO.md` for follow-ups whose fix lies outside the
`tests/pblas/` subtree.

### Multifloats PBLAS C linkage

The multifloats target compiles migrated C sources as C++ (to pick up
`float64x2` operator overloading). PBLAS's Fortran-callable entry
points (`pddgemm_`, `pzzherk_`, …) are defined without explicit
`extern "C"`, so under C++ compilation every entry point is
name-mangled and `target_pblas.f90` can't resolve them at link time.
Running `ctest` on a multifloats staging build will therefore skip
`tests/pblas` (`CMakeLists.txt` returns early when `C_AS_CXX` is on).

Fixing this requires lockstep `extern "C"` wrappers on both the
function declarations (in `PBpblas.h` / `PBtools.h`) and the per-file
entry-point definitions generated by the C migrator. Until that lands,
the kind10 and kind16 targets exercise the migrated PBLAS code paths
in full.


