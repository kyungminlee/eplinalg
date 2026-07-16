# Test suites

This directory holds runtime correctness and performance tests for the
libraries produced by eplinalg. They exercise the *outputs*
of migration — the compiled extended-precision libraries — not the
migration engine itself (engine tests live under `test/unit/`).

Tests here are all of one form: **migrated vs Netlib reference compiled
at REAL(KIND=16)** (the `{blas,lapack,blacs,pblas,pbblas,scalapack,mumps,xblas,ptzblas}/` subdirectories here).
Differential precision: run the same operation through the migrated
extended-precision library and through Netlib source compiled with
`-freal-8-real-16`, compare in `REAL(KIND=16)`. JSON reports under
`<build>/precision_reports/`, aggregated by
`scripts/precision_report.py`.

The overlay-vs-migrated test scheme (epblas-parallel, epblas-openblas) moved
to the [epblas-parallel](../../epblas-parallel/) sibling repo along with the
overlays themselves.

## Migrated vs Netlib reference @ KIND=16 (differential precision)

Tests the *migrated* libraries against the original Netlib source
compiled with `-freal-8-real-16`. All comparisons happen in
`REAL(KIND=16)` regardless of which migrated target is under test.

Test data is generated at `KIND=16`, converted down to the target's
precision before calling the migrated routine, converted back to
`KIND=16` after, and compared against the reference call which stayed
at `KIND=16` throughout.

```
tolerance = safety_factor * n_ops * target_eps
pass      = max_rel_err <= tolerance
```

| Domain | Test files | Reference | Notes |
|---|---|---|---|
| `blas/` | 25 L1 + 33 L2 + 17 L3 = **75** | `refblas_quad` (Netlib BLAS @ KIND=16, built from `_refblas_src/`) | Canonical scheme; mirrored by every other domain |
| `xblas/` | 15 L1 + 14 L2 + 15 L3 = **44** | Netlib XBLAS @ extended precision | Tests the `BLAS_qgemv_x` etc surface |
| `lapack/` | **669** under `auxiliary/`, `eigenvalue/`, `factorization/`, `linear_solve/` | `reflapack_quad` | Largest domain |
| `ptzblas/` | **35** | tz-prefix triangle/packed-zero BLAS variants | |
| `blacs/` | **22** | BLACS ref @ KIND=16 | Distributed; uses MPI |
| `pblas/` | **62** across `level1/`, `level2/`, `level3/` | serial quad BLAS on global matrix | Distributed; BLACS process grids |
| `pbblas/` | **21** | PBBLAS ref | Distributed pbblas |
| `scalapack/` | **178** under `auxiliary/`, `eigenvalue/`, `factorization/`, `linear_solve/` | serial quad LAPACK on global matrix | Distributed |
| `mumps/` | **25** under `fortran/` and `c/` | `refmumps_quad` (sparse direct) | Both Fortran and C entrypoints |

Plus infrastructure:

- `blas/refblas/` — builds `refblas_quad.a` from `_refblas_src/`
- `common/` — `prec_kinds`, `compare`, `prec_report`, `test_data` modules used across domains
- `common/target_<name>/` — per-target conversion wrappers (cast `KIND=16` ↔ target precision)

### Testing model

Each test runs the same operation through two implementations and
compares the results in REAL(KIND=16) space:

- **Reference**: the vendored Netlib source compiled at quad precision
  via gfortran's `-freal-8-real-16` flag, producing a `refblas_quad`
  static library. Routine names stay as the canonical D/Z prefixes
  (`dgemm`, `ddot`, `zaxpy`, …) but operate at REAL(KIND=16).
- **Target under test**: the migrated extended-precision library
  (`qxblas`/`eyblas`/`mwblas`) for the active build target.

Inputs are generated at REAL(KIND=16). The reference call is direct;
the target call goes through a per-target wrapper that converts the
inputs to the target's precision (REAL(KIND=10) for kind10,
TYPE(real64x2) for multifloats) and converts the result back to quad
for comparison. Comparison precision is **always REAL(KIND=16)**,
regardless of which target is being tested.

Each test program writes one JSON file per `(routine, target)` pair to
`<build>/precision_reports/`. A Python aggregator
(`scripts/precision_report.py`) collates all JSON into a Markdown
summary table.

### Differential testing principle

The reference is *not* an oracle. It is the same Netlib BLAS source
the migrator processed, just compiled at quad precision instead of
being type-promoted by the migrator. **A failure that appears on both
sides is not a test failure** — `max_rel_err = |ref − got| / |ref|`
naturally returns ~0 when the two implementations diverge from
"truth" in the same way, which is what we want. The test fires only
when the migrated routine *diverges from* the reference algorithm.

For routines whose specification leaves some output slots unset
(e.g. `drotmg`'s `param(2)`/`param(5)` are only meaningful when the
flag is `1`), the test skips those slots — comparing them is a false
positive, not a real divergence.

### Directory layout (BLAS example)

```
test/integration/            (staged into the build tree as tests/)
├── README.md                 ← you are here
├── CMakeLists.txt            ← top-level dispatcher (added to staging build)
└── blas/
    ├── README.md
    ├── CMakeLists.txt        ← BLAS test build (selects target wrapper, links libs)
    ├── refblas/
    │   └── CMakeLists.txt    ← builds refblas_quad from _refblas_src/
    ├── common/               ← target-agnostic helpers (used by every test)
    │   ├── prec_kinds.f90    ← parameter ep = 16
    │   ├── compare.f90       ← rel_err helpers (scalar / vector / matrix, real / complex)
    │   ├── prec_report.f90   ← report_init / report_case / report_finalize (JSON writer)
    │   ├── test_data.f90     ← seeded RNG, gen_{vector,matrix}_{quad,complex}
    │   └── ref_quad_blas.f90 ← explicit interfaces to refblas_quad (D/Z names @ KIND=16)
    ├── target_kind10/
    │   └── target_blas.f90   ← wraps e/y-prefix routines, casts ep ↔ kind=10
    ├── target_kind16/
    │   └── target_blas.f90   ← wraps q/x-prefix routines (passthrough — kind=16 already)
    ├── target_multifloats/
    │   └── target_blas.f90   ← wraps dd/zz-prefix routines, splits ep ↔ TYPE(real64x2)
    ├── level1/               ← 17 BLAS Level 1 tests (dasum, daxpy, …, zscal, dzasum)
    ├── level2/               ← 13 BLAS Level 2 tests (dgemv, dgbmv, …, zhemv, zgerc)
    └── level3/               ← 10 BLAS Level 3 tests (dgemm, dsymm, …, ztrsm)
```

The tests are wired into the migrator's build pipeline through two
hooks:

1. `cmake/CMakeLists.txt` (the unified build) ends with
   `add_subdirectory(tests)` if `tests/` is present in the staging
   directory.
2. `codegen/migrator/__main__.py:cmd_stage` copies this directory
   (as `tests/`) and
   `extern/lapack-3.12.1/BLAS/SRC/` (renamed to `_refblas_src/`)
   into the staging directory, making it self-contained.

So the only command needed to test a target is:

```bash
python -m migrator stage /tmp/staging --target kind16 --libraries blas
cmake -S /tmp/staging -B /tmp/staging/build -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/staging/build -j8
ctest --test-dir /tmp/staging/build
python scripts/precision_report.py /tmp/staging/build/precision_reports/
```

### Per-target wrapper contract

Every `target_<target>/target_blas.f90` exposes the same module
`target_blas` with the same set of public symbols:

```fortran
module target_blas
    public :: target_name, target_eps      ! identifying constants
    public :: target_ddot, target_dgemm, … ! one wrapper per routine
end module
```

- `target_name` is a string used in JSON filenames
  (`<routine>.<target_name>.json`).
- `target_eps` is the target's machine epsilon, used by tests to size
  their tolerance.
- Each `target_<routine>` takes/returns REAL(KIND=16) (or
  COMPLEX(KIND=16)) and internally converts to/from the target's
  native type before/after calling the migrated routine.

The wrapper is the only place that knows the target's prefix
(`q`/`e`/`dd`) or precision. Test programs call `target_<routine>` and
stay precision-agnostic — the same test source runs unchanged against
all three targets.

### Adding a new test

For a routine `xroutine` that's already declared in
`common/ref_quad_blas.f90` and wrapped in every `target_*/target_blas.f90`:

1. Drop a new program into the appropriate `level{1,2,3}/`
   directory:

   ```fortran
   program test_xroutine
       use prec_kinds,    only: ep
       use prec_report,   only: report_init, report_case, report_finalize
       use compare,       only: max_rel_err_vec      ! or whichever shape fits
       use test_data,     only: gen_vector_quad      ! and any matrix helpers
       use target_blas,   only: target_name, target_eps, target_xroutine
       use ref_quad_blas, only: xroutine
       implicit none
       ! … generate inputs at REAL(KIND=16), call ref + target, compare,
       !   call report_case() per problem shape, report_finalize() at end.
   end program
   ```

2. The CMakeLists.txt picks up `test_*.f90` via `file(GLOB CONFIGURE_DEPENDS …)`,
   so no manual registration is needed.

3. For a routine **not** yet declared, add an explicit interface to
   `common/ref_quad_blas.f90` and a wrapper to *every*
   `target_<target>/target_blas.f90`. The wrapper signature is
   uniformly REAL(KIND=16) / COMPLEX(KIND=16) on the outside.

### Adding a new library subdirectory

The same scaffold extends to LAPACK, MUMPS, PBLAS:

1. Create a `<library>/` subdirectory here with the same internal
   layout as `blas/`: `CMakeLists.txt`, `common/`, `target_*/`, level
   subdirectories appropriate to the library.
2. The `CMakeLists.txt` at this directory's root already iterates over
   `${STAGED_LIBRARIES}` and adds any subdirectory whose
   `CMakeLists.txt` exists — no edits needed there.
3. The reference for LAPACK/etc. follows the same recipe as
   `blas/refblas/` — pull the vendored sources from
   `extern/lapack-3.12.1/SRC/` into the staging dir (extend
   `cmd_stage` similarly), and compile with `-freal-8-real-16`.

### Reference selection

`blas/CMakeLists.txt` exposes one cache option:

```
-DBLAS_REF=netlib   # default; build refblas_quad from vendored sources
-DBLAS_REF=system   # link the system -lblas (capped at ~1e-16, double precision)
```

`netlib` is the default because (a) it provides true quad-precision
reference values (so kind16 tests show bit-exact agreement and
multifloats tests can resolve their full ~32 digits), and (b) it
removes the dependency on whatever BLAS happens to be installed.
`system` is a fallback for environments where gfortran's
`-freal-8-real-16` is unavailable.

### Reading the report

Per-routine summary, sample (kind10 / kind16 / multifloats):

| routine | kind10 | kind16 | multifloats |
|---------|--------|--------|-------------|
| ddot    | 17.59  | exact  | 30.68       |
| dgemm   | 18.50  | exact  | 31.22       |
| zgemm   | 18.63  | exact  | 31.33       |

- "exact" → bit-identical (`max_rel_err = 0`). For kind16 this is the
  expected outcome on every routine: reference and target are both
  REAL(KIND=16) computing the same algorithm.
- A number → the worst-case digits of agreement
  (`-log10(max_rel_err)`) across that routine's problem shapes.
- A `✗` next to the number → at least one case exceeded its
  tolerance. The per-case detail section identifies which one.

Theoretical ceilings per target:

- kind10 (80-bit x87 extended): ~19 decimal digits
- kind16 (IEEE binary128): ~33 decimal digits — but capped at "exact"
  here because reference and target are computed identically
- multifloats (double-double): ~32 decimal digits

---

## How to run

| Test type | Command |
|---|---|
| Differential precision (BLAS) | `ctest -R '^blas_'` — writes JSON to `precision_reports/` |
| All differential-precision tests | `ctest` from the staging build dir |
| Aggregate precision JSON → Markdown | `uv run scripts/precision_report.py` |

For the overlay-vs-migrated test suite (epblas-parallel, epblas-openblas),
see [epblas-parallel/tests/](../../epblas-parallel/tests/).

## Quick mental model

- **precision tests** = migrated vs Netlib@KIND=16, all 9 library
  domains, JSON-aggregated via `precision_report.py`. Lives here.
- **overlay tests** = hand-written kernels vs migrated baseline at
  the same precision. Lives in `../../epblas-parallel/tests/`.

## Related docs

- `../../epblas-parallel/doc/optimization-findings.md` — overlay perf findings
- `../../scripts/precision_report.py` — aggregates per-test JSON reports across domains
