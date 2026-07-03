# tests/mumps

Differential precision tests for the migrated MUMPS sparse direct solver.
Mirrors the tests/blas / tests/lapack pattern: every check happens at
`REAL(KIND=16)`, the migrated `${LIB_PREFIX}mumps` archive
(`qmumps` / `emumps` / `mmumps` for kind16 / kind10 / multifloats) is
exercised both from Fortran (`call qmumps(id)`) and from C
(`qmumps_c(&id)` via the bridge described in
[`c/include/qmumps_c.h`](c/include/qmumps_c.h)), and per-case JSON
reports are written under `<build>/precision_reports/` for the standard
aggregator. All three targets pass 26/26 mumps ctests under
`linux-impi` (B4 RESOLVED 2026-05-04).

## How to run

```bash
cd /home/kyungminlee/code/eplinalg/src
uv run python -m migrator stage /tmp/stg-q --target kind16 --parser gfortran
cmake -S /tmp/stg-q -B /tmp/stg-q/build --preset=linux-impi
cmake --build /tmp/stg-q/build -j8
ctest --test-dir /tmp/stg-q/build -R '^mumps_' --output-on-failure
```

Swap `--target kind16` for `kind10` or `multifloats` to run the same
suite on the other targets — no other flags change.

If you restrict the staging via `--libraries`, the mumps tests need
`scalapack_c` in the list (it ships the precision-promoted C clones
`*lamov_` / `p*gemr2d_` that the test executables link against).
Concretely:

```bash
uv run python -m migrator stage /tmp/stg-q --target kind16 \
    --libraries blas blacs ptzblas pbblas pblas \
                scalapack scalapack_c lapack mumps
```

(The default — no `--libraries` flag — stages everything and is fine.)

The `recipes/` and `cmake/` trees live in this single repo
(`eplinalg`) — the historical fm-mumps split was retired when
the mumps work merged into `tests` (see TODO.md B7).

Tests are wrapped via `mpiexec -n 1` because the migrated qmumps archive
calls MPI primitives unconditionally (see TODO.md B3 for the
libmpiseq-based no-mpiexec alternative — currently a known incomplete
path).

The suite also runs at `np ≥ 2` (default `MUMPS_TEST_NPROCS=2`). Several
MUMPS API conventions only bite once a slave rank exists — host-only sparse
RHS, the distributed-solution `INFO(23)` slice size, and `-j1` serialization
— documented in
[`../../doc/mumps-implicit-constraints.md`](../../doc/mumps-implicit-constraints.md).

## Layout

```
tests/mumps/
├── CMakeLists.txt        — gates + bridge build + test registration
├── README.md             — this file
├── TODO.md               — open issues, design plan, deferred coverage
├── common/               — shared helpers (prec_kinds, compare,
│                           prec_report, test_data_mumps,
│                           ref_quad_lapack_solve, target_mumps_body.fypp)
├── target_kind16/        — kind16 fypp shim setting prefix=q/x
├── target_kind10/        — kind10 fypp shim setting prefix=e/y
├── target_multifloats/   — multifloats fypp shim setting prefix=m/w
├── fortran/              — test_*mumps_*.f90 drivers
└── c/
    ├── include/          — quad-precision header overrides for the bridge
    └── test_*mumps_c_*.c — C drivers
```

Supplementary libmpiseq C-side stubs live alongside the Fortran ones at
`cmake/mpiseq_c_stubs.c` (folded into the `mpiseq` target when
`USE_LIBMPISEQ=ON` — see the `linux-libmpiseq` preset).

## Coverage plan

See [TODO.md](TODO.md) for the full design and the per-test parameter
matrix. The first landed tests cover the `JOB=-1 → JOB=6 → JOB=-2`
roundtrip on a small unsymmetric problem; subsequent tests will sweep
SYM, ICNTL ordering, JOB phasing, NRHS, error paths, and the C-side
parity equivalent.
