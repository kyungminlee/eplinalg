# tests/mumps

Differential precision tests for the migrated MUMPS sparse direct solver.
Mirrors the tests/blas / tests/lapack pattern: every check happens at
`REAL(KIND=16)`, the migrated `${LIB_PAIR_PREFIX}mumps` archive
(`qxmumps` / `eymumps` / `mwmumps` for kind16 / kind10 / multifloats) is
exercised both from Fortran (`call qmumps(id)`) and from C
(`qmumps_c(&id)` via the bridge described in
[`c/include/qmumps_c.h`](c/include/qmumps_c.h)), and per-case JSON
reports are written under `<build>/precision_reports/` for the standard
aggregator. All three targets pass 26/26 mumps ctests under
`linux-impi` (see the B4 entry in `CHANGELOG.md`).

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

If you restrict the staging via `--libraries`, the mumps tests still
need `scalapack_c` in the list — it supplies the precision-promoted C
clones `*lamov_` / `p*gemr2d_`, which are now folded into the
`scalapack` archive (no separate `scalapack_c` library is built, but
its sources must still be staged for the fold). Concretely:

```bash
uv run python -m migrator stage /tmp/stg-q --target kind16 \
    --libraries blas blacs ptzblas pbblas pblas \
                scalapack scalapack_c lapack mumps
```

(The default — no `--libraries` flag — stages everything and is fine.)

The `codegen/recipes/` and `cmake/` trees live in this single repo
(`eplinalg`) — the historical fm-mumps split was retired when
the mumps work merged into `tests` (see the B7 entry in `CHANGELOG.md`).

Each test is built twice: linked against real MPI (wrapped via
`mpiexec -n 1`, since the migrated qxmumps archive calls MPI primitives
unconditionally) and linked against the in-tree `mpiseq` archive (plain
binary, suffix `_seq`, no mpiexec needed). Both variants pass 26/26 on
all three targets with bit-identical JSON precision reports — see the
libmpiseq entries in `CHANGELOG.md` for the linkage mechanics.

The suite also runs at `np ≥ 2` (default `MUMPS_TEST_NPROCS=2`). Several
MUMPS API conventions only bite once a slave rank exists — host-only sparse
RHS, the distributed-solution `INFO(23)` slice size, and `-j1` serialization
— documented in
[`../../doc/user/api/mumps.md`](../../doc/user/api/mumps.md).

## Layout

```
tests/mumps/
├── CMakeLists.txt        — gates + bridge build + test registration
├── README.md             — this file
├── TODO.md               — open issues (currently the D1 input-validation watch)
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
`runtime/mpiseq/mpiseq_c_stubs.c` (folded into the `mpiseq` target when
`USE_LIBMPISEQ=ON` — see the `linux-libmpiseq` preset).

## Coverage

The 26 tests cover the `JOB=-1 → JOB=6 → JOB=-2` roundtrip, SYM
variants, ICNTL ordering choices, JOB phasing, multiple NRHS, error
paths (see TODO.md D1), and Fortran/C parity — each in real (`d`-named
drivers) and complex (`z`-named) form against the migrated archive.
Resolved design history lives in [CHANGELOG.md](CHANGELOG.md); the only
open item is TODO.md.
