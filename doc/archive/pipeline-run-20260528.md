# Pipeline run — full sweep, 2026-05-28

**Branch:** `develop` @ `8bbc65b6` (matches `origin/develop`)
**Targets exercised:** kind16, kind10, multifloats
**Pipeline:** `migrator stage` → `cmake -S … --preset linux-impi` →
`cmake --build` → `ctest`

## Result

| target       | stage | build  | ctest         | green? |
|--------------|-------|--------|---------------|--------|
| kind16       | 45s   | 6m16s  | 1152/1152     | ✓ (with fabric workaround) |
| kind10       | 37s   | 5m22s  | 1151/1152     | ✗ (1 pre-existing std-precision failure) |
| multifloats  | ~35s  | 9m04s  | 1152/1152     | ✓ (with fabric workaround + matched gcc-15) |

All three targets `migrator stage` and `cmake --build` cleanly on first
pass once the environmental issues below are resolved. No migrator
regressions surfaced.

## Findings

### F1 — Intel MPI fabric init fails on Arch (env, all targets)

Out-of-the-box on this box (Arch Linux + Intel oneAPI MPI 2021.14), every
MPI-using test (`blacs_test_*`, `pbblas_test_*`, `scalapack_test_*`,
`mumps_test_*`) aborts in `MPIDI_OFI_mpi_init_hook(1594)` with
`Unknown error class`. Affects ~306/1152 tests per target.

**Workaround:**
```
FI_PROVIDER=tcp I_MPI_FABRICS=shm ctest …
```
With those env vars set, the same builds pass cleanly.

**Followup:** Consider exporting the workaround as defaults in
`tests/CTestCustom.cmake` (or in the `linux-impi` preset's `environment`
block) for Arch developers. Other Linux distros with bare-metal Intel MPI
may need similar.

### F2 — multifloats LTO version mismatch when toolchain isn't pinned

`linux-impi` preset pins `CMAKE_Fortran_COMPILER=gfortran-15` (required
for Intel MPI .mod compatibility on Arch GCC-16 systems) but leaves
`CMAKE_C_COMPILER`/`CMAKE_CXX_COMPILER` unpinned. On Arch where
`cc`/`c++` are GCC 16, the multifloats sub-build emits a GCC-15.1
LTO bytecode in `libmultifloatsf.a` (Fortran side) while the C++ side
links via GCC-16's `lto-wrapper`, which rejects the older bytecode:

```
lto1: fatal error: bytecode stream in file '_mf/fsrc/libmultifloatsf.a'
generated with LTO version 15.1 instead of the expected 16.0
```

Affects only the `multifloats` target (kind10/kind16 don't pull
multifloats).

**Workaround:**
```
cmake -S /tmp/stage-m --preset linux-impi \
    -DCMAKE_C_COMPILER=gcc-15 -DCMAKE_CXX_COMPILER=g++-15
```

**Followup:** Either (a) pin `CMAKE_C_COMPILER=gcc-15`,
`CMAKE_CXX_COMPILER=g++-15` in the `linux-impi` preset so the whole
toolchain travels together, or (b) disable LTO when the resolved Fortran
and C++ compiler majors disagree (probably (a) — same root cause as the
gfortran pin).

### F3 — `lapack_test_zgeqp3rk` overflow under kind10 only

Test #408 passes under kind16 and multifloats but fails under kind10
with floating-point exceptions and `max_rel_err=********`:

```
zgeqp3rk [m=16,n=8]  max_rel_err=************  digits=****** FAIL
zgeqp3rk [m=24,n=12] max_rel_err=************  digits=****** FAIL
IEEE_INVALID_FLAG IEEE_OVERFLOW_FLAG IEEE_DENORMAL
```

The `z*` test runs against the unmigrated standard-precision reference
LAPACK, so this isn't a migration defect — the same source builds and
tests fine under kind16/multifloats but blows up under kind10. Likely
exposed by the kind10 build flags or by the `REAL(KIND=10)` x87 80-bit
ABI interacting badly with `zgeqp3rk`'s internal scaling. **Not a
blocker for kind10 production builds** (the migrated `y*` rank-revealing
path is the one users actually call), but worth a separate look.

**Followup:** Bisect compiler flags between kind10 and kind16 builds;
isolate the kind10-specific flag that destabilizes the std-precision z
QR-with-pivoting reference test.

## Cleanup

- `/tmp/stage-{q,e,m}` removed after each target's run.
- `external/openblas/` gitignored (commit `8bbc65b6`).
- Code-review handoff `BUGFIX-20260527.md` archived to
  `doc/archive/code-review-20260527.md` (commit `0942aef6`).
