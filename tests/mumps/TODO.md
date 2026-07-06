# tests/mumps — TODO

Resolved items have moved to `CHANGELOG.md`. The original 2026-04-29
"deferred test implementation plan" is closed: tests/CMakeLists.txt
ships, all 26 mumps ctests pass on every target under
`cmake --preset=linux-impi`, and `README.md` documents the run flow.

## D1 — MUMPS 5.8.2 doesn't validate bad inputs (sidestepped, monitoring)

Upstream MUMPS 5.8.2 crashes (SIGSEGV inside analysis or factorization)
on inputs the manual claims are rejected via `INFOG(1) < 0`:

- `id%N = -1` → SIGSEGV before any diagnostic.
- `id%IRN(1) = N + 5` or `id%JCN(1) = N + 3` → SIGSEGV inside
  `qmumps_validate_input_` or callees.
- Mismatched `NNZ` (e.g. `NNZ = 0` with non-empty IRN/JCN/A) → SIGSEGV
  during reformatting.

`target_mumps` exports `check_dmumps_input` / `check_zmumps_input`
plus `MIC_*` codes (see `tests/mumps/common/target_mumps_body.fypp`).
The wrapper layers above MUMPS without touching `id%INFOG(1)`. Tests
at `tests/mumps/fortran/test_{d,z}mumps_errors.f90` exercise BAD_N
(neg + zero), BAD_NNZ, BAD_IRN/JCN (high + zero), SIZE_MISMATCH, plus
a final valid-input pass that reaches MIC_OK and factors via JOB=6.

**Reopen condition**: when MUMPS upstream tightens validation and
returns the documented `INFOG(1) = -16` / `-6` codes, drop the wrapper
and let the tests check `INFOG(1)` directly.

## libmpiseq linkage — landed on all three targets

Each test source under `tests/mumps/{fortran,c}/` is built twice — once
linked against real MPI (`mpiexec -n 1`), once linked against the
in-tree `mpiseq` archive (plain binary, suffix `_seq`). libmpiseq is
Intel-MPI ABI compatible: runtime/mpiseq/CMakeLists.txt overlays the configured
MPI's `mpi.h` / `mpif.h` onto `_mpiseq_src/`, libseq's bundled `mpic.c`
is replaced by `runtime/mpiseq/mpiseq_c_stubs.c` (compiled against Intel's
signatures), and `migrator stage` patches libseq's `mpi.f` to extend
`MUMPS_COPY` with the precision-specific datatype cases:

- `MPI_REAL16` / `MPI_COMPLEX32` (kind16)
- `MPI_LONG_DOUBLE` / `MPI_C_LONG_DOUBLE_COMPLEX` (kind10 — 80-bit
  extended types map to MPI's long-double tokens, no `MPI_REAL10`
  exists in standard MPI)
- `268435472` / `268435488` (multifloats — sentinel handles for
  `MPI_FLOAT64X2` / `MPI_COMPLEX64X2`, see below)

**Multifloats sentinel scheme**: `multifloats_mpi.cpp` registers its
custom datatypes via `MPI_Type_contiguous(count, MPI_DOUBLE, &out)`. The
libmpiseq C stub for `MPI_Type_contiguous` encodes the total byte size
in low bits with a high tag (`0x10000000 | nbytes`), so float64x2 → 16
bytes → handle `0x10000010`, complex64x2 → 32 bytes → handle
`0x10000020`. `MPI_Type_c2f` is a passthrough cast in Intel mpi.h, so
the Fortran handle is the same value, and `MUMPS_COPY` dispatches on
those constants to a generic byte-block copy helper. `MPI_Op_create`
returns distinct non-null sentinels above `MPI_OP_NULL` (0x18000000);
the user-op callback never fires under single-rank ALLREDUCE. The
0x10000000 / 0x18000000 ranges are well clear of Intel's 0x4c00****
datatype constants.

**Status (all three targets)**: 26/26 `_seq` ctests pass; per-test JSON
precision reports are bit-identical to the impi-linked runs (verified
via `diff -q`).
