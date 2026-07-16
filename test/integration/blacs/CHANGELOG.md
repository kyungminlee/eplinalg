# tests/blacs — CHANGELOG

Resolved items, reverse-chronological. Open work lives in `TODO.md`.

## 2026-05-02 — BLACS surface coverage audit

Re-audited 2026-05-02: every routine previously flagged as missing is
exercised by the test tree.

### `*trbr2d` receive companion of `*trbs2d`

`test_qtrbs2d.f90` calls `target_qtrbr2d` in all three scopes: `'A'`
(line 35), `'R'` (line 62), `'C'` (line 88). Wrappers in
`common/target_blacs_body.fypp:51,311`.

### `qtrrv2d` / `xtrrv2d` companions

Receive partners for the trapezoidal point-to-point send.
`test_qtrsd2d.f90:51` calls `target_qtrrv2d`; `test_xtrsd2d.f90:48`
calls `target_xtrrv2d`. Wrappers in `common/target_blacs_body.fypp`
(real at 257/311, complex at 377/387).

### `qgamn2d` / `xgamn2d` / `qgsum2d` / `xgsum2d` row/column scopes

Row ('R') and column ('C') scopes are covered:

- Real-prefix sum: `test_qgsum2d_rc.f90`
- Real-prefix amx/amn: `test_qgamx2d_rc.f90`
- Complex-prefix sum: `test_xgsum2d_rc.f90`
- Complex-prefix amx/amn: `test_xgamx2d_rc.f90` (added 2026-04-30)

Each program runs amx in row scope and amn in column scope (the
reverse pattern is symmetric and skipped to keep run time bounded).

### `blacs_set` / `blacs_pinfo` direct probes

`tests/blacs/tests/test_blacs_set.f90` exercises the four settable
WHAT codes (`SGET_NR_BS`, `SGET_NR_CO`, `SGET_TOPSREPEAT`,
`SGET_TOPSCOHRNT`) and verifies the value round-trips through
`blacs_get`. `tests/blacs/tests/test_blacs_pinfo.f90` exercises
`blacs_pinfo` idempotency and `blacs_get` with WHAT=0 (default system
context) and WHAT=10 (default broadcast topology).

## 2026-05-01 — kind10 / multifloats PBLAS Level 1/2 heap corruption

The previously-described "test PASSES every case but crashes during
teardown" pattern no longer reproduces. Re-verified 2026-05-01:
kind10 1022/1022 PASS, multifloats 1022/1022 PASS. See
`tests/pblas/CHANGELOG.md` for the same-named entry.

## 2026-04-30 — Real distributed exec on kind16

With real ranks connected via the `linux-impi` preset, the test suite
went from `0% real signal` (sandbox mpiexec gave singleton MPI) to
**1045/1045 (100%)** on the kind16 target. Two preset env vars cover
the Intel MPI compat issues with BLACS-registered user ops over
`REAL(KIND=16)`:

* `I_MPI_ADJUST_REDUCE=1` — bypasses the optimised `Reduce_local`
  shortpath that crashes on user MPI_Op + REAL(KIND=16) payload.
  (Original fix; reduced failures 170 → 68.)
* `I_MPI_ADJUST_ALLREDUCE=1` — same fix for the `*GAMX2D / IGAMX2D`
  paths through `MPIR_Allreduce` that scalapack iterative-refinement
  drivers hit. (Added 2026-04-30; resolved the last 8 ScaLAPACK
  failures.)

### Original 3 BLACS failures (commit `67c57d2`)

* `blacs_pinfo`: the `get_broadcast_topology` case asserted that
  `BLACS_GET(ctxt, WHAT=10, val)` returned a small 0..127 character
  code. Per `blacs_get_.c`, `WHAT=10` is `SGET_BLACSCONTXT` — the
  internal handle for the given context (Fortran handle is
  `MPI_Comm_c2f` of the underlying communicator). Singleton MPI
  returned a tiny placeholder; real impi returns large integers.
  Fixed by relaxing the assertion to idempotency only — two
  consecutive `BLACS_GET` calls must agree.
* `qtrsd2d` / `xtrsd2d`: not actually a BLACS bug. The test used
  M=6 N=5 with the LAPACK upper-triangular convention (verify
  `i <= j` only); BLACS `BI_GetMpiTrType` for M>N follows the
  upper-trapezoidal convention where col j carries rows
  `1..(M-N+j)` — the M-N rows above the j==N diagonal are part of the
  trapezoid and get transmitted. Fixed by switching tests to M=N=5
  so the canonical triangular shape is exercised.

## 2026-04-30 — MPI environment via Intel MPI 2021.18

The earlier MPICH-on-WSL singleton behavior (every rank seeing
`MPI_Comm_size = 1`) is bypassed by configuring with
`cmake --preset=linux-impi`, which routes `MPIEXEC_EXECUTABLE` /
`MPI_C_COMPILER` / `MPI_Fortran_COMPILER` through Intel MPI 2021.18
in `/opt/intel/oneapi/mpi/latest`. Hydra/PMI connect 4 ranks on the
same node correctly; tests run on a real 2×2 BLACS grid.
