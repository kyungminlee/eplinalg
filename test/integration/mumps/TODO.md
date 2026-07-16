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
