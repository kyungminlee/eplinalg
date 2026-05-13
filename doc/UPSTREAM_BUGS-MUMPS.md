# Upstream bugs: MUMPS 5.8.2

*Last catalogued: 2026-05-11. See `UPSTREAM_BUGS.md` for the cross-library
index, audit methodology, bug summary table, and how fixes are carried.*

This file collects MUMPS 5.8.2 bugs in the vendored
`external/MUMPS_5.8.2/` source. LAPACK and ScaLAPACK bugs are
catalogued in `UPSTREAM_BUGS-LAPACK.md` and
`UPSTREAM_BUGS-ScaLAPACK.md` respectively.

## MUMPS 5.8.2: bad-input cases SIGSEGV instead of returning `INFOG(1) < 0`

**Symptom.** The MUMPS user guide documents that input-validation
failures populate `id%INFOG(1)` with a negative error code (e.g.
`-16` for an invalid `N`, `-6` for IRN/JCN out of range) and return
control to the caller. In practice, MUMPS 5.8.2 SIGSEGVs on every
invalid-input case the differential-precision suite tried, before
any diagnostic is written. Concretely:

- `id%N = -1` (or `id%N = 0`) at JOB=1 — SIGSEGV inside analysis,
  no `INFOG(1)` set.
- `id%NNZ = -1` (negative count) — SIGSEGV during pattern setup.
- `id%IRN(1) = N + 5` (out-of-range row index, high) — SIGSEGV
  inside the input-reformatting path.
- `id%IRN(1) = 0` (out-of-range row index, low — zero is invalid
  under 1-based indexing) — same crash signature.
- `id%JCN(1) = N + 3` / `id%JCN(1) = 0` — symmetric to IRN.
- `id%NNZ` mismatched against the actual `IRN/JCN/A` sizes
  (e.g. `NNZ = 0` with non-empty arrays, or `NNZ` larger than
  what the arrays hold) — SIGSEGV during reformatting.

The crashes are unconditional: the analysis / pattern-setup code
doesn't bounds-check the user's IRN/JCN/N/NNZ before indexing into
internal work arrays sized from those same values, so a wrong value
either dereferences past the end of a fresh allocation or hits an
arena boundary.

**Root cause (inferred — not from upstream).** `*MUMPS_ANA_DRIVER`
and the IRN/JCN reformatters trust the assembled-format inputs the
user sets in the `id` struct. There is a `MUMPS_ANA_F_CHECK` (and
similar) helper internally, but it isn't invoked early enough — by
the time it would run, the analysis code has already indexed past
end-of-buffer using N or IRN values it never validated. So the
manual's promised `INFOG(1) < 0` return is a documented-but-unwired
contract: the validation pass exists in name but doesn't gate the
indexing it's supposed to protect.

**Affected files.** Symptom is observed at the `?MUMPS` driver
entry point; the actual crashing site varies per case (analysis
driver, IRN/JCN reformatter, pattern-set work-array sizing). All
under `external/MUMPS_5.8.2/src/`. Not narrowed to a specific
file because the workaround sidesteps the issue at the caller
rather than patching MUMPS internals.

**Workaround (in-tree, not an override).** Per
`tests/mumps/TODO.md` D1, the differential-precision wrapper at
`tests/mumps/common/target_mumps_body.fypp` exports
`check_dmumps_input` / `check_zmumps_input` that mirror the
documented validation contract before any `?MUMPS` call. They
return integer codes (`MIC_BAD_N`, `MIC_BAD_NNZ`, `MIC_BAD_IRN`,
`MIC_BAD_JCN`, `MIC_SIZE_MISMATCH`, `MIC_OK`). Tests at
`tests/mumps/fortran/test_{d,z}mumps_errors.f90` exercise each
class plus a final valid-input pass that reaches `MIC_OK` and
factors via `JOB=6` to confirm the wrapper isn't over-rejecting.
This is a *test-side* contract simulator — it lets us assert the
documented behavior without crashing the test harness, but does
not modify MUMPS itself. A real caller still has to pre-validate
their inputs externally to avoid the SIGSEGV.

**Why upstream's tests miss it.** MUMPS's own test drivers
(`*examples/`) feed valid matrices only. The validation contract
documented in the user guide is asserted by the prose, not by any
shipped test that passes a deliberately bad `id` struct and
checks `INFOG(1)`. So the gap between "documented" and "wired" is
invisible to the upstream test cycle.

**Upstream report.** Not yet filed.
