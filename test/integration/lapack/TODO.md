# tests/lapack — TODO

Resolved items have moved to `CHANGELOG.md`. Remaining work: parked
items with explicit reopen conditions.

## sb2st kernels (2 routines) — PARKED

`dsb2st_kernels`, `zhb2st_kernels` — SBR (Successive Band Reduction)
inner kernels, normally driven by `*sb2st` / `*hb2st`. Testing them in
isolation requires constructing a partial band-reduction state that
matches the kernel's mid-stream invariants — non-trivial.

These are the only two user-facing routines without a test driver per
the audit. The orchestrator paths (`*sb2st` / `*hb2st`) ARE tested
end-to-end, and they exercise the kernels on every call, so coverage
is effectively transitive.

**Reopen condition**: a regression or precision divergence
attributable to the kernels that the orchestrator tests fail to catch.
Until then, standalone drivers aren't worth the state-construction
work.

## kind4 / kind8 baseline missing the orbdb3 override — RESOLVED 2026-05-09

Closed by the pipeline refactor (`doc/archive/refactor-20260509.md`,
Phase E). `migrator stage --target kind{4,8}` now stages baseline
sources from `build/staged-sources/<lib>/` (which carries the
patched bodies) instead of `extern/<vendor>/.../SRC/`. The L01
LDX21 fix (and the rest of the patch catalogue) rides into
`_reflapack_src/` for every baseline build.

Verification: `grep -n "CALL DROT" /tmp/stage-baseline-k4/_reflapack_src/dorbdb3.f`
shows `LDX21` at the INCY arg; upstream still has `LDX11`. The
former parked failures (`lapack_test_dorbdb3`, `lapack_test_zunbdb3`)
should clear on the next baseline ctest run.
