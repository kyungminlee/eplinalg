# Convergence

## What it is

Upstream BLAS / LAPACK / ScaLAPACK / MUMPS ship every routine in two
co-family halves: `S`/`D` for real, `C`/`Z` for complex. After type
migration, both halves of a co-family are supposed to map to the same
target file — `sgemm.f` and `dgemm.f` both become `qgemm.f` for
`kind16`. **Convergence** is the property that the two migrated bodies
are textually identical (up to whitespace, comments, and a handful of
documented cosmetic equivalences). When they aren't, that is a
**divergence** — a sign that the migrator missed something, or that
upstream's two halves genuinely differ in a way the migrator cannot
reconcile.

The migrator emits the D/Z half as canonical and discards the S/C half
silently when they converge. Divergences are surfaced by the `diverge`
subcommand (an earlier separate `converge` on-disk check was merged
into it).

## The report

```bash
uv run python -m migrator diverge codegen/recipes/lapack.yaml --target kind16
```

`diverge` re-migrates every co-family pair in memory and compares the
halves with `_canonicalize_for_compare` — the heavy normalizer
(uppercase folding, `KIND=` suffix stripping, `REAL`/`CMPLX` cast
unwrapping, prefix-token folding). It needs only the recipe; no
on-disk output directory is involved.

`migrator run` chains migrate → diverge → verify → build, so
convergence is checked on every full-pipeline invocation and exits
non-zero when any pair diverges. CI gates on the exit code.

The report itself is returned as a structured list by
`pipeline.run_divergence_report` — see `codegen/migrator/pipeline.py`
— so any custom rollup (per-target matrix, per-library counts,
machine-readable export) can be assembled from that return value
without re-running the migration.

## What gets suppressed

`diverge` deliberately tolerates a small set of cosmetic
upstream-drift patterns so the report focuses on real bugs:

- **Comments** — stripped before comparison; D/Z and S/C halves
  routinely paraphrase the same prose differently.
- **Semantic equivalences** — the `_canonicalize_for_compare` suite:
  `KIND=` suffix stripping, `REAL`/`CMPLX` cast unwrapping,
  declaration-list sorting, and friends. Each rule is documented in
  the source.
- **Pure precision-prefix local-variable drift** — if a hunk differs
  only at S↔D / C↔Z positions inside identifier tokens,
  `pipeline._filter_precision_drift` silently drops it (these are
  unclassified locals, not symbol-table entries — any
  properly-classified symbol has the same migrated name in both
  halves).

## When divergences are expected

Some pairs are known to diverge for upstream reasons (e.g. a routine
where the C-half intentionally calls a different helper than the
Z-half). Recipes whitelist these via the `expected_divergences:` list
(per-file) or `defer_all_divergences:` (whole recipe) — see
`config.py`. Pass `--no-whitelist` to `diverge` to see the report
with the whitelist bypassed. Alternatively, a recipe can pin one half
via `prefer_source:` so the drifting sibling never participates.

## Related

- [`architecture.md`](architecture.md) — overall pipeline.
- `pipeline._canonicalize_for_compare` — the normalization function,
  with extensive in-source comments on every rule and why it exists.
- `pipeline._filter_precision_drift` — the S↔D / C↔Z prefix-folding
  step applied by `diverge`.
