# fortran-migrator documentation

## Guide

User- and contributor-facing reference for the migration tool itself.

- [architecture.md](guide/architecture.md) — system overview, components, data flow
- [usage.md](guide/usage.md) — running the `migrator` CLI
- [recipes.md](guide/recipes.md) — writing library recipes (YAML)
- [developer.md](guide/developer.md) — contributing, build, repo layout
- [intrinsics.md](guide/intrinsics.md) — Fortran intrinsics reference for KIND=16

## Migrator output

Catalogues describing what the migration pipeline emits.

- [procedures.md](output/procedures.md) — auto-generated routine cross-reference (Netlib → kind10/kind16/multifloats), produced by `tools/gen_procedures.py`
- [convergence.md](output/convergence.md) — what "convergence" means and how it's measured
- [kind16-divergences.md](output/kind16-divergences.md) — per-routine KIND=16 divergence analysis

## Upstream bugs

Bugs in the vendored Netlib sources that the migrator works around without editing `external/`.

- [upstream-bugs/](upstream-bugs/README.md) — index, methodology, cross-library summary
  - [lapack.md](upstream-bugs/lapack.md)
  - [scalapack.md](upstream-bugs/scalapack.md)
  - [mumps.md](upstream-bugs/mumps.md)
  - [lapack-nits.md](upstream-bugs/lapack-nits.md) — cosmetic dead declarations (catalogued, not patched)

## Parallel BLAS overlay

Branch-scoped sub-project: a hand-written, thread-parallel BLAS overlay for the extended-precision targets.

- [parallel-blas/design.md](parallel-blas/design.md) — overlay design and phased plan
- [parallel-blas/optimization-findings.md](parallel-blas/optimization-findings.md) — what was tried, what worked, what didn't
- [`../reports/overlay-coverage.md`](../reports/overlay-coverage.md) — auto-generated coverage + speedup table

## Test suites

- [`../tests/README.md`](../tests/README.md) — map of every test family (overlay-vs-migrated and migrated-vs-Netlib schemes)

## Archive

- [archive/](archive/) — historical surveys and timestamped reports
