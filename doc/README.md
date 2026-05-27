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

The hand-written extended-precision BLAS overlay moved to a separate
project: [epparablas](../../epparablas/). Migrator-side, the public
`${LIB_PREFIX}blas` target is now the plain serial migrated archive;
the overlay-equipped composite (`epparablas::eblas` and friends)
ships from the epparablas package.

## Test suites

- [`../tests/README.md`](../tests/README.md) — map of every test family (overlay-vs-migrated and migrated-vs-Netlib schemes)

## Archive

- [archive/](archive/) — historical surveys and timestamped reports
