# Test

Four layers, from cheapest to most expensive.

## 1. Migrator unit tests

```bash
uv run pytest        # repo root
```

Covers the transform passes, C migrator, staging, and CLI plumbing.
Run after any migrator change.

## 2. Codegen-level checks

- **Patch symmetry** (CI check): `uv run python -m migrator
  verify-patches <recipe>` — every patch touching a precision-prefixed
  file must touch all four siblings unless declared asymmetric.
- **Convergence**: `uv run python -m migrator converge <recipe> <out>
  --target <t>` — the authoritative dual-origin check: both precision
  halves of every routine must migrate to identical output. Methodology
  and current status: [convergence.md](convergence.md); command details:
  [codegen.md](codegen.md).

## 3. Stack tests (ctest)

Configure the staging tree with `MIGRATOR_BUILD_TESTS=ON` (the
default), then:

```bash
ctest --test-dir /tmp/stage-q/build -j2 --timeout 300
```

Keep `-j` low: oversubscribed MPI tests can deadlock indefinitely. The
`*_seq` test executables and BLACS/PBLAS p2p tests are only verified
against Intel MPI.

The suites are *differential*: migrated output is compared against
Netlib reference results or overlay implementations at the target
precision. Map of every test family: [../../test/integration/README.md](../../test/integration/README.md).

## 4. Release validation coverage

What a release has actually been checked against:

- **CI** (`.github/workflows/release.yml`): per-target staging with a
  convergence report; full build matrix 3 targets ×
  {gfortran-12, gfortran-15} × {openmpi, mpich, intelmpi, seq}; `ctest`
  in one slot (kind10 + gfortran-15 + intelmpi).
- **Off-CI**, per major change: a 774-configuration MUMPS sweep — all
  six extended arithmetics across orderings/ICNTL options/np≤4, with an
  MPFR backward-error check against each family's epsilon — run under
  MKL LP64 + Intel MPI in the hostile MKL-first link order, in both
  static-archive and repackaged-`.so` form.
