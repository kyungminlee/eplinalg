# Development Workflow

Configure → build → test → debug → release, for developers of eplinalg
itself. Migrator internals: [developer.md](developer.md). Migrator CLI
reference: [usage.md](usage.md). Consuming released binaries:
[binary-releases.md](binary-releases.md).

## Configure and build the full stack

```bash
# migrate every library into a self-contained CMake tree
cd src
uv run python -m migrator stage /tmp/stage-q --target kind16

# configure + build (the preset file is copied into the staging tree)
cmake -S /tmp/stage-q --preset linux-impi
cmake --build /tmp/stage-q/build -j8
```

Presets (`cmake/CMakePresets.json`):

| preset | use |
|---|---|
| `linux-impi` | canonical: Intel oneAPI MPI (`source /opt/intel/oneapi/setvars.sh` first) |
| `archlinux-impi-gfortran15` | Intel MPI with pinned gfortran-15 (Arch host) |
| `linux-system-mpi` | whatever `find_package(MPI)` finds; fallback when a preset's pinned compilers aren't installed |

The primary correctness target is **MKL LP64 + Intel MPI on Linux**;
OpenMPI/MPICH builds are convenience validation.

Caveats:

- `stage` snapshots `tests/CMakeLists.txt` and the preset file into the
  staging tree. After editing them in the repo, re-stage — rebuilding an
  old tree silently uses the stale copies.
- A full build tree per target is large; building all targets means
  three full linalg stacks. Remove each build tree before starting the
  next on small disks.

## Test

- Migrator unit tests: `uv run pytest` (repo root).
- Patch symmetry (CI check): `uv run python -m migrator verify-patches <recipe>`.
- Convergence: `uv run python -m migrator converge <recipe> <out> --target <t>`
  — the authoritative dual-origin check (see usage.md).
- Stack tests: configure with `-DMIGRATOR_BUILD_TESTS=ON`, then

  ```bash
  ctest --test-dir /tmp/stage-q/build -j2 --timeout 300
  ```

  Keep `-j` low: oversubscribed MPI tests can deadlock indefinitely.
  The `*_seq` test executables and BLACS/PBLAS p2p tests are only
  verified against Intel MPI.

## Release validation coverage

What a release has actually been checked against:

- **CI** (`.github/workflows/release.yml`): per-target staging with a
  convergence report; full build matrix 3 targets ×
  {gfortran-12, gfortran-15} × {openmpi, mpich, intelmpi, seq}; `ctest`
  in one slot (kind10 + gfortran-15 + intelmpi).
- **Off-CI**, per major change: a 774-configuration MUMPS sweep — all
  six extended arithmetics across orderings/ICNTL options/np≤4, with an
  MPFR backward-error check against each family's epsilon — run under
  MKL LP64 + Intel MPI in the hostile MKL-first link order, in both
  static-archive and repackaged-`.so` form (v0.13.1: 774/774 both).

## Known pitfalls

- **Intel MPI user-op reduce**: Intel MPI's shm reduce can hand a
  commutative user-op callback a buffer aligned to 8-but-not-16 bytes
  above the short-message cutoff; aligned-SSE callbacks (`__float128`)
  fault. The quad/multifloats reduce ops therefore register with
  `commute=0` under an Intel MPI compile-time guard (`commute=1`
  elsewhere — OpenMPI verified unaffected). Don't "optimize" it back
  without re-testing under Intel MPI.
- **MPI init contract**: extended collectives/MUMPS require
  `multifloats_mpi_init()` / `quad_mpi_init()` after `MPI_Init`
  (kind10 needs nothing). Missing init is the usual cause of
  "Intel MPI crashes in a quad reduction".
- **Shared-library links**: every produced `.so` needs
  `-Wl,--no-define-common` (Intel MPI COMMON sentinels fork per-`.so`
  otherwise → PT-Scotch parallel analysis dies, `INFOG(1)=-9980`) and
  `-Wl,-z,now` (lazy PLT resolution corrupts live vector-register FP
  state → flaky wrong double-double limbs). Details and the one
  exception: [mkl-coexistence.md](mkl-coexistence.md).
- **Scotch vs `_FORTIFY_SOURCE`**: fortify false-aborts Scotch's based
  arrays ("buffer overflow detected" at `ICNTL(7)=3`); the vendored
  build compiles scotch/esmumps with `-D_FORTIFY_SOURCE=0`.
- **Never rename shipped symbols**: MKL coexistence is solved by the
  source-level `ep_` privatization of the family-independent engine,
  not by link tricks; if a symbol exists in MKL, MKL must provide it.

## Releasing

1. Merge `develop` → `main`; create an annotated tag `vX.Y.Z`; push the
   tag.
2. The `v*` tag triggers `.github/workflows/release.yml`: stage each
   target, emit convergence reports, build the matrix, run the test
   slot, and publish 25 assets (one combined archive + 24 per-combo)
   with generated release notes.
3. List versions with `git tag --sort=-v:refname` — plain `git tag`
   sorts lexically (`v0.13.0` before `v0.8.0`).
