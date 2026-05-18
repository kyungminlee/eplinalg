# fortran-migrator

Automated type-migration pipeline for the classical numerical-linear-algebra
stack — BLAS, BLACS, LAPACK, PBLAS, ScaLAPACK, MUMPS — retargeted from the
standard `DOUBLE PRECISION` / `DOUBLE COMPLEX` working type to a wider
arithmetic of your choice.

Three targets ship out of the box:

| target        | real type              | complex type           | real prefix | complex prefix |
|---------------|------------------------|------------------------|-------------|----------------|
| `kind10`      | `REAL(KIND=10)`        | `COMPLEX(KIND=10)`     | `e`         | `y`            |
| `kind16`      | `REAL(KIND=16)`        | `COMPLEX(KIND=16)`     | `q`         | `x`            |
| `multifloats` | `TYPE(float64x2)`      | `TYPE(complex64x2)`    | `m`         | `w`            |

`kind10` / `kind16` rely on the compiler's native extended-precision modes
(`__float80` / `__float128` via gfortran). `multifloats` uses the
double-double library at <https://github.com/kyungminlee/multifloats>,
fetched via CMake `FetchContent`.

All three retain the co-family structure of the original source: where
upstream has `dgemm.f` and `zgemm.f`, the migrator emits a single
`qgemm.f` (or `ddgemm.f` / `egemm.f`) whose body is identical regardless
of which half of the family it came from. Divergences — files where
`s/d` and `c/z` halves disagree after migration — are reported.

## Quick start

```bash
# Install deps (Python 3.11+, gfortran ≥ 11, MPI, CMake ≥ 3.20)
uv sync

# Migrate and build the full stack for a target into a staging dir
uv run python -m migrator stage /tmp/stage-q --target kind16
cmake -S /tmp/stage-q -B /tmp/stage-q/build -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/stage-q/build -j8
```

### MPI: use the `linux-impi` preset

The canonical MPI for builds and tests is Intel oneAPI MPI. `cmake/CMakePresets.json`
ships two presets:

| preset             | when to use                                                          |
|--------------------|----------------------------------------------------------------------|
| `linux-impi`       | default — points CMake at `/opt/intel/oneapi/mpi/latest` wrappers    |
| `linux-system-mpi` | CI / machines without Intel MPI; uses whatever `find_package(MPI)` finds (OpenMPI, MPICH) |

The staging tree includes the preset file, so:

```bash
cmake -S /tmp/stage-q --preset linux-impi
cmake --build /tmp/stage-q/build -j8
ctest  --preset linux-impi --test-dir /tmp/stage-q/build
```

`linux-system-mpi` works for compiling the libraries themselves, but the
`*_seq` test executables (which link `libmpiseq` while including a system
`mpi.h`) and BLACS/PBLAS p2p tests are only verified against Intel MPI;
OpenMPI/MPICH may hang or fail to link.

Resulting archives: `libqblas-<tag>.a`, `libqlapack-<tag>.a`,
`libqscalapack-<tag>.a`, etc.

Link against them from a Fortran program:

```fortran
program demo
  implicit none
  real(16) :: x(3), r
  real(16), external :: qasum
  x = [1.0_16, 2.0_16, 3.0_16]
  r = qasum(3, x, 1)
  print *, r
end program
```

```bash
gfortran demo.f90 -o demo \
    -Wl,--start-group /tmp/stage-q/build/lib*.a -Wl,--end-group \
    -lmpi -lmpicxx -lmpifort -lquadmath -lstdc++ -lpthread
```

## How it works

The migrator parses each source file with a compiler-based front-end
(Flang or GFortran) to extract structural facts — symbol kinds, call
graph, EXTERNAL/INTRINSIC tables — then applies regex-based rewrites
guided by those facts. The hybrid keeps the transform syntactically
aware while preserving formatting, comments, and preprocessor
directives. C libraries in the stack (BLACS, PBLAS) are handled by
template-based cloning with mechanical type substitution rather than
full parsing.

See [`doc/README.md`](doc/README.md) for the full documentation
index — including the developer guide, architecture overview, and
upstream-bug catalogues.

## CLI

```
migrator migrate   <recipe> <out>   # migrate source files only
migrator build     <recipe> <out>   # generate CMake + build one library
migrator run       <recipe> <work>  # full: migrate + build + verify
migrator stage     <dir>            # migrate all libraries into a unified CMake tree
migrator diverge   <recipe> <out>   # report s/d vs c/z divergences
migrator converge  <recipe> <out>   # whitespace-tolerant convergence check
migrator verify    <recipe> <out>   # post-migration verification
```

All commands accept `--target {kind10,kind16,multifloats}` or a path to a
target `.yaml`. `stage` is the usual entry point — it produces a
self-contained directory that builds with plain CMake.

## Repository layout

```
recipes/           # per-library YAML recipes (blas.yaml, lapack.yaml, …)
targets/           # per-target YAML configs (kind10, kind16, multifloats)
src/migrator/      # the migrator
external/          # vendored upstream sources (LAPACK, MUMPS, ScaLAPACK,
                   # Intel MPI headers, multifloats MPI-bridge companion)
cmake/             # staging build system
```

### recipes vs targets

A **recipe** describes a library: where its source lives, which files to
skip, which modules to copy verbatim, which routine families exist.
A **target** describes a retarget: the type system, prefix conventions,
and compile-time overlays (extra Fortran module helpers, C++ bridge
header, MPI datatype names). The migrator is the cartesian product of
the two.

Recipes live in `recipes/*.yaml` with per-library sidecar directories
(`recipes/<lib>/`) for line-level manifests and hand-written override
modules too library-specific to put in YAML. See `recipes/README.md`
for the sidecar conventions.

### External dependencies

| item                   | source                                      |
|------------------------|---------------------------------------------|
| LAPACK 3.12.1          | vendored under `external/lapack-3.12.1/`    |
| MUMPS 5.8.2            | vendored under `external/MUMPS_5.8.2/`      |
| ScaLAPACK 2.2.3        | vendored under `external/scalapack-2.2.3/`  |
| Intel MPI headers      | vendored under `external/impi-headers/` (compile-time only — link/run against Intel oneAPI MPI via the `linux-impi` preset; other MPIs are best-effort) |
| multifloats            | fetched at CMake time from GitHub (`FetchContent`) |
| `multifloats-mpi`      | `external/multifloats-mpi/` — MPI bridge (datatype + reduction ops) |

To pin a specific multifloats release, pass
`-DMULTIFLOATS_GIT_TAG=v0.2.3` to CMake. For an offline build, pass
`-DMULTIFLOATS_DIR=/path/to/local/multifloats`.

## Tests

```bash
uv run pytest
```

The fast unit tests cover the migrator's regex/AST transforms per target.
The end-to-end build-and-link path is exercised manually via `migrator
stage` + CMake (see **Quick start**).

## Status

- `kind10` (e/y), `kind16` (q/x), and `multifloats` (m/w) all build
  the full blas/xblas/blacs/lapack/pbblas/pblas/ptzblas/scalapack/scalapack_c
  archives. The differential precision suite runs 1 125 tests per
  target end-to-end and **all three targets pass 1 125 / 1 125** —
  the previously-persistent `pzdbtrsv` `MPI_Finalize` crash was
  closed by extending the upstream `PDDBTRS`/`PZDBTRS` `LWMIN`
  override to the `*trsv` oracle call site (commit `eec88b3`).
  See `tests/REPORT.md` for the per-target breakdown.
- MUMPS kind10 / kind16 build and pass all 26 MUMPS tests
  (Fortran drivers + C-bridge parity / basic / sym). Keep-kind
  manifest + EP bridge modules handle the DP-stable shared modules;
  see `recipes/README.md` for the mechanism.

## License

Vendored upstream libraries retain their original licenses (see the
`LICENSE*` files under `external/*/`). The migrator itself is
unlicensed for now.
