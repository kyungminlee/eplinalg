# Configure

CMake configures a *staging tree* produced by `migrator stage`, not the
repo itself — see [build.md](build.md) for the staging step. The CMake
system (`cmake/CMakeLists.txt`, presets, toolchains) is copied into the
staging tree by `stage`.

## Presets

`CMakePresets.json`:

| preset | use |
|---|---|
| `linux-impi` | canonical: Intel oneAPI MPI (`source /opt/intel/oneapi/setvars.sh` first) |
| `archlinux-impi-gfortran15` | Intel MPI with pinned gfortran-15 (Arch host) |
| `linux-system-mpi` | whatever `find_package(MPI)` finds; fallback when a preset's pinned compilers aren't installed |

```bash
cmake -S /tmp/stage-q --preset linux-impi
```

The primary correctness target is **MKL LP64 + Intel MPI on Linux**;
OpenMPI/MPICH builds are convenience validation.

## Options

| option | default | effect |
|---|---|---|
| `MIGRATOR_BUILD_TESTS` | `ON` | build the differential-precision test suite (`tests/`, adds ctest targets) |
| `MUMPS_LIBSEQ_RELEASE` | `OFF` | build the sequential libmpiseq release: tag MPI-dependent archives `-seq` instead of `-<mpi>`, force PT-Scotch off, install libmpiseq |

Target/prefix configuration is not a cache option — it comes from the
staging step (`--target`) via the generated target-config file in the
staging tree.

## Multifloats acquisition (multifloats target only)

Three modes, selected in `cmake/CMakeLists.txt`:

1. **Default — upstream binary release.** Downloads the pinned
   multifloats release (currently `v0.8.2`) and uses its non-LTO static
   library; the Fortran module is compiled from the release's
   pre-expanded source as an eplinalg-owned target.
   - `MULTIFLOATS_VERSION` — release tag; `latest` resolves live via the
     GitHub API (falls back to the pin).
   - `MULTIFLOATS_ARCHIVE` — local path / custom URL to a prefetched
     bundle (offline builds); overrides the GitHub asset URL.
2. **`MULTIFLOATS_FROM_SOURCE=ON`** — git-source FetchContent build
   (`MULTIFLOATS_GIT_REPO` / `MULTIFLOATS_GIT_TAG`).
3. **`MULTIFLOATS_DIR=/path`** — local source checkout via
   `add_subdirectory`.
