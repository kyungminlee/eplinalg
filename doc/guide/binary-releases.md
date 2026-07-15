# Binary Releases (Consumer Guide)

Every tagged release publishes prebuilt static-library archives for
Linux x86-64 on the GitHub Releases page.

## Choosing an archive

Per-combination archives are named

```
<target>-<compiler>-<mpi>-linux-x86_64.tar.gz
```

plus one combined `eplinalg-vX.Y.Z-linux-x86_64.tar.gz` containing every
combination (ABI tags in the library filenames keep them apart).

| axis | values | pick by |
|---|---|---|
| target | `kind16` (q/x), `kind10` (e/y), `multifloats` (m/w) | the arithmetic you want (see below) |
| compiler | `gfortran-12`, `gfortran-15` | your Fortran compiler — `.mod` files and the Fortran ABI are compiler-version-specific |
| mpi | `openmpi`, `mpich`, `intelmpi`, `seq` | your MPI; `seq` bundles the `libmpiseq` stub for single-process use with no MPI installed |

Targets:

| target | real | complex | prefixes |
|---|---|---|---|
| `kind16` | `REAL(KIND=16)` | `COMPLEX(KIND=16)` | `q` / `x` |
| `kind10` | `REAL(KIND=10)` | `COMPLEX(KIND=10)` | `e` / `y` |
| `multifloats` | `TYPE(real64x2)` | `TYPE(cmplx64x2)` | `m` / `w` |

The multifloats types are double-double. In Fortran they come from the
`multifloats` module as `real64x2` / `cmplx64x2`; the C/C++ bridge
(`multifloats_bridge.h`) exposes the same layouts as `float64x2` /
`complex64x2`.

Only MPI-facing archives (BLACS, MUMPS + its glue, PT-Scotch, the MPI
bridges) are MPI-ABI-bound; the rest depend only on the compiler. The
tags are baked into the filenames — `libqxblas-gfortran-15.a` vs
`libqxmumps-gfortran-15-intelmpi-<ver>.a` — so several flavors can be
installed into one prefix.

## Installing

```bash
tar xzf kind16-gfortran-15-openmpi-linux-x86_64.tar.gz -C /usr/local
```

Layout: archives in `lib/`, Fortran modules in `lib/fmod/<mod-tag>/`,
CMake package configs in `lib/cmake/`.

## CMake consumers

Each library is its own package — `find_package(qxblas)`,
`find_package(eylapack)`, `find_package(mwmumps)`, … — exporting
`eplinalg::`-namespaced targets (`eplinalg::qxblas`,
`eplinalg::mwmumps`, and a full-closure convenience target
`eplinalg::<pair>mumps_full`). The generated Config detects the
consumer's Fortran compiler and MPI and loads the matching ABI-tagged
Targets file, so multiple installed flavors coexist and the right one is
selected automatically. Each package chains its own link closure via
`find_dependency`; list every package whose routines you call directly.

```cmake
find_package(qxmumps REQUIRED)
target_link_libraries(app PRIVATE eplinalg::qxmumps_full)
```

## Manual link lines

The archives have circular references — group them:

```bash
gfortran main.f90 \
  -Wl,--start-group /usr/local/lib/lib*-gfortran-15*.a -Wl,--end-group \
  -lquadmath -lstdc++ -lpthread    # plus your MPI's Fortran/C libraries
```

With MKL in the same executable, follow
[mkl-coexistence.md](mkl-coexistence.md): both link orders work; wrap
the MKL group in `-Wl,--push-state,--no-as-needed … -Wl,--pop-state`;
let MKL provide the standard-precision symbols instead of the shipped
Netlib archives; mixed dz+extended MUMPS also needs
`libesmumps_mumps.a` and `libmetis_mumps.a`.

## Runtime contract

Once per process, after `MPI_Init` and before the first extended call:

- multifloats (`m`/`w`): `multifloats_mpi_init()`
- kind16 (`q`/`x`): `quad_mpi_init()`
- kind10 (`e`/`y`): nothing

These register the custom MPI datatypes and reduction ops; skipping the
call crashes the first extended collective.

## Shared libraries

Releases ship static archives only, but every object is compiled
`-fPIC`. To repackage as `.so`s use
`scripts/repackage_shared_libs.sh <pair> <install-lib-dir> <out-dir>`,
which applies the mandatory link rules (`-Wl,--no-define-common`,
`-Wl,-z,now`, export control) — see the repackaging section of
[mkl-coexistence.md](mkl-coexistence.md) for why each rule exists.
