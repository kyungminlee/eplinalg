# Installation

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

Next: [usage.md](usage.md) — linking against the installed archives and
the runtime init contract.

## Building from source

If no prebuilt combination fits (different compiler version, another
platform, patched sources), build the stack yourself: see the developer
guide, starting at [../dev/build.md](../dev/build.md).
