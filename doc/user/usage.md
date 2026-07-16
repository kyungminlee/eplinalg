# Using the libraries

How to link and call the installed archives. Archive selection and
install layout: [installation.md](installation.md).

## Naming scheme

Routine names follow the Netlib convention with the precision prefix
swapped: `DGEMM` → `QGEMM` (kind16 real), `PZGESV` → `PWGESV`
(multifloats complex), `DMUMPS` → `MMUMPS`, and so on:

| target | real prefix | complex prefix |
|---|---|---|
| `kind16` | `q` | `x` |
| `kind10` | `e` | `y` |
| `multifloats` | `m` | `w` |

Library names carry the pair prefix (`libqxlapack.a`,
`libmwscalapack.a`); precision-independent engine code lives in
`*_common` archives with `ep_`-prefixed entry points (grid/descriptor
management etc. — see [mkl-coexistence.md](mkl-coexistence.md)).

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

## API reference

[api/](api/index.md) — routine cross-reference, MUMPS API constraints,
intrinsic surface of the extended types.
