# User guide

Obtaining and using the extended-precision libraries.

- [installation.md](installation.md) — release archives, choosing a variant, install layout
- [usage.md](usage.md) — naming scheme, linking (CMake and manual), runtime MPI init contract
- [mkl-coexistence.md](mkl-coexistence.md) — MKL LP64 in the same executable; shared-library repackaging rules
- [api/](api/index.md) — API reference

The hand-written parallel BLAS overlay lives in the separate
[epblas-parallel](https://github.com/kyungminlee/epblas-parallel)
project; the `<prefix>blas` archives here are the plain serial
migrated BLAS.
