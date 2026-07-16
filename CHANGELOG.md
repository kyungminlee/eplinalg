# Changelog

Release notes live on the [GitHub releases page](https://github.com/kyungminlee/eplinalg/releases);
this file summarizes each tagged version. The current version is in
[VERSION](VERSION).

## v0.15.0

- Repository restructured to the standard package layout: `extern/`
  (was `external/`), `codegen/{migrator,recipes,targets}/`,
  `src/` (first-party runtime, was `runtime/`),
  `test/{unit,integration}/`, `example/`, root `CMakePresets.json`.
  The staged-tree layout is unchanged, so consumers of the release
  archives are unaffected.
- Version single-sourced from the root `VERSION` file (read by
  `pyproject.toml` and the staged CMake build).
- New root files: `LICENSE` (MIT), `CHANGELOG.md`, `CONTRIBUTING.md`,
  thin root `CMakeLists.txt` with an optional Sphinx `doc` target,
  `.editorconfig`, `.clang-format`, `.clang-tidy`.

## v0.14.0

- Multifloats MPI reduce ops renamed `MPI_MM_*` / `MPI_WW_*`
  (`MPI_DD_*` / `MPI_ZZ_*` kept as deprecated aliases for one release
  cycle).
- Documentation reorganized into `doc/user/` and `doc/dev/`.

## v0.13.x

- v0.13.1: shared-library repackaging requires `-Wl,-z,now`; reference
  repackaging script.
- v0.13.0: position-independent archives — static libraries can be
  repackaged as shared objects.

## v0.12.x

- v0.12.1: extended-precision PBLAS typeset SIGSEGV fix
  (hidden-CHARACTER-length ABI trampolines).
- v0.12.0: MKL-coexistence `ep_` symbol privatization + Intel MPI
  non-commutative reduce-op guard.

## v0.11.0

- Netlib-pristine public C headers via precision-prefixed siblings.

## v0.10.0

- Multifloats acquisition defaults to the upstream binary release.

## v0.9.x

- v0.9.1: `float64x2` diagnostic-print robustness + seq datatype
  defaults.
- v0.9.0: vendored MUMPS upgraded 5.8.2 → 5.9.0; all 10 arithmetics ×
  openmpi/intelmpi/seq validated.

## v0.8.0

- 10-precision linalg stacks, seq consumer path, clean module install,
  ordering-library headers.

## Earlier

v0.1.0 – v0.7.1: incremental build-out of the migrator, the kind16 /
kind10 / multifloats stacks, packaging, and `find_package` support. See
the tag annotations (`git tag -n9`) and the releases page.
