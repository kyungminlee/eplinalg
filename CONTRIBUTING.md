# Contributing

Developer documentation lives under [doc/dev/](doc/dev/index.md):

- [configure.md](doc/dev/configure.md) — CMake options, presets, multifloats acquisition
- [build.md](doc/dev/build.md) — the stage → cmake two-step
- [codegen.md](doc/dev/codegen.md) — running and modifying the Python migrator
- [test.md](doc/dev/test.md) — pytest, patch/convergence checks, ctest, release validation
- [architecture.md](doc/dev/architecture.md) — internal design and module map
- [debugging.md](doc/dev/debugging.md) — known pitfalls
- [release.md](doc/dev/release.md) — cutting a release

Quick start:

```sh
uv sync                       # editable install of codegen/migrator
uv run pytest                 # unit tests (test/unit/)
uv run python -m migrator stage staged --target kind16
cmake -S staged -B staged/build -DCMAKE_BUILD_TYPE=Release
cmake --build staged/build -j
```

The single source of truth for the version is the [VERSION](VERSION)
file; `pyproject.toml` and the staged CMake build both read it.
Vendored sources under `extern/` are never edited — upstream fixes are
applied as recipe patches (`codegen/recipes/<lib>/patches/`).
