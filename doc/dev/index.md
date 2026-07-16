# Developer guide

Developing eplinalg itself. Consuming the produced libraries:
[../user/](../user/index.md).

Lifecycle:

- [configure.md](configure.md) — CMake presets, toolchains, options
- [build.md](build.md) — staging (Python codegen) + building the stack
- [codegen.md](codegen.md) — the migrator: CLI, transformation rules, extending it
- [recipes.md](recipes.md) — library recipe (YAML) reference
- [test.md](test.md) — unit tests, convergence, stack tests, release validation
- [debugging.md](debugging.md) — known pitfalls and debugging knowledge
- [release.md](release.md) — cutting a release
- [architecture.md](architecture.md) — internal design and module map

Reference material:

- [convergence.md](convergence.md) — the dual-origin convergence methodology
- [kind16-divergences.md](kind16-divergences.md) — per-routine KIND=16 divergence analysis
- [../upstream-bugs/](../upstream-bugs/README.md) — vendored-source bug catalogue and workarounds
