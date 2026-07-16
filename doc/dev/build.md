# Build

Building is two steps: Python codegen (staging), then CMake.

## Prerequisite: uv

The migrator's Python dependencies are managed with `uv`:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

## Stage, then build

```bash
# 1. codegen: migrate every library into a self-contained CMake tree
cd src
uv run python -m migrator stage /tmp/stage-q --target kind16

# 2. configure + build (presets: see configure.md)
cmake -S /tmp/stage-q --preset linux-impi
cmake --build /tmp/stage-q/build -j8
```

The staging tree carries everything the build needs (migrated sources,
CMake system, presets, tests); the repo itself is never configured.
Migrator CLI details: [codegen.md](codegen.md).

## Generated Fortran modules

`.mod` files land in `fmod/<FORTRAN_MOD_COMPAT_TAG>/` in the build tree
and install under `<libdir>/fmod/<mod-tag>/`. They are
compiler-version-specific — a consumer with a different gfortran major
version needs a rebuild (this is why release archives carry a compiler
tag).

## Caveats

- `stage` **snapshots** `tests/CMakeLists.txt` and the preset file into
  the staging tree. After editing them in the repo, re-stage —
  rebuilding an old tree silently uses the stale copies.
- A full build tree per target is large; building all targets means
  three full linalg stacks. Remove each build tree before starting the
  next on small disks.
