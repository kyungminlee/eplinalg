# Codegen: the migrator

The Python engine (`migrator`, in `src/`) rewrites the vendored Netlib
Fortran/C sources to extended precision. It is the code generator of the
project: everything under a staging tree is its output. Design and
module map: [architecture.md](architecture.md). Recipe schema:
[recipes.md](recipes.md).

All commands are run from `src/` with `uv` (see
[build.md](build.md) for the prerequisite).

## Target selection

`--target` accepts a built-in name or a path to a target YAML
(`codegen/targets/*.yaml`):

| Target | Types | Prefix (Real/Complex) |
| :--- | :--- | :--- |
| `kind10` | `REAL(KIND=10)` / `COMPLEX(KIND=10)` | `E` / `Y` |
| `kind16` (default) | `REAL(KIND=16)` / `COMPLEX(KIND=16)` | `Q` / `X` |
| `multifloats` | `TYPE(real64x2)` / `TYPE(cmplx64x2)` | `M` / `W` |

## CLI commands

### `stage`
Migrates **all** libraries into one staging directory with a unified
CMake build — the usual entry point for building the full stack.

```bash
uv run python -m migrator stage /tmp/stage-q --target kind16
```
*   `--libraries`: Subset of libraries to migrate (default: all).
*   `--parser` / `--parser-cmd`: As for `migrate`.
*   The staging tree is self-contained (sources, CMake system, presets,
    tests); build it with plain CMake. Note it *snapshots*
    `test/integration/CMakeLists.txt` and the presets — re-stage after editing them.

### `migrate`
Source-to-source rewriting for one recipe.

```bash
uv run python -m migrator migrate ../codegen/recipes/blas.yaml output/ --target kind16
```
*   `recipe`: Path to the library's YAML recipe.
*   `output_dir`: Where to write the migrated files.
*   `--target`: Target name or path to target YAML (default: `kind16`).
*   `--dry-run`: Show what would be changed without writing files.
*   `--parser`: Parser backend (`flang` or `gfortran`); omit for regex-only.
*   `--parser-cmd`: Explicit path to the parser compiler binary.

### `verify`
Heuristic checks on migrated sources: residual precision types
(`DOUBLE PRECISION` that wasn't converted), residual `D`-exponent
literals (`1.0D0`), fixed-form column-width overflows (> 72).

```bash
uv run python -m migrator verify output/
```

### `diverge`
Reports every co-family pair whose migrated text differs (both halves
migrated in memory with the merged heavy canonicalizer). Honors the
recipe's `expected_divergences:` whitelist; `--no-whitelist` bypasses
it. This is the authoritative convergence check — the earlier separate
`converge` subcommand was merged into it.

```bash
uv run python -m migrator diverge ../codegen/recipes/blas.yaml --target kind16
```
*   `--grep` / `--exclude`: Regex filters on the diff text.
*   `--context`: Max diff lines per entry (default: 8); `--full` for all.
*   `--no-whitelist`: Ignore `expected_divergences:` entries.

Methodology: [convergence.md](convergence.md).

### `build`
Compiles one recipe's migrated files into static archives
(`lib<pair><library>.a` + `lib<library>_common.a`). `--fc` selects the
Fortran compiler. Mostly superseded by `stage` + CMake.

### `run`
Full single-recipe pipeline: migrate → diverge → verify → build.

### `prepare`
Stages upstream sources for a recipe and applies its patch list
(normally run implicitly by the other commands). `--rebuild` wipes and
re-stages even if the cache stamp is fresh.

### `verify-patches`
CI check: every patch touching a precision-prefixed file must touch all
four siblings, unless listed under `asymmetric_patches:` in the recipe.

## What gets transformed (Fortran)

The tool preserves comments, preprocessor directives, and fixed-form
column layout exactly; only precision-bearing constructs change:

1. **Type declarations** — `REAL`, `REAL*8`, `DOUBLE PRECISION`,
   `COMPLEX*16`, `DOUBLE COMPLEX`, … all map to the target type. Both
   the single- (`s`/`c`) and double-precision (`d`/`z`) variants of a
   routine map to the *same* target.
2. **Routine names** — the precision prefix is swapped (`DGEMM` →
   `QGEMM`, `PZGESV` → `PWGESV`) in definitions, `CALL` sites,
   `EXTERNAL` declarations, expression references, and string literals
   passed to `XERBLA`. Non-prefixed type-independent helpers (`LSAME`,
   `XERBLA`, `ILAENV`) are never renamed; the symbol scanner
   (`symbol_scanner.py`) + prefix classifier (`prefix_classifier.py`)
   only rename symbols with confirmed precision-variant siblings.
3. **Literals** — `1.0D+0` → `1.0E+0_16` (KIND targets) or a
   constructor/module-constant form (multifloats). The `known_constants`
   mechanism replaces local `PARAMETER` constants with module-provided
   values where the target defines them.
4. **Intrinsics** — type-specific intrinsics become generic or
   KIND-explicit forms (`DBLE(x)` → `REAL(x, KIND=16)`, `DIMAG` →
   `AIMAG`, `DABS` → `ABS`, …). Catalog: `codegen/migrator/intrinsics.py`
   and [../user/api/intrinsics.md](../user/api/intrinsics.md).
5. **Machine parameters** — `DLAMCH`-family routines are regenerated,
   not renamed: the constants they return are precision-dependent.
6. **File names** — `dgemm.f` → `qgemm.f`, `pdgesv.f` → `pmgesv.f`.

### Hybrid parser + regex strategy

Pure regex is brittle on fixed-form Fortran (continuation lines,
comments, string literals); pure compiler unparse (flang
`-fdebug-unparse`) destroys formatting. So: parse with flang or gfortran
(`--parser`) to locate transformation sites with provenance, then apply
the edits as source-level text patches. When a replacement overflows
column 72, whitespace is squeezed or a continuation line is inserted in
the file's own continuation style. Without `--parser` the regex-only
fallback runs.

### C libraries (BLACS, PBLAS C layer)

Migrated by clone-and-substitute in `c_migrator.py`: mechanical text
substitution on C types, MPI datatypes, and function-name prefixes, with
`[^a-zA-Z_0-9]` boundaries so replacements don't re-match inside
already-substituted names. Public Netlib headers stay pristine;
extended-precision declarations go in precision-prefixed sibling
headers. Fragile spots to know about: `PBTYP_T` function-pointer tables
(a dead entry once caused np≥2 SIGSEGVs) and hidden CHARACTER-length
ABI at Fortran↔C callback boundaries (solved by per-callback C
trampolines, `PBcharshim.h`).

## Edge cases worth knowing

- `EQUIVALENCE` creates type aliasing that can break when type sizes
  change — flagged for manual review.
- Fortran `INCLUDE` is a statement, not a preprocessor directive;
  included files are migrated too.
- Code inside `#ifdef` blocks is transformed; the directives themselves
  are never touched.
- Some upstream routines use `DOUBLE PRECISION` with two intents in the
  same family ("working precision" in `d*` files, "timing/MPI
  interface" in `s*` files) — the divergence machinery is what catches
  such semantic asymmetries; see [convergence.md](convergence.md).

## Modifying the migrator

- Per-library file lists, patches, and exceptions: the recipe YAMLs
  (`codegen/recipes/*.yaml`, schema in [recipes.md](recipes.md)).
- Per-target types, literals, name mappings, `c_interop` symbol maps:
  `codegen/targets/*.yaml`.
- Intrinsic conversion table: `codegen/migrator/intrinsics.py`.
- Fortran rewrite passes: the `codegen/migrator/fortran/` package (one
  module per concern: declarations, literals, renames, MPI calls, …).
- Run the unit suite (`uv run pytest`) and a `diverge` pass after any
  transform change; see [test.md](test.md).
