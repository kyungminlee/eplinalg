# Usage Guide

**fortran-migrator** provides a multi-step pipeline to migrate numerical Fortran/C libraries to extended precision.

The tool is primarily used through the Python engine (`migrator`) in `src/`.

## Prerequisite: uv

This project uses `uv` for Python dependency management. Ensure you have it installed:
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

## Basic Workflow

The most common way to use the tool is the `run` command, which executes the full pipeline:
1.  **Migrate**: Rewrite source files to the target precision.
2.  **Converge**: Run convergence checks (dual-origin verification).
3.  **Verify**: Check the output for residual types, literals, or column-width overflows.
4.  **Build**: Create static libraries from the migrated source.

```bash
cd src
uv run python -m migrator run ../recipes/blas.yaml work/ --target kind16
```

## Target Selection

The `--target` option selects the precision target. It accepts either a built-in name or a path to a target YAML file:

| Target | Description | Prefix (Real/Complex) |
| :--- | :--- | :--- |
| `kind10` | 80-bit extended precision (x86) | `E` / `Y` |
| `kind16` (default) | 128-bit quad precision | `Q` / `X` |
| `multifloats` | Double-double (`float64x2`) | `M` / `W` |

Example:
```bash
uv run python -m migrator migrate ../recipes/blas.yaml output/ --target multifloats
```

## CLI Commands

All commands are run from the `src/` directory.

### `migrate`
Performs the source-to-source rewriting step.

```bash
uv run python -m migrator migrate ../recipes/blas.yaml output/ --target kind16
```
*   `recipe`: Path to the library's YAML recipe.
*   `output_dir`: Where to write the migrated files.
*   `--target`: Target name or path to target YAML (default: `kind16`).
*   `--dry-run`: Show what would be changed without writing files.
*   `--parser`: Parser backend (`flang` or `gfortran`); omit for regex-only.
*   `--parser-cmd`: Explicit path to the parser compiler binary.

### `verify`
Performs heuristic checks on migrated source files to identify potential issues.

```bash
uv run python -m migrator verify output/
```
Checks for:
*   Residual precision types (e.g., `DOUBLE PRECISION` that wasn't converted).
*   Residual `D`-exponent literals (e.g., `1.0D0`).
*   Column-width overflows in fixed-form code (lines exceeding 72 characters).

### `build`
Compiles and archives the migrated files into static libraries. It produces two libraries:
1.  `lib<prefix><library>.a`: Contains precision-specific routines (e.g., `libqblas.a`).
2.  `lib<library>_common.a`: Contains precision-independent routines.

```bash
uv run python -m migrator build ../recipes/blas.yaml output/ --target kind16 --fc gfortran
```
*   `--fc`: Path to the Fortran compiler (defaults to `gfortran`).

### `diverge`
Reports every co-family pair whose migrated text differs (both halves migrated in memory with a heavy canonicalizer). Useful for preliminary exploration.

```bash
uv run python -m migrator diverge ../recipes/blas.yaml --target kind16
```
*   `--grep`: Regex to only show entries with diff matching this pattern.
*   `--exclude`: Regex to drop entries whose diff matches this pattern.
*   `--context`: Max diff lines per entry (default: 8).
*   `--full`: Print full diff per entry.

### `converge`
Post-migration verification: reads each pair's canonical from disk, re-migrates the S/C sibling in memory, and compares with a light normalizer. This is the authoritative convergence check.

```bash
uv run python -m migrator converge ../recipes/blas.yaml output/ --target kind16
```
*   Same filtering options as `diverge` (`--grep`, `--exclude`, `--context`, `--full`).

See [Convergence Testing](#convergence-testing) for the underlying methodology.

### `run`
Runs the full pipeline: migrate → converge → verify → build.

```bash
uv run python -m migrator run ../recipes/blas.yaml work/ --target kind16
```

## Convergence Testing

Convergence testing leverages the dual-origin nature of numerical libraries. Almost every routine exists in both single-precision (S/C) and double-precision (D/Z) variants. When both are migrated to the same target, the output should be identical:

```
sgemm.f ──migrate(kind16)──→ qgemm_from_s.f ─┐
                                               ├── diff
dgemm.f ──migrate(kind16)──→ qgemm_from_d.f ─┘
```

The migrator normalizes all precision-dependent constructs (type declarations, literals, intrinsics, routine names), so identical output confirms correct handling on both paths.

**What mismatches reveal:**
*   **Algorithmic differences**: Genuinely different implementations between S and D halves (intentional upstream).
*   **Hardcoded constants**: Precision-dependent tolerances (e.g., `1.0E-6` vs `1.0D-12`).
*   **Missing conversions**: A residual `DBLE()` or unrenamed routine name reveals a migrator bug.

The `converge` command automates this for all co-family pairs and generates a summary report.
