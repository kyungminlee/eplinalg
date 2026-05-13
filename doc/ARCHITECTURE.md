# Architecture

**fortran-migrator** is an automated source-to-source translation tool
for numerical libraries. It produces extended-precision (kind10,
kind16) and double-double (multifloats) variants of the standard
Netlib BLAS / LAPACK / ScaLAPACK / BLACS / PBLAS / PBBLAS / PTZBLAS
ecosystem, plus XBLAS and MUMPS, by rewriting source while preserving
formatting, comments, and preprocessor directives. Output libraries
co-exist with the unmigrated standard-precision archives in a single
build.

## High-level pipeline

```
recipes/<lib>.yaml           targets/<target>.yaml
        │                            │
        └────────────┬───────────────┘
                     ▼
            migrator stage <dir> --target <target> --libraries ...
                     │
                     ▼
   ┌─────────────────────────────────────────────────────────┐
   │  /tmp/stage/                                             │
   │  ├── <lib>/src/        (migrated sources, target prefix)│
   │  ├── _<lib>_src/       (read-only copy of upstream)     │
   │  ├── tests/            (differential precision suite)   │
   │  ├── CMakeLists.txt   (generated)                      │
   │  └── target_config.cmake (target-specific defines)     │
   └─────────────────────────────────────────────────────────┘
                     │
                     ▼
            cmake -S … -B build  &&  cmake --build build
                     │
                     ▼
   ┌─────────────────────────────────────────────────────────┐
   │  build/                                                  │
   │  ├── lib<prefix>blas-<compiler>.a                        │
   │  ├── lib<prefix>lapack-<compiler>.a                      │
   │  ├── lib<prefix>scalapack-<compiler>.a                   │
   │  ├── …                                                  │
   │  ├── libblas-<compiler>.a   (standard-precision archive)│
   │  ├── libreflapack_quad.a    (kind16 reference, tests)   │
   │  └── tests/<lib>/test_*                                 │
   └─────────────────────────────────────────────────────────┘
```

`<prefix>` is one of `q` / `e` / `m` (real) for kind16 / kind10 /
multifloats, paired with `x` / `y` / `w` (complex). See
**Targets** below.

## Hybrid migration strategy

Source-level regex rewriting alone is fragile (Fortran's
fixed-form vs free-form, continuation-line semantics, intrinsic vs
local-variable shadowing). A full parse-and-emit toolchain would
destroy the upstream's column layout, comments, and preprocessor
macros, which is unacceptable for a project that wants its diff
against upstream to be readable.

The migrator combines both:

1. **Parse with a real Fortran compiler.** `flang-new -fc1
   -fdebug-dump-parse-tree` (preferred) or `gfortran
   -fdump-fortran-original` (fallback) produces a parse-tree dump.
   The migrator extracts:
   - type declarations (`REAL`, `DOUBLE PRECISION`, `COMPLEX`,
     `KIND=…`)
   - routine definitions (SUBROUTINE / FUNCTION names + scope)
   - call sites (so name-rewriting doesn't accidentally hit
     character-literal lookalikes)
   - intrinsic vs user-defined name shadowing
   - floating-point literal positions and exponent letters
2. **Apply targeted regex rewrites** guided by those facts. The
   parser output is the oracle; the regex pass is the actuator.
   Comments, whitespace, and `#ifdef` blocks are left exactly as
   they appear upstream.
3. **Fixed-form reformatting.** Any line longer than 72 columns
   after substitution gets re-split with continuation markers in
   column 6. Continuation lines are joined into single statements
   for paren-walkers, then re-emitted in upstream layout.

If neither parser is available, the migrator falls back to
regex-only mode (less precise but still produces correct output for
the well-conditioned upstream sources we use).

## Core components

### Engine — `src/migrator/`

Twelve modules, ~9,000 LOC. Acyclic, unidirectional dependency
graph: lower-level modules (`intrinsics`, `target_mode`,
`symbol_scanner`) have no upward dependencies; higher-level
modules (`pipeline`, `__main__`) orchestrate.

| Module | Purpose |
|---|---|
| `__main__.py` | CLI entry. Subcommands: `migrate`, `verify`, `stage`, `build`, `run`, `diverge`, `converge`. `cmd_stage` (305 lines) orchestrates the per-library staging into a build-ready directory; `_generate_cmake` (233 lines) produces the top-level `CMakeLists.txt` template wired to the multifloats `FetchContent`, optional MPI, and per-target prefix selection. |
| `pipeline.py` | Multi-file migration orchestration. `run_migration` is the entry; `_collect_all_symbols` walks the recipe's `depends:` graph to assemble the symbol universe; `_apply_extra_renames` is the single post-classification rename hook (called from four sites: `run_divergence_report`, `run_convergence_report`, `run_c_convergence_report`, `run_migration`). |
| `fortran_migrator.py` | The Fortran source rewriter. Type-declaration rewriting (`DOUBLE PRECISION` → `REAL(KIND=16)` etc.), literal rewriting (`1.0D0` → `1.0E0_16`), intrinsic rewriting (`DBLE(x)` → `REAL(x, KIND=16)`), routine-name rewriting (case-preserving), `XERBLA('DGEMM ', …)` string rewriting, and fixed-form continuation handling. Largest module by far (~2,800 LOC) because the rules are many. |
| `c_migrator.py` | The C source rewriter and clone-template engine, used for BLACS / PBLAS / PBBLAS / XBLAS / scalapack_c. Case-preserving multi-character prefix expansion (legacy DD/ZZ era); single-character rename in the current M/W era. Header-patch and pointer-cast machinery for kind-dependent C type aliasing. |
| `prefix_classifier.py` | The empirical S/D/C/Z family-discovery engine. `_find_families_single_pass` (lines 215-280) iterates **every** position in each upstream symbol where a precision letter could sit, then groups symbols by tagged-pattern equivalence. **Position-agnostic by construction**: as long as both the S-version and D-version of a name (or C and Z) exist in the symbol universe, the family is discovered regardless of where the precision letter sits. This is why `BLAS_SGBMV_X` ↔ `BLAS_DGBMV_X` (precision letter at position 5) renames cleanly without special handling — see `recipes/lapack.yaml`'s declared `xblas` dep. |
| `intrinsics.py` | Pure-data lookup table: `INTRINSIC_MAP` (~250 entries) maps type-specific intrinsic names (`DABS`, `DCONJG`, `DSQRT`, …) to their generic equivalents and flags whether a `KIND=…` argument is needed at the call site. Updated whenever a new generic appears in modern Fortran. |
| `flang_parser.py` | Subprocess wrapper around `flang-new -fc1`'s JSON parse-tree dump. Extracts symbols, types, and call-site positions. |
| `gfortran_parser.py` | Subprocess wrapper around `gfortran -fdump-fortran-original`'s text tree-dump. Used when LLVM Flang isn't available or when its parse fails. |
| `symbol_scanner.py` | Scans source trees (Fortran or C) for SUBROUTINE / FUNCTION / `void foo(...)` definitions; also can run `nm` over compiled archives if a recipe specifies `symbols.method: nm_library`. |
| `config.py` | YAML recipe loader. `RecipeConfig` dataclass defines the per-library schema: `library`, `language`, `source_dir`, `extensions`, `depends`, `skip_files`, `copy_files`, `prefer_source`, `extra_renames`, `local_renames`, `extra_migrate_files`, `extra_c_dirs`, `extra_fortran_dirs`, `c_return_types`, `c_type_aliases`, `c_pointer_cast_aliases`, `header_patches`, `overrides`, `keep_kind_lines`, `module_renames`, `patches`, `asymmetric_patches`, `one_sided_cleanup`. |
| `target_mode.py` | YAML target loader. `TargetMode` dataclass: prefix map, intrinsic mode (`add_kind` for kind10/kind16, `wrap_constructor` for multifloats), known constants, la_constants_map, module_type_names, MPI datatype names. |
| `__init__.py` | Package marker. |

### Recipes — `recipes/<lib>.yaml`

Eleven libraries. Each recipe declares the source layout, file
extensions, dependency chain, skip lists, patches, and
target-mode-specific knobs.

| Recipe | depends | one-line description |
|---|---|---|
| `blas.yaml` | (none) | BLAS Levels 1-3 from Netlib LAPACK 3.12.1's `BLAS/SRC/`. Two patches on `dnrm2`/`dznrm2` (multifloats safe-summation). |
| `xblas.yaml` | blas | Extra-Precise BLAS from Netlib XBLAS 1.0.248. Auto-generated 402-entry `skip_files` list excludes mixed-precision variants unreachable from LAPACK. Header patches inject quad/multifloats typedefs into `blas_enum.h` (the umbrella header `blas_extended.h` and every `*-f2c.c` bridge transitively include it). |
| `lapack.yaml` | blas, **xblas** | LAPACK 3.12.1 `SRC/`. Three `extra_renames` (ILAENV/ILAENV2STAGE/IPARAM2STAGE → `_EP` variants) route around the upstream gate-bug that returns uninitialized for migrated names. The xblas dep (added in commit `9938731`) brings `BLAS_{S,D,C,Z}*_X` symbols into the classifier's universe so position-5 precision-letter renames Just Work. |
| `blacs.yaml` | (none) | BLACS C library from ScaLAPACK 2.2.3. Multifloats-specific overrides (compiled as C++ for operator overloading) live in `recipes/blacs/mfc_overrides/`. |
| `pblas.yaml` | blas, blacs | Parallel BLAS C library. Includes PBLAS PTOOLS sources via `extra_c_dirs`. |
| `pbblas.yaml` | blas, blacs | Parallel-banded BLAS Fortran helpers. |
| `ptzblas.yaml` | blas | Parallel transpose-on-demand BLAS Fortran helpers. Pulls `ZZDOTC`/`ZZDOTU` from `_scalapack_tools_src/` via `extra_symbol_dirs`. |
| `scalapack.yaml` | lapack, blacs, pblas | ScaLAPACK Fortran. Patches for `pdlanhs`/`pzlanhs` (NPROW=1 norm bug, see `doc/UPSTREAM_BUGS-ScaLAPACK.md`). Four `extra_renames` (`PJLAENV`/`PILAENVX` and `PDLAIECTB`/`PDLAIECTL`). |
| `scalapack_c.yaml` | lapack, blacs, pblas, scalapack | ScaLAPACK C-side wrappers. Mirrors `scalapack`'s `PDLAIECTB`/`PDLAIECTL` renames. |
| `mumps.yaml` | blas, lapack, scalapack | MUMPS 5.8.2 sparse direct solver. 32-entry `skip_files` (C headers, GPU code), 8-entry `copy_files` (integer-only helpers), per-line `keep_kind_lines` manifest preserving `DOUBLE PRECISION` declarations that mean "wall-clock seconds" rather than "working precision". |

### Targets — `targets/<target>.yaml`

Three precision modes. Each carries the type names, literal-form
constructor, intrinsic-rewrite mode, prefix map, and C/MPI interop
constants.

| Target | Real prefix | Complex prefix | Real type | Complex type | Module needed | Notes |
|---|---|---|---|---|---|---|
| `kind10` | E | Y | `REAL(KIND=10)` | `COMPLEX(KIND=10)` | none | 80-bit x87 extended; literals end `_10`; intrinsics get `KIND=10` arg appended |
| `kind16` | Q | X | `REAL(KIND=16)` | `COMPLEX(KIND=16)` | none | 128-bit IEEE quad; the differential-precision reference for the test suite |
| `multifloats` | M | W | `TYPE(real64x2)` | `TYPE(cmplx64x2)` | `multifloats` | double-double; ~60 generic-name overloads + 8 operators; custom MPI types `MPI_FLOAT64X2` etc. |

The previous prefix conventions `t/v` and `dd/zz` were used during
multifloats bring-up and are now retired (commit `823e838`); only
`m/w` is active.

### Targets, prefixes, and the symbol classifier

The reason the engine handles non-trivial cases like
`BLAS_DGBMV_X` (precision letter at position 5, after a literal
`BLAS_` prefix) without special handling: the classifier's family
discovery is **position-agnostic**. Given the symbol universe — a
union of the recipe's own source plus everything its `depends:` graph
brings in — the classifier finds *any* position where swapping S↔D
or C↔Z produces a name that also exists in the universe, and
reports both as a precision family. The rename map then maps the
matched letter to the target's prefix at that exact position.

This is why `recipes/lapack.yaml`'s `depends:` includes both
`blas` and `xblas`: the LAPACK iterative-refinement routines
(`*la_*rfsx_extended.f`) call XBLAS's extra-precise routines, and
without xblas in the symbol universe the classifier would not
discover the position-5 family.

A second class of bugs surfaces only in the kind16 *reference*
build: gfortran's `-freal-8-real-16` flag promotes `REAL(KIND=8)` and
`DOUBLE PRECISION` to quad, but not the second argument of
`CMPLX(real, imag, KIND=8)` — that intrinsic still truncates its
output to `REAL(KIND=8)` even when the surrounding storage is 16
bytes. Bugs of this shape live in `tests/lapack/reflapack/overrides/`
and replace the upstream Netlib source in the **reference** library
only (the migrated target, with explicit `KIND=16` everywhere, is
unaffected). One override file — `zgedmd.f90` — currently exists;
the mechanism is documented in `tests/lapack/reflapack/CMakeLists.txt`.

## C migration

C-based libraries (BLACS, PBLAS, PBBLAS, XBLAS, scalapack_c) follow
a clone-and-substitute pattern:

1. Identify the upstream "template" — usually the double-precision
   variant (e.g. `Cdgesd2d.c`).
2. Apply `c_type_aliases` (e.g. `double` → `float128`,
   `float64x2`) and `c_pointer_cast_aliases`
   (e.g. `(double*)` → `(float128*)`).
3. Apply `header_patches` to inject typedefs / extern "C" guards
   into headers that the migrated source includes.
4. Rename routine names case-preservingly (e.g. `Cdgesd2d` →
   `Cqgesd2d` for kind16 real, `Cwgesd2d` for multifloats complex).
5. Update MPI datatype constants (e.g. `MPI_DOUBLE` → `MPI_REAL16`
   or `MPI_FLOAT64X2`).
6. For multifloats only: emit the file as C++ (so operator
   overloading on `real64x2` works); apply `extern "C"` guards
   around any blocks the Fortran callers expect to see.

## Build system — `cmake/`

The migrator does not ship CMakeLists.txt for individual libraries.
Instead, `migrator stage` writes a top-level `CMakeLists.txt` (from
the template in `__main__.py:_generate_cmake`) into the staging
directory, plus a target-specific `target_config.cmake` that sets
the prefix, type names, and MPI-datatype symbols.

The shared infrastructure lives in `cmake/CMakeLists.txt` (~870
lines) and `cmake/FortranCompiler.cmake` (~410 lines). The major
sections:

1. Compiler probes (KIND=10, KIND=16 support) and global flag
   setup.
2. Multifloats `FetchContent` — pinned via `-DMULTIFLOATS_GIT_TAG=…`
   or substituted with `-DMULTIFLOATS_DIR=/local/path` for offline
   builds.
3. `add_standard_fortran_library` / `add_standard_c_library` — build
   the unmigrated upstream archive.
4. `add_migrated_fortran_library` / `add_migrated_c_library` — build
   the per-target migrated archive, with PUBLIC linkage to the
   upstream archive so the migrated routines can call into S/D/C/Z
   helpers (`LSAME`, `XERBLA`, `IEEECK`, etc.) that don't carry a
   precision letter.
5. Per-library wiring: BLAS, XBLAS, BLACS, LAPACK, PTZBLAS, PBBLAS,
   PBLAS, ScaLAPACK, scalapack_c, MUMPS.
6. **libmpiseq** — sequential MPI stub built from MUMPS 5.8.2's
   `libseq/`. Provides `mpi_init_`, `mpi_send_`, `blacs_pinfo_`,
   etc. so single-rank tests can link cleanly without a real MPI.
   Conditionally compiled if `_mpiseq_src/` is present in the stage.
7. Multifloats helper modules: `la_constants_mf`, `la_xisnan_mf`,
   `multifloats_mpi`. Wrapped in `$<BUILD_INTERFACE:>` generator
   expressions to avoid leaking into downstream consumers.
8. Install rules — per-target export sets, Config.cmake files
   keyed on compiler tag (gfortran-13, ifx-2024.x, etc.) so
   downstream `find_package(qlapack)` fails clearly when the
   consumer's compiler doesn't match an installed build.

## Test infrastructure

`tests/<lib>/` holds differential-precision suites. Each test
constructs the same input matrix at quad precision, runs both the
**reference** (Netlib upstream compiled with `-freal-8-real-16`,
producing `libreflapack_quad.a`) and the **migrated target** library
on it, and compares outputs. The reference and target should agree
to within ~33-34 decimal digits for kind16, ~19 digits for kind10,
~30 digits for multifloats.

Test wrappers (`tests/<lib>/common/target_*.fypp`) bridge between
the standard reference signature and the target-specific symbol
names. The wrappers are emitted once per target and used by every
test driver.

## Out of band: upstream bug tracking

Upstream Netlib bugs that require in-tree workarounds are tracked
in `doc/UPSTREAM_BUGS-{LAPACK,ScaLAPACK,MUMPS}.md`. Two mechanisms
route around them:

- `patches` (recipe field): a unified-diff `<basename>.patch` file
  under `recipes/<lib>/patches/` is applied to the staged upstream
  source before migration. The patched body still goes through the
  normal pipeline so it emits the per-target migrated name. Used
  for ScaLAPACK `pdlanhs.f` / `pzlanhs.f` (NPROW=1 norm bug),
  LAPACK `cheevx.f` (LWKOPT min guard), and dozens of others.
- `tests/lapack/reflapack/overrides/<file>.f90`: replace the
  upstream source in the **reference** build only. Used for
  `zgedmd.f90` (gfortran `-freal-8-real-16` `CMPLX(…,KIND=8)`
  truncation).

## CI — `.github/workflows/release.yml`

Triggered on `v*` tag push or manual dispatch. Three jobs:

1. **stage** (per target) — runs `migrator stage` and uploads the
   staged tree as an artifact.
2. **build** (45 combinations: 3 targets × 5 compilers ×
   3 MPI implementations, minus impossible combos like kind10+ifx
   and kind16+flang) — downloads the staging artifact, configures
   with CMake, builds, installs, packages a tarball.
3. **release** — only on tag push; collects all build tarballs into
   a GitHub Release.

The CI pipeline does **not** currently run ctest. Migration
correctness is verified only on developer machines via local
ctest runs.

## Repository layout

```
fortran-migrator/
├── src/migrator/        # The migration engine (12 modules, ~9k LOC)
├── recipes/             # 11 library recipes
├── targets/             # 3 precision targets
├── tests/               # 10 differential-precision suites (1051 tests)
├── external/            # Vendored upstream: LAPACK 3.12.1, ScaLAPACK 2.2.3,
│                        # MUMPS 5.8.2, XBLAS 1.0.248, multifloats, multifloats-mpi,
│                        # impi-headers
├── cmake/               # Shared CMake infrastructure
├── doc/                 # User and architectural documentation
│   ├── ARCHITECTURE.md  # This file
│   ├── USAGE.md         # CLI reference and worked examples
│   ├── RECIPES.md       # Recipe YAML schema
│   ├── PROCEDURES.md    # Generated routine cross-reference
│   ├── INTRINSICS.md    # Generic-intrinsic table (manual reference)
│   ├── UPSTREAM_BUGS*.md# Tracked upstream bugs (LAPACK / ScaLAPACK / MUMPS)
│   ├── DEVELOPER.md     # Developer onboarding
│   ├── NOTE.md          # Misc design notes
├── scripts/             # Manual helpers (compile_*.sh, sweep tools)
├── tools/               # Build artifacts (gen_procedures.py)
├── pyproject.toml       # uv-managed; runtime: pyyaml, tqdm; dev: fypp, pytest
└── README.md
```

## Design principles

- **Don't fight the upstream's layout.** The migrated source's diff
  against the upstream should be small and obvious (precision-related
  only). This rules out reformatting passes that re-flow comments or
  rearrange continuation lines.
- **Preserve the build chain.** Migrated archives co-exist with the
  unmigrated archive; PUBLIC linkage handles symbol visibility for
  helpers like `LSAME` and `XERBLA` that aren't precision-renamed.
- **Push fixups out of the engine, into recipes.** When upstream has
  a bug that affects only one half of a precision pair, prefer a
  `patches/<basename>.patch` over engine special-cases. When a
  routine needs a different rename (e.g., extended-precision suffix
  `_EP`), prefer `extra_renames` over hard-coding.
- **Discovery, not configuration.** The classifier discovers
  precision families by scanning the symbol universe. Recipes
  rarely need to declare per-routine prefixes — the classifier
  finds them.
- **Reference-side patches are a last resort.** The
  `tests/lapack/reflapack/overrides/` mechanism exists for the
  narrow class of bugs that surface only under
  `-freal-8-real-16` promotion of upstream code. The migrator's
  output is correct; the reference build is the one that needs
  patching.
