# Handoff — Comprehensive code review findings for eplinalg

**Date:** 2026-05-28
**Repo:** `/home/kyungminlee/code/eplinalg`
**Branch:** `develop` @ `c344cb19` (matches `origin/develop`)
**Status:** Review complete; nothing fixed yet. Awaiting per-item authorization to apply.

## What just happened

A comprehensive correctness + simplification review was run, separate from
the rename-residue cleanup earlier today (see git log and the now-empty
`/tmp/handoff-pre-split-residue-20260528.md`). Five parallel subagents
covered (1) Fortran/C migrator core, (2) migrator support modules,
(3) CLI/staging, (4) CMake build system, (5) tests/CI/stubs. Each finding
below was either spot-verified by reading the cited line or empirically
reproduced — verification status is recorded per item.

**Out of review scope:** `external/`, `build/`, pre-migrated outputs under
`blas/`, `lapack/`, `scalapack*/`, `blacs/`, `mumps/`, `ptzblas/`,
`pbblas/`, `pblas/`, `xblas/src/` (these are generated artifacts, not
hand-written code).

**Verification key:** ✓ verified by code read or empirical Python/CLI
reproduction · △ plausible mechanism, runtime impact unconfirmed.

## Pending unrelated state (carry forward)

- 10 rename-residue edits already applied earlier this session, uncommitted
  on working tree. See `git status` — covers README.md title, doc/guide/*
  bold names, src/migrator/__init__.py docstring, tests/mumps/README.md,
  tools/gen_procedures.py ROOT path (now `Path(__file__).resolve().parents[1]`),
  scripts/run_perf_n512.sh (deleted per user's choice).
- Project norm: NO commits without explicit user authorization.

---

## CORRECTNESS — Tier 1 (build/numerical correctness)

### C1 — MPI byte-size sentinels wrong for kind10 — `mpiseq_c_stubs.c:215,220` ✓

```c
if (t == MPI_LONG_DOUBLE)            return 10;   /* x87 80-bit */
if (t == MPI_C_LONG_DOUBLE_COMPLEX)  return 20;
```

Comment says "x87 80-bit" — but on x86_64 Linux `sizeof(long double)` is
16 bytes (verified: compiled `printf("%zu", sizeof(long double))` →16).
The Fortran-side `REAL(KIND=10)` storage matches that 16-byte layout.

`targets/kind10.yaml:32-33` maps:
```yaml
mpi_real:    MPI_LONG_DOUBLE
mpi_complex: MPI_C_LONG_DOUBLE_COMPLEX
```

The function's return value is encoded into the libmpiseq sentinel
(`0x10000000 | total_size_in_bytes`, comment at L227-229) which
`MUMPS_COPY` in libseq's `mpi.f` uses to walk byte buffers. Under-reporting
size 10 vs actual 16 means MUMPS_COPY truncates each element by 6 bytes
on every walk → silent data corruption at kind10.

**Fix:** return 16 and 32 to match the C ABI. Update the comment to
clarify that the value is byte-stride, not bit-precision.

### C2 — mumps_full LINK_GROUP uses wrong namespace — `cmake/CMakeLists.txt:1083` △

`project(${TARGET_NAME}_libraries …)` at L29 sets `PROJECT_NAME` to
`kind10_libraries` (or kind16/multifloats). The LINK_GROUP at L1083 uses
`${PROJECT_NAME}::emumps_c` / `${PROJECT_NAME}::emumps`. But the install
side (`_install_library_pair`, L1381) writes the targets file with
`NAMESPACE eplinalg::`. The serialized Targets file then references
`eplinalg::emumps_full` whose LINK_GROUP body names `kind10_libraries::`
imports that don't exist on the consumer side.

**Fix:** literal `eplinalg::` (or thread NAMESPACE as a variable through
both code paths). Verify with a downstream `find_package(emumps)` smoke
test.

### C3 — multifloats_mpi PUBLIC-links unexported target — `cmake/CMakeLists.txt:235` △

```cmake
target_link_libraries(multifloats_mpi … PUBLIC multifloats)
```

`multifloats` comes from `FetchContent_Populate` at L203 with
`EXCLUDE_FROM_ALL`, so it has no install rule. Compare with L405 which
correctly wraps `$<BUILD_INTERFACE:multifloatsf>` to keep the
build-private dependency out of the exported INTERFACE_LINK_LIBRARIES.

When `fortran_install_library(multifloats_mpi … MPI)` runs at L1445,
CMake errors with `install(EXPORT) includes target multifloats_mpi
which requires target multifloats that is not in any export set`.

**Fix:** wrap the L235 link in `$<BUILD_INTERFACE:multifloats>` mirroring
the L405 pattern.

### C4 — `verify` subcommand has no `--target` flag — `src/migrator/__main__.py:1717-1720` ✓

```python
p = sub.add_parser('verify', help='Verify migrated output')
p.add_argument('output_dir', type=Path)
p.set_defaults(func=cmd_verify)
```

No call to `_add_target_args(p)`. `cmd_verify` (L388) calls
`_get_target_mode(args)`, which at L40 falls back to `'kind16'` when
`args.target` is unset. Result: `migrator verify stage-kind10/` runs
kind16 rules, false-passing on multifloats `_wp` residuals and
false-failing on kind10 column-overflow expectations.

**Fix:** add `_add_target_args(p)` between L1719 and L1720.

### C5 — std-dir patching asymmetry — `src/migrator/__main__.py:1455-1491` ✓

`cmd_stage` stages 10 standard-precision sibling source dirs
(`_blacs_src`, `_pblas_src`, `_ptzblas_src`, `_pbblas_src`,
`_scalapack_src`, `_scalapack_tools_src`, `_scalapack_redist_src`,
`_mpiseq_src`, `_mumps_upstream_src`, `_mumps_upstream_include`) directly
from `external/` with no patch application:

```python
src = proj_root / 'external' / rel_src
…
shutil.copytree(src, dst)
```

But `_stage_baseline` (L1586-1599) for kind4/kind8 routes those same
dirs through `run_prepare` when a recipe exists, applying recipe patches.

Result: 7-of-10 std archives (those with recipes — blacs, pblas, ptzblas,
pbblas, scalapack, scalapack_tools, scalapack_c, mumps) ship unpatched
inside the migrated stage but patched inside the kind4/kind8 baseline
stage. Any "symmetric" recipe patch (the default contract per
`cmd_verify_patches`) silently fails to apply to std archives that share
a precision target.

**Fix:** route `_std_dirs` through `run_prepare` for entries that have
a matching recipe; pass them through `_staged_or_external` like the
baseline does. Or document in `asymmetric_patches` that the migrated
stage's std archives should not get the patch.

### C6 — CPP guard stack ignores `#else`/`#elif` — `src/migrator/fortran_migrator.py:1898` △

`convert_parameter_stmts`' `cpp_stack` (L1896-1905) tracks `#if`/`#ifdef`
push, `#endif` pop, but not `#else`/`#elif`. A PARAMETER inside the
`#else` branch of an `#if`/`#else`/`#endif` gets its converted runtime
assignment wrapped under the outer `#if` guard only — fires in the wrong
branch.

**Trigger:** any LAPACK/MUMPS file with
```fortran
#if FOO
      DOUBLE PRECISION ZZ
#else
      DOUBLE PRECISION ZZ
      PARAMETER(ZZ = 1.0D0)
#endif
```
Post-rewrite: `#if FOO \n ZZ = 1.0D0 \n #endif` — assignment in the
`#if FOO` true branch where the original PARAMETER did not exist.

**Fix:** push a "flipped" sentinel onto `cpp_stack` on `#else`/`#elif`
so the head guard gets negated; pop it on the corresponding `#endif`.

### C7 — `break` after unbalanced-paren skips all later intrinsic calls — `src/migrator/fortran_migrator.py:866-871` △

In `replace_intrinsic_calls`, the else-branch on the unbalanced-paren
case renames the intrinsic in-place then `break`s the inner while loop,
abandoning any subsequent occurrences of that intrinsic on the same line.

**Trigger:** lines where the first intrinsic call has malformed parens
(rare in clean source; observed in fixed-form lines bisected by the col-72
cutoff and in comment-mangled inputs).

**Fix:** replace `break` with `search_start = m.end(); continue`.

### C8 — `_rewrite_int_of_complex` emits `real()` without USE-only entry — `src/migrator/fortran_migrator.py:1499` △

`_rewrite_int_of_complex` rewrites `INT(zsomething)` → `INT(real(inner))`
for multifloats targets, but `_build_use_only_clause` only predictively
adds `dble` (L2502-2509) when INT-on-real64x2 is detected. It does NOT
add `real` when INT co-occurs with any of the procedure's complex names.

**Trigger:** any `z*.f` LAPACK/BLAS file calling `INT(ZSCALE)` or similar
under multifloats target. Post-rewrite: gfortran refuses because `real`
is not a USE-imported generic.

**Fix:** extend `_build_use_only_clause` to add `real` when both `int(`/
`nint(` and a complex_name co-occur in the procedure body.

### C9 — Array constructor's first literal is dropped — `src/migrator/prepare.py:44-46` ✓

```python
_BARE_LITERAL_RE = re.compile(
    r'(?<![\w.\[])(\d+\.\d*|\d*\.\d+)(?![DdEe\w]|_\d|_[A-Za-z])'
)
```

The `[` in the lookbehind character class makes the first literal in a
F90 array constructor `[5.0, 6.0]` invisible to the regex. Empirically:
`re.findall(_BARE_LITERAL_RE, '[5.0, 6.0]')` → `['6.0']` (NOT both).

The D-half normalizer then rewrites that to `[5.0, 6.0D+0]` — asymmetric
kind-suffix promotion within the same constructor.

**Fix:** drop `\[` from the lookbehind char class. `(` isn't excluded
either and `FOO(1.0)` already works fine; F90 array constructors `[…]`
have no syntactic difference that justifies the exclusion.

### C10 — Empty-stderr used as success signal — `scripts/compile_lapack.sh:41` (also `compile_pblas.sh:55`, `compile_scalapack.sh:68`) △

```sh
out=$(gfortran … 2>&1 >/dev/null)
if [ -z "$out" ]; then
    pass=$((pass+1))
else
    echo "$f" >> /tmp/lapack_failures.txt
fi
```

Detects compile success by empty stderr, not by `$?`. gfortran routinely
emits warnings (`-Wunused`, `-Wmaybe-uninitialized`, legacy-feature) on
Netlib-style fixed-form code, so warning-only files are logged as
failures. C/C++ compile blocks in the same scripts use the correct
`if [ $? -eq 0 ]` pattern.

**Fix:** check `$?` exit code, log stderr separately if non-empty.

---

## CORRECTNESS — Tier 2 (silent edge cases, low frequency)

### C11 — `END  SUBROUTINE FOO` (2 spaces) matches as new definition — `src/migrator/symbol_scanner.py:16-20` ✓

```python
_FORTRAN_DEF_RE = re.compile(
    r'(?<!\w)(?<!END\s)(?<!end\s)'
    r'(?:SUBROUTINE|FUNCTION|MODULE…)\s+([A-Za-z]\w*)',
    re.IGNORECASE,
)
```

Python's negative lookbehind is fixed-width; `(?<!END\s)` only catches a
single whitespace char. Files in `external/MUMPS_5.8.2/src/cfac_*.F`,
`dfac_*.F`, `tools_common.F:107`, `mumps_pivnul_mod.F:71` contain literal
`END  SUBROUTINE NAME` with two spaces. Empirically:
`_FORTRAN_DEF_RE.findall('END  SUBROUTINE FOO')` → `['FOO']`.

In MUMPS the affected sites happen to name routines that exist anyway,
so symbol set is correct by coincidence — but the regex would also
match `END  MODULE OTHER` typos and inject non-existent symbols.

**Fix:** explicitly strip `^\s*END\s+` then match the bare keyword form,
or use a structural pre-scan to skip lines starting with `END`.

### C12 — `_USE_STMT_RE` misses F2003 forms — `src/migrator/gfortran_parser.py:142` ✓

`_USE_STMT_RE` is `^\s+USE\s+(\w+)`. Empirically:
- `'    USE foo'` → matches ✓
- `'    USE :: foo'` → `None` ✗
- `'    USE, intrinsic :: iso_c_binding'` → `None` ✗

If any gfortran dump preserves `USE ::` form (or any input file uses the
F2003 syntax that survives into the dump), those USE associations are
dropped from `facts.use_stmt_ranges`.

**Fix:** `^\s*USE(?:\s*,[^:]*)?\s*(?:::)?\s*(\w+)`.

### C13 — `TYPE IS` / `TYPE DEFAULT` matched as derived-type definitions — `src/migrator/symbol_scanner.py:28-31` △

`_FORTRAN_TYPE_DEF_RE` allows trailing `\s*(?:!.*)?$` after the captured
name. For `TYPE IS` (continuation onto next line) or `TYPE DEFAULT`
inside a SELECT TYPE block, it captures `IS` or `DEFAULT` as the type
name.

**Fix:** explicitly exclude `IS|DEFAULT` from the capture group.

### C14 — `--libraries <typo>` silently produces empty stage — `src/migrator/__main__.py:1188-1192` ✓

```python
if args.libraries:
    lib_set = set(args.libraries)
    libraries = [(n, r) for n, r in LIBRARY_ORDER if n in lib_set]
```

No diagnostic for `set(args.libraries) - {n for n,_ in LIBRARY_ORDER}`.
`--libraries blss` produces `libraries=[]`, then the rest of `cmd_stage`
still copies cmake/helpers/tests and writes `target_config.cmake` with
`STAGED_LIBRARIES=''`.

**Fix:** validate the input set; raise on unknown names with a list of
valid names.

### C15 — Subset `--libraries` overwrites STAGED_LIBRARIES — `src/migrator/__main__.py:1311` ✓

`staged_list = ';'.join(staged)` only reflects libraries iterated in
**this** run. Line 1371 unconditionally rewrites `target_config.cmake`,
clobbering the previous run's listing. The just-removed libraries are
still on disk (manifest.cmake + src/) but excluded from the unified
build by the new STAGED_LIBRARIES.

**Fix:** merge `staged` with the previous file's STAGED_LIBRARIES list
before writing, or document re-stage as full-stage-only.

### C16 — `file(REMOVE …)` mutates source tree at configure time — `cmake/CMakeLists.txt:1131,1134` △

```cmake
file(REMOVE "${CMAKE_CURRENT_SOURCE_DIR}/_mpiseq_src/mpif.h")
file(REMOVE "${CMAKE_CURRENT_SOURCE_DIR}/_mpiseq_src/mpi.h")
```

Deletes libseq's bundled headers when MPI is found, mutating the staging
directory. If the user later reconfigures with MPI disabled (or wipes
the build dir alone), the headers are gone and the "use libseq's bundled
copies" fallback is unreachable.

**Fix:** `configure_file` into `CMAKE_CURRENT_BINARY_DIR` and put that
on the include path conditionally, instead of deleting from the source
tree.

### C17 — Error message names a non-existent `name_overrides:` YAML key — `src/migrator/prefix_classifier.py:122,159` ✓

The cross-family collision error tells the user to add a
`name_overrides:` entry in `targets/<target>.yaml` (with
`targets/multifloats.yaml` as an example). Grep across the repo: the
string appears only in the comment and the error message itself.
`target_mode.py` has no `name_overrides` field; no YAML reader looks
for it.

**Fix:** implement the field, or rewrite the error to reference the
actual mechanism (target prefixes + recipe `extra_renames`).

### C18 — `tests/CMakeLists.txt` LIB_PREFIX gate collapses when empty — `tests/CMakeLists.txt:47-50` △

```cmake
if(TARGET ${LIB_PREFIX}${_lib} AND EXISTS ${_libdir}/CMakeLists.txt)
```

If `LIB_PREFIX` is empty for any reason outside `BASELINE_BUILD` (recipe
load failure, etc.), the gate degenerates to `TARGET <lib>` — the std
archive accidentally satisfies it, and `tests/<lib>/CMakeLists.txt` runs
its non-baseline branch against the std archive.

**Fix:** add an explicit `if(NOT BASELINE_BUILD AND LIB_PREFIX STREQUAL
"") message(FATAL_ERROR "…")` guard before the loop.

---

## SIMPLIFICATION (aggressive, no behavior change)

### S1 — Delete duplicate `DetectExtendedPrecision.cmake` ✓

Root-level `DetectExtendedPrecision.cmake` is byte-identical (md5
`8f40c09817e45f25914d7f03b9dd3618`) to `cmake/DetectExtendedPrecision.cmake`.
All consumers reach the `cmake/` copy:
- `cmake/CMakeLists.txt:70` → `cmake/DetectExtendedPrecision.cmake`
- `src/migrator/__main__.py:955,1378,1554` → `cmake/DetectExtendedPrecision.cmake`

**Saving:** one file (51 lines), zero risk.

### S2 — `INTRINSIC_DECL_MAP` is derivable in one line — `src/migrator/intrinsics.py:265-357` △

```python
INTRINSIC_DECL_MAP = {k: v[0] for k, v in INTRINSIC_MAP.items() if k != v[0]}
```

Every key in `INTRINSIC_DECL_MAP` exists in `INTRINSIC_MAP` with matching
value (modulo the tuple flatten); the only divergence is the CMPLX self-
rename which the comprehension filters.

**Saving:** ~90 lines → 1 line.

### S3 — Flang variants of `migrate_fixed_form` / `migrate_free_form` are ~95% duplicates — `src/migrator/fortran_migrator.py:3987-4124` △

`_migrate_fixed_form_flang` (L3987-4062) vs `migrate_fixed_form`
(L3042-3110). Same 4-call pipeline (`_dedup_intrinsic_stmts` ×2 +
per-stmt transforms + `specialize_use_module`). Flang variants only gate
`has_float_types` / `has_real_literals` checks. Similar for
`_migrate_free_form_flang` (L4065-4124) vs `migrate_free_form`
(L3238-3347).

**Saving:** ~200 lines if parameterized with default-True flags.

### S4 — mpiseq Fortran stubs are a 25-line generator — `mpiseq_qx_stubs.f` + `mpiseq_mw_stubs.f90` △

30 near-identical stubs: 5 routines (PQGETRF/PQGETRS/PQPOTRF/PQPOTRS/
PQTRTRS) × {Q,X,E,Y} = 20 in qx_stubs.f + 5 × {M,W} = 10 in mw_stubs.f90.
Every body is `WRITE(*,*) '… should not be called.' / STOP`. Type varies
only by precision.

**Saving:** ~360 lines → ~25-line generator (Python script or fypp
template). Also removes the qx/mw split which is artificial (kind10 and
kind16 are independent precisions but bundled in qx_stubs).

### S5 — `convert_parameter_stmts` has two near-identical branches — `src/migrator/fortran_migrator.py:1838-2101` △

Combined-form (`TYPE, PARAMETER :: name = val`, L1918-2036) and
standalone-form (`PARAMETER(name = val)`, L2039-2101) independently run
the same 6-clause `is_cx_value` check, the same drop-known-constant
logic, the same convert-to-assignment, and the same cpp_stack wrap.

**Saving:** ~80 lines via `_entry_to_assignment(name, val, …) ` helper.

### S6 — `replace_intrinsic_decls` is redundant — `src/migrator/fortran_migrator.py:964-1107` △

`_dedup_intrinsic_stmts` already runs both before AND after the per-line
loop in each `migrate_*` function (e.g. L3052 + L3108). The per-line
`replace_intrinsic_decls` (L964-1012) reimplements the single-line case
of `_dedup_intrinsic_stmts`.

**Saving:** ~50 lines — either inline into `_dedup_intrinsic_stmts`'s
single-line branch, or drop the per-line call from L3084, L3313, L4036,
L4090.

### S7 — Fixed-form continuation detection open-coded in 4 places — `src/migrator/fortran_migrator.py:1036, 1055, 3018, 3487` △

Canonical `is_continuation_line` at L2880 already exists. The four
ad-hoc copies have subtle differences; L3487 even omits the cols 1-5
blank-check, opening a false-positive risk.

**Saving:** ~15 lines + L3487 risk closed.

### S8 — Per-library wiring blocks repeat 4-step pattern × 8 libs — `cmake/CMakeLists.txt:561-895` △

8 sections (blas L561-605, xblas L650-698, lapack L700-716, ptzblas
L743-757, pbblas L766-799, pblas L810-836, scalapack L845-895,
scalapack_c) all do: `add_standard_*` (optional) → `if(NOT BASELINE_BUILD)
add_migrated_*` → foreach deps → link std.

**Saving:** ~150 lines via `_wire_migrated_lib(name DEPS … [STD_TARGET
<bare>])` helper.

### S9 — Test-helper .f90s are byte-identical across 5-6 test trees — `tests/<lib>/common/*` △

Per md5sum (agent-confirmed): `prec_kinds.f90` identical across
blas/lapack/pblas/pbblas/scalapack/blacs (all 6); `compare.f90` across 5
(lapack diverges with extra sort helper); `test_data.f90` across 5
(lapack diverges with SPD/HPD/banded helpers); `pblas_distrib.f90`
across pbblas+scalapack.

Each `tests/<lib>/CMakeLists.txt` compiles its own copy into
`<lib>_test_common`. Per-tree fmod subdirs exist solely to keep the
identical .mod files from colliding (see e.g.
`tests/lapack/CMakeLists.txt:54-71`).

**Saving:** ~5× redundant compile, deletes per-tree fmod workarounds.
Lift truly-shared helpers into `tests/common/` alongside `target_conv`.

### S10 — Two recipe-YAML keys are declared-but-unread — `src/migrator/config.py:23` ✓

`_KNOWN_RECIPE_KEYS` advertises `symbols.method` and `prefix` but:
- `config.py` only reads `symbols_cfg.get('library_path')`; ignores
  `symbols.method`.
- `RecipeConfig` has no `prefix` field.

Every recipe ships `symbols: method: scan_source` (11 recipes) and
`prefix: style: …` (12 recipes) — pure config noise.

**Saving:** drop from `_KNOWN_RECIPE_KEYS` + remove from every
`recipes/*.yaml`, or implement the documented behavior.

### S11 — Five ParseTreeFacts fields are write-only — `src/migrator/flang_parser.py:213` + `gfortran_parser.py` △

`data_stmts`, `parameter_stmts`, `char_literals`, `procedure_boundaries`,
`intrinsic_names` — grep across `src/` shows two writers (in each parser)
and **zero readers**. Both parsers redundantly populate them; the
`DataInfo`/`ParameterInfo`/`ProcInfo` dataclasses backing them are
unused.

**Saving:** drop fields + writers + backing dataclasses.

### S12 — Redundant `find_package(MPI)` in generated Config — `cmake/FortranCompiler.cmake:357-359` △

```cmake
find_dependency(MPI COMPONENTS C Fortran)
find_package(MPI QUIET COMPONENTS C)
```

`find_dependency` already populates `MPI_C_FOUND`, `MPI_C_HEADER_DIR`,
`MPI::MPI_C`. The second call is a no-op.

**Saving:** delete L359.

---

## REFUTED / DEMOTED

- **Agent finding "`_helpers_default` fallback dead"** (`__main__.py:713`) —
  finding stands as written but `cmd_stage` explicitly copies helpers at
  L1313-1320 regardless of lib_name, so the broken fallback only fires
  in `cmd_build`'s single-recipe path. Lower impact.
- **Agent finding "case-insensitive FS double-count"** (`__main__.py:390`,
  `glob('*.f') + glob('*.F')`) — only matters on macOS HFS+; on Linux
  (the project's CI target) these are distinct files. Low priority.
- **Agent finding "compile_blacs.sh `set -e` missing"** — `set -uo
  pipefail` is present; partial protection. Tier-2 not Tier-1; rolled
  into the broader compile-script audit (#C10).

---

## Recommended fix order

The cheapest highest-impact path:

1. **#C1** — change `MPI_LONG_DOUBLE`→16, `MPI_C_LONG_DOUBLE_COMPLEX`→32
   in `mpiseq_c_stubs.c`. One file, two literals, prevents silent kind10
   data corruption.
2. **#S1** — delete duplicate root `DetectExtendedPrecision.cmake`. Zero
   risk.
3. **#C4** — add `_add_target_args(p)` to `verify` subparser. One line.
4. **#C9** — drop `\[` from `_BARE_LITERAL_RE` lookbehind. One char.
5. **#C11** — fix `_FORTRAN_DEF_RE` to handle multi-space `END  SUB`.
6. **#C12** — extend `_USE_STMT_RE` for F2003 `USE ::` form.
7. **#S10** — drop unused `symbols.method` / `prefix.style` from
   `_KNOWN_RECIPE_KEYS` + the YAML files (or implement them).
8. **#S2** — collapse `INTRINSIC_DECL_MAP` to comprehension.
9. **#C5** — std-dir patching asymmetry — needs design decision (apply
   patches in `cmd_stage` vs document asymmetry).
10. **#C2, #C3** — namespace / FetchContent export issues — verify with
    a clean install + downstream `find_package` smoke test.

The remaining Tier-1 correctness items (#C6, #C7, #C8, #C10) and the
larger simplifications (#S3, #S4, #S5, #S6, #S8, #S9) are independent
and can be sequenced freely.

## Verification commands (reproduce findings)

```bash
# C1 — confirm long-double size
echo '#include <stdio.h>
int main(){printf("%zu\n", sizeof(long double));}' | gcc -xc - -o /tmp/sz && /tmp/sz   # → 16

# C9 — array constructor literal drop
uv run python -c "from migrator.prepare import _BARE_LITERAL_RE; \
  print(_BARE_LITERAL_RE.findall('[5.0, 6.0]'))"   # → ['6.0']

# C11 — END double-space
uv run python -c "from migrator.symbol_scanner import _FORTRAN_DEF_RE; \
  print(_FORTRAN_DEF_RE.findall('END  SUBROUTINE FOO'))"   # → ['FOO']

# C12 — F2003 USE forms
uv run python -c "from migrator.gfortran_parser import _USE_STMT_RE; \
  print(_USE_STMT_RE.match('    USE :: foo'))"   # → None

# C17 — name_overrides nonexistent
grep -rn 'name_overrides' src/ recipes/ targets/   # → only the error msg

# S1 — confirm cmake duplicate
md5sum DetectExtendedPrecision.cmake cmake/DetectExtendedPrecision.cmake
```

## Project norms (preserve)

- **No commits without explicit user authorization.** Earlier this session
  the rename commits were each authorized one-by-one; the residue
  cleanup was applied to working tree, not committed. Don't stage or
  commit any of these fixes unless explicitly asked.
- **Auto Mode is on** — bias toward action, no unnecessary clarifying
  questions; but irreversible actions (public repo creation, pushes,
  force-pushes, destructive git ops) still need explicit confirmation.
- **`migrator` is the Python CLI module name** and stays unchanged — not
  packaging residue, not subject to the eplinalg rename.
- **Per-precision package names** (qblas, eblas, mblas, etc.) stay —
  not residue.

## Suggested skills

For the next agent picking this up:

- **`grill-with-docs`** — before fixing #C5 (std-dir patching asymmetry)
  or #C16 (`file(REMOVE)`), walk the design through against CONTEXT.md
  and any ADRs in `doc/adr/` — these are architectural calls, not
  mechanical edits.
- **`code-review`** with `--fix` (or just `simplify`) — after a batch
  of mechanical fixes, run on the working diff to catch regression-style
  mistakes before committing.
- **`verify`** — for #C2 / #C3 namespace+export changes, drive a clean
  `migrator stage` + `cmake --install` + downstream `find_package`
  smoke test to confirm the fix.
- **`run`** — for #C1 (mpiseq byte sizes), exercise a kind10 MUMPS test
  before and after to confirm the corruption it was producing.

## Reference paths

- Conversation transcript:
  `/home/kyungminlee/.claude/projects/-home-kyungminlee-code-eplinalg/f565fd74-cef4-4be0-9682-6cf8ccede299.jsonl`
- Subagent raw outputs:
  `/tmp/claude-1000/-home-kyungminlee-code-eplinalg/f565fd74-cef4-4be0-9682-6cf8ccede299/tasks/{ae79d82ca61dd1305,a1c318a90f1fafe92,a6c9a5ae4c7806fea,acfefbf19409cc20b,ab58714c99cfbe02d}.output`
- Earlier rename-residue handoff (already applied, now empty):
  `/tmp/handoff-pre-split-residue-20260528.md`
