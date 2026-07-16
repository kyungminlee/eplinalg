# tests/mumps — CHANGELOG

Resolved items, reverse-chronological. Open work lives in `TODO.md`.

## 2026-05-05 — libmpiseq linkage extended to multifloats

Multifloats `_seq` tests now pass 26/26, with bit-identical JSON
precision reports vs the impi-linked runs (matches the kind16/kind10
result already on develop). The blocker — `multifloats_mpi.cpp`
registering its custom datatypes via `MPI_Type_create_struct` /
`MPI_Op_create`, which the previous libmpiseq stubs returned as
`MPI_DATATYPE_NULL` / `MPI_OP_NULL` — is closed by giving libmpiseq's
derived-type stubs an Intel-disjoint sentinel ABI:

- `cmake/mpiseq_c_stubs.c`: `MPI_Type_contiguous` / `MPI_Type_vector` /
  `MPI_Type_indexed` / `MPI_Type_create_struct` now compute the total
  byte-size from the input arguments (via a small `mpiseq_base_type_bytes`
  switch on Intel's named-type constants) and return
  `0x10000000 | nbytes` as the new handle. The 0x10000000 tag is below
  `MPI_OP_NULL` (0x18000000) and well clear of Intel's 0x4c00****
  datatype range, so the sentinel is unambiguous on either side of the
  Fortran/C bridge. `MPI_Op_create` returns
  `0x18000000 + counter` — distinct non-null handles above
  `MPI_OP_NULL`, since the user-op callback never fires under
  single-rank ALLREDUCE (it collapses to memcpy in `MUMPS_COPY`).
- `src/migrator/__main__.py` (`_patch_libseq_mpi_f`): adds explicit
  `MUMPS_COPY` dispatch cases for the two multifloats sentinels
  (`268435472` → 16-byte real64x2 elements via `MUMPS_COPY_FLOAT64X2`;
  `268435488` → 32-byte complex64x2 elements via
  `MUMPS_COPY_COMPLEX64X2`). The new helpers are byte-equivalent block
  copies through `REAL(KIND=8)` arrays of the right multiplicity.

`MPI_Type_c2f` and `MPI_Op_c2f` are passthrough casts in Intel mpi.h
(`#define MPI_Type_c2f(d) (MPI_Fint)(d)`), so the Fortran handle that
multifloats stores in `mf_mpi_float64x2_f` etc. is the same sentinel
the Fortran-side `MUMPS_COPY` checks against.

Verified: kind16 / kind10 / multifloats all 26/26 on both impi and
`_seq` variants; `diff -q` over per-test JSON reports yields zero
differences on every stage.

## 2026-05-05 — End-to-end libmpiseq linkage on kind16 + kind10

Each `tests/mumps/{fortran,c}/test_*.{f90,c}` source now produces a
second executable suffixed `_seq` that links the in-tree `mpiseq`
archive in place of `MPI::*` / `mblacs` / `mscalapack` / `mpblas` /
`mpbblas` / `mptzblas`. Same compiled archives, two link variants —
exercises the libmpiseq path downstream consumers use when they want
a no-MPI executable. The `_seq` ctests run as plain binaries (no
mpiexec), about ~150× faster per test than the impi-wrapped variant.

**kind16 + kind10 result**: 26/26 `_seq` tests pass on each target;
per-test JSON precision reports are bit-identical to the impi-linked
runs (`md5sum` confirms all 26 files match on both stages).

Implementation breakdown:

- **Header overlay** (cmake/CMakeLists.txt) — when `MPI_*_FOUND`,
  replace `_mpiseq_src/{mpi.h,mpif.h}` with the configured MPI's copies
  via `file(REMOVE)` + `target_include_directories(mpiseq BEFORE ...)`,
  so libmpiseq's stubs see the **same constants and COMMON-block
  layout** (`MPI_DOUBLE_PRECISION = 1275070495`, `MPI_IN_PLACE` in
  `/MPIPRIV1/`, ...) as the migrated archives that consume them.
- **mpic.c replacement** (cmake/mpiseq_c_stubs.c) — libseq's bundled
  `mpic.c` declares `MPI_Comm_f2c` as a function but Intel's `mpi.h`
  defines it as a function-like macro; the two collide before any
  compile-flag tricks. Folded the four lifecycle functions
  (`MPI_Init` / `MPI_Comm_rank` / `MPI_Finalize` / `MPI_Comm_f2c` /
  `MPI_Wtime` + `MUMPS_CHECKADDREQUAL` aliases) into our supplementary
  C stubs file; libseq's `mpic.c` is dropped from the mpiseq target
  source list. C signatures match Intel's exactly (`const` qualifiers,
  no `LIBSEQ_CALL` macro, opaque types from Intel's `mpi.h`).
- **MUMPS_COPY extension** (src/migrator/__main__.py:_patch_libseq_mpi_f) —
  migrator stage patches the staged copy of `_mpiseq_src/mpi.f` to add
  `MPI_REAL16` / `MPI_COMPLEX32` (kind16) and `MPI_LONG_DOUBLE` /
  `MPI_C_LONG_DOUBLE_COMPLEX` (kind10 — there's no `MPI_REAL10` in
  standard MPI, so the migrator emits the long-double tokens) cases
  plus matching `MUMPS_COPY_REAL{10,16}` / `MUMPS_COPY_COMPLEX{20,32}`
  helpers, so reductions on extended-precision buffers dispatch
  correctly. Upstream `external/MUMPS_5.8.2/libseq/mpi.f` stays
  read-only.
- **PUBLIC → PRIVATE on MPI deps** — `multifloats_mpi`,
  `${LIB_PREFIX}mumps`, `mumps_common`, `mumps_c_bridge` flipped from
  `PUBLIC MPI::*` to `PRIVATE MPI::*`. The `_seq` link path then
  doesn't get Intel MPI dragged back in via STATIC-archive PRIVATE-dep
  propagation; the impi test exes were already adding `MPI::*`
  explicitly via `_mumps_apply_link_variant`.
- **Manual link group** (tests/mumps/CMakeLists.txt) — the `_seq`
  archive list is hand-rolled inside a single `-Wl,--start-group`
  block via `$<TARGET_FILE:>`, bypassing cmake's link-graph cycle
  detector (mmumps PUBLICly chains scalapack / blacs / pblas etc.
  that we deliberately replace with libmpiseq).

**Multifloats**: `_seq` binaries currently STOP at
`MPI_ALLREDUCE, DATATYPE = 201326592` (Intel's `MPI_DATATYPE_NULL`)
because the C++ `multifloats_mpi.cpp` registers `MPI_FLOAT64X2` /
`MPI_DD_SUM` etc. via `MPI_Type_create_struct` / `MPI_Op_create` and
our stubs return `MPI_DATATYPE_NULL`. Multifloats-specific follow-up,
documented in TODO.md.

## 2026-05-05 — Multifloats (M/W) libmpiseq stubs

`cmake/mpiseq_mw_stubs.f90` adds the 10 multifloats-precision ScaLAPACK
stub forwarders (`p[mw]{getrf,getrs,potrf,potrs,trtrs}_`) symmetric to
the existing Q/X/E/Y set in `mpiseq_qx_stubs.f`. Free-form Fortran so
the multifloats module's `TYPE(real64x2)` / `TYPE(cmplx64x2)` can be
USEd for the array-argument shadows. Wired into the `mpiseq` target in
`cmake/CMakeLists.txt` only when `NEEDS_MULTIFLOATS AND TARGET
multifloatsf` (so plain kind10/kind16 builds don't pull multifloats).

Verified the migrated `libmmumps-gfortran-13.a` references exactly
`p[mw]{getrf,getrs,potrf,potrs}_` (TRTRS goes through the sequential
`mtrtrs_`/`wtrtrs_` from lapack), all of which the new stubs export.
The full sequential-link path (mmumps relinking against `mpiseq`
instead of `MPI::MPI_Fortran`) is still not exercised by any test —
this lifts the symbol gap so it can be without further migrator work.

## 2026-05-04 — B4: full multi-target support

All three targets (kind10 / kind16 / multifloats) build, link, and run
the full mumps test set under `cmake --preset=linux-impi`: 26/26 mumps
ctests pass on multifloats; the kind10 + kind16 runs match. Recipe
`overrides:` for all three targets in `recipes/mumps.yaml` ship the
matching `mumps_memory_mod_ep.F` / `mumps_lr_stats_ep.F`; the
multifloats variant uses `STORAGE_SIZE`-driven byte accounting so the
double-double layout (16-byte real, 32-byte complex) is computed at
compile time rather than hardcoded.

Verified `INFOG(20) = 1024 = n²` exactly on every target for the fixed
(n=32, SYM=0, JOB=6, seed=4099) problem in
`tests/mumps/fortran/test_dmumps_infog20.f90`. The test was tightened
from the structural `[n², 50·n²]` window to ±5% against the captured
`infog20_baseline = 1024` — drift that would fly under the loose bound
(off by an EP-record-size factor) is now flagged immediately.

## 2026-05-03 — B9c: multifloats MUMPS `MPI_SUM` rewriter (commit `eb83d5a`)

A token-context MPI reduction-op rewriter was added to the Fortran
migrator at `src/migrator/fortran_migrator.py:3291-3320`
(`_rewrite_mpi_sum`), wired into the post-migration pipeline
immediately after `_rewrite_mpi_datatypes`. It matches each
`MPI_(I?All)Reduce[_Scatter]` call as a unit and only swaps
`MPI_SUM/MPI_MAX/MPI_MIN → MPI_DD_*` (real) or `→ MPI_ZZ_*` (complex)
inside calls whose argument list contains the rewritten floating
datatype token (`MPI_FLOAT64X2` / `MPI_COMPLEX64X2`). Integer
reductions in the same translation unit are left alone.

Original symptom: after Group A (B9b) unblocked the multifloats mumps
test harness, runs failed with `internal_Allreduce: MPI_Op operation
not defined for this datatype` (USER<contig>, MPI_SUM). Root cause:
`_rewrite_mpi_datatypes` rewrote the buffer-type argument
(`MPI_DOUBLE_PRECISION → MPI_FLOAT64X2`) but not the matching
reduction-op argument. Sites originally observed at
`mfac_scalings.F:372,376`, `mana_dist_m.F:206,605,626`, and
`[mw]sol_distrhs.F`. The fix addresses every site uniformly without
per-file recipe overrides.

## 2026-05-02 — Group A / B9b: multifloats test harness

C-side and Fortran tests are now target-aware across all three
precision targets:

- **C-side tests (kind10, kind16)** — `tests/mumps/c/include/test_real_compat.h`
  picks `__float128` for kind16 and `long double` for kind10 along with
  matching literal suffix (`q` / `L`), math (`fabsq` / `fabsl`,
  `sqrtq` / `sqrtl`), epsilon, and snprintf wrapper. CMake glue
  (`add_mumps_c_test()`) wires `TARGET_REAL_HEADER` /
  `TARGET_REAL_MUMPS_C` / `TARGET_REAL_STRUC_C` (plus complex
  counterparts) macros.
- **C-side tests (multifloats)** — `test_real_compat.h` gained a
  `TEST_TARGET_MULTIFLOATS` branch: `test_real` is plain `double` (host
  precision), `test_complex` is `{double r, i}`, and
  `tr_widen` / `tr_narrow` (real) / `tc_widen` / `tc_narrow` (complex)
  translate at the bridge boundary to `mumps_float64x2` /
  `mumps_complex64x2`. The 4×4 hand-built test problems don't strain
  double precision; the C tests verify the bridge wires correctly, not
  multifloats arithmetic in C.
- **Multifloats Fortran tests** — drivers now
  `use target_mumps, only: ..., q2t_r, t2q_r` (or `_c` for the complex
  side) and wrap `id%A` / `id%RHS` / `id%A_loc` / `id%RHS_*`
  assignments with `q2t_r(...)` and readbacks with `t2q_r(...)`. The
  `target_conv` library from `tests/common` provides identity
  converters on kind16, intrinsic narrow + widen on kind10, and
  double-double split on multifloats. `target_mumps_body.fypp`
  re-exports the converters through host association.
- **Parity tests** — `test_dmumps_c_parity.f90` /
  `test_zmumps_c_parity.f90` rewritten as `.fypp` and run on all three
  targets. The bind(C) interface declares the per-target type
  (`real(16)` / `real(10)` / `type(real64x2)`); a per-target staging
  block converts via `q2t_r` / `t2q_r` before/after the C call. CMake
  processes them with `-Dtarget=${TARGET_NAME}`.
- **Skip-gate removal** — the `_MUMPS_TESTS_SKIP_EXECUTABLES` skip guard
  at the head of `tests/mumps/CMakeLists.txt`, the kind16-only
  `_c_parity\.f90$` filter, and the
  `NOT TARGET_NAME STREQUAL "multifloats"` clause on the C-test loop
  were removed. All three targets register the same Fortran + C-side
  test sets.

## 2026-05-02 — UNK-3 / UNK-4 verification

- `NO_SAVE_RESTORE` macro is **undefined** (`recipes/mumps.yaml` does
  not define it; build does not pass `-D`), so the SAVE/RESTORE branches
  are compiled in. The F77 prototype the C bridge calls
  (`${LIB_PREFIX}mumps_f77_`) is the SAVE/RESTORE-extended one,
  matching upstream's default.
- Zero unresolved PORD symbols in the migrated archive
  (`nm libqmumps-gfortran-13.a | grep -i pord | grep ' U '` → 0). PORD
  is reached through the C-side `mumps_pord.c` (compiled in
  `mumps_c_bridge`), not Fortran.

## 2026-05-02 — D1: invalid-input wrapper

`target_mumps` exports `check_dmumps_input` / `check_zmumps_input` plus
`MIC_*` codes (`tests/mumps/common/target_mumps_body.fypp`). Restored
tests at `tests/mumps/fortran/test_{d,z}mumps_errors.f90` exercise
BAD_N (neg + zero), BAD_NNZ, BAD_IRN (high + zero), BAD_JCN (high +
zero), SIZE_MISMATCH, plus a final valid-input pass that reaches
MIC_OK and factors via JOB=6.

Background: MUMPS 5.8.2 doesn't validate `id%N = -1`,
`id%IRN(1) = N + 5`, etc. — those crash with SIGSEGV inside analysis
or factorization rather than returning `INFOG(1) < 0` per the manual.
The wrapper layers above MUMPS without touching `id%INFOG(1)`.

## 2026-05-02 — G2: per-target wrapper directories

Both wrapper directories ship on disk:
- `tests/mumps/target_kind10/target_mumps.fypp`
- `tests/mumps/target_multifloats/target_mumps.fypp`

Each ships a precision-aware `c/include/mumps_c_types.h` shadow plus
the per-target `*mumps_c.h` wrappers.

## 2026-05-01 — B8 series: multifloats migrator pipeline fixes

Multiple latent migrator bugs blocked the multifloats `mmumps` build.
Resolved (all in commit `bf33cde` + follow-ups):

- **B8** — `USE multifloats` injector tracked parenthesis depth
  through CPP-interspersed argument continuation; helper
  `_count_open_parens` returns net paren delta per line, ignoring
  quoted strings and inline comments. CPP `#if/#endif/#else/#elif/
  #define` directives now treated as transparent.
- **B8b** — `int(real64x2_expr, KIND_INT)` had no multifloats
  overload. `_rewrite_int_kind_on_real64x2` now wraps the inner with
  `dd_to_double(...)` (later replaced with `dble`/`real` per B8f) so
  the standard Fortran `INT` intrinsic handles the kind selector.
- **B8c** — `MODULE PROCEDURE` lines inside INTERFACE blocks were
  matched as procedure-header starts. Negative lookahead
  `(?!\s+PROCEDURE\b)` after `MODULE` in `_PROC_HEADER_RE` and the
  local `proc_header_re` inside `insert_use_multifloats`.
- **B8d** — `dble(...)` masked by keep-kind sentinel was missing from
  USE only-list. `_build_use_only_clause` now detects
  `_KK_DBLE_SENTINEL` in the proc body text and unconditionally adds
  `dble` to the referenced set so multifloats's
  `dd_dble(real64x2) -> real(dp)` overload is imported.
- **B8e** — `MPI_FLOAT64X2` was undeclared in migrated MUMPS Fortran.
  The `insert_use_multifloats_mpi_f` post-pass walks each procedure
  header (with proper depth tracking through INTERFACE-block inner
  SUBROUTINEs) and emits `USE multifloats_mpi_f` whenever the body
  references one of the handle names. Scope: real
  `MPI_FLOAT64X2` + complex `MPI_COMPLEX64X2` + reduction ops
  (`MPI_DD_*`, `MPI_ZZ_*`).
- **B8f** — Pre-existing migrator bugs surfaced by B8e:
  `insert_use_multifloats`'s decl-block walker missing `INCLUDE` and
  keep-kind sentinels; `convert_parameter_stmts` losing `#if defined`
  guards; `_rewrite_int_kind_on_real64x2` / `_rewrite_int_of_complex`
  emitting non-existent multifloats symbols (now emit `dble`/`real`);
  `end_proc_re` falsely matching `END IF`/`END DO`.
- **B8g** — Migrator wrapped PARAMETER RHS literal even when LHS type
  was keep-kind DOUBLE PRECISION. `replace_literals` now early-returns
  unchanged in `constructor` mode whenever any keep-kind sentinel
  appears on the line.
- **B8h** — More pre-existing bugs surfaced by B8g: combined-form
  `TYPE, PARAMETER ::` declarations now split into plain decl + runtime
  assignment; file-scope COMPLEX/REAL ambiguity for known-constant
  names handled via per-file `complex_names` set; renames overwriting
  LHS of converted-PARAMETER assignments now suppressed; keep-kind
  sentinel masking declaration detection fixed; `int(var, KIND_INT)`
  over real64x2 variables now triggers; INTERFACE blocks now treated
  as transparent for runtime-assignment placement.

## 2026-05-01 — B9: per-target C bridge

The bridge cmake glue (`tests/mumps/CMakeLists.txt`) and the
per-arithmetic header shadow
(`tests/mumps/target_${TARGET_NAME}/c/include/`) follow `${LIB_PREFIX}`
/ `${LIB_PREFIX_COMPLEX}`. Each target ships its own:

- `mumps_c_types.h` — precision-aware `DMUMPS_REAL` / `ZMUMPS_COMPLEX`
  typedefs (kind16: `__float128` / two-`__float128` complex; kind10:
  `long double` / two-`long double` complex; multifloats:
  `mumps_float64x2` POD / two such with `.r/.i`).
- `${LIB_PREFIX}mumps_c.h` / `${LIB_PREFIX_COMPLEX}mumps_c.h` — thin
  wrappers around upstream `dmumps_c.h` / `zmumps_c.h` with per-target
  macro renames.

Bridge OBJECT-library names are now
`${LIB_PREFIX}mumps_c_bridge` / `${LIB_PREFIX_COMPLEX}mumps_c_bridge`,
with `-Ddmumps_c=${LIB_PREFIX}mumps_c` etc. driving the per-arithmetic
dispatch. `c_parity_helpers.c` is templated through the C macros
`TARGET_REAL_HEADER` / `TARGET_REAL_MUMPS_C` / `TARGET_REAL_STRUC_C`
(plus complex counterparts) so it builds unchanged against any target.
`migrator stage` writes `LIB_PREFIX_COMPLEX` to `target_config.cmake`
alongside the existing `LIB_PREFIX`.

## 2026-05-01 — B7: cmake/recipe tree consolidation

The fm-mumps split is gone. The mumps work that briefly lived in a
sibling tree was merged into fortran-migrator on the `tests` branch
(commits `9ec8357`, `0b92df6`). All three originally-diverged items
sit in one place:

- `cmake/CMakeLists.txt` — section 9 (MUMPS) at line 747; section 10
  (libmpiseq) at line 779.
- `recipes/mumps.yaml` — `FAC_FUTURE_NIV2_MOD` wired at line 89.
- `src/migrator/__main__.py` — `('mumps', 'mumps.yaml')` in
  `LIBRARY_ORDER` at line 724.

`migrator stage` reads from this repo with no `--project-root`
override needed.

## 2026-04-30 — B3 partial: libmpiseq stub buildable; Q/X path-b stubs (2026-05-02)

`libmpiseq.a` is buildable. `external/MUMPS_5.8.2/libseq/` copied to
`_mpiseq_src/` by `migrator stage` (added to `_std_dirs` in
`src/migrator/__main__.py`); `cmake/CMakeLists.txt` (section 10)
declares `add_library(mpiseq STATIC mpi.f mpic.c elapse.c)` with the
required `-DAdd_` flag. Exports 157 symbols.

Supplementary C-MPI stubs live in `tests/mumps/c/mpiseq_c_stubs.c`
(~30 additional `MPI_*` functions referenced by standard BLACS / PBLAS
C archives). Compiled separately and linked alongside libmpiseq.

**Full sequential link is non-trivial.** The migrated qmumps archive
depends on Q-prefixed Fortran routines (`pqgetrs_`, `pqgetrf_`,
`pqpotrf_`, `pqpotrs_`) that libmpiseq provides only as upstream
D-prefixed equivalents.

**2026-05-02 follow-up (path b, partial):** Q/X/E/Y per-precision
stubs landed in `cmake/mpiseq_qx_stubs.f` and are appended to the
`mpiseq` target via `cmake/CMakeLists.txt:835`. Twenty new stubs
exported (`p[qxey]{getrf,getrs,potrf,potrs,trtrs}_`), each cloning
upstream `mpi.f`'s "should not be called / STOP" body. Multifloats
(M/W) prefixes intentionally NOT covered yet.

Linking a test executable against `mpiseq` (instead of
`MPI::MPI_Fortran`) plus `libqscalapack` / `libqblacs` should now
resolve all symbols for kind10 + kind16, but no test in this repo
exercises that path today (mpiexec -n 1 stays the default).

## 2026-04-30 — G1: Fortran/C cross-language parity test

`fortran/test_dmumps_c_parity.f90` and `test_zmumps_c_parity.f90` plus
`c/c_parity_helpers.c` (the Fortran-callable bridge into qmumps_c /
xmumps_c). Each test drives the same generated problem through both
paths and verifies bit-identical solutions. Catches silent
struct-extraction corruption that the per-side basic tests can't
(since they only check `result ≈ x_true` to ~33 digits).

Avoided the `__float128` ISO_C_BINDING gymnastics by using
`real(16)` / `complex(16)` in the Fortran interface block — gfortran
and gcc agree on the by-reference layout. Verified
`fortran-vs-c = 0.0e+00` on both D and Z paths.

## 2026-04-29 — B2: header-override C bridge

Recipe `skip_files` excludes every `*MUMPS_C` and `MUMPS_C_TYPES`
header, and `recipes/mumps.yaml` declares `language: fortran`. The C
side is provided in-tree by `tests/mumps/c/`:

- `tests/mumps/c/include/mumps_c_types.h` — precision-aware shadow of
  upstream's mumps_c_types.h (placed FIRST on the include path).
- `tests/mumps/c/include/mumps_int_def.h` — pin to `MUMPS_INTSIZE32`.
- `tests/mumps/c/include/${LIB_PREFIX}mumps_c.h` and
  `${LIB_PREFIX_COMPLEX}mumps_c.h` — wrappers around upstream
  `dmumps_c.h` / `zmumps_c.h` with macro renames.

Build glue at `tests/mumps/CMakeLists.txt`: type-agnostic C runtime
(upstream sources verbatim) + per-arithmetic `mumps_c.c` compiled twice
via `-D` macro renames. Renames preserve the `dmumps_assign_*` /
`dmumps_set_tmp_ptr_c_` C-side symbols that migrated Fortran calls
under the original D name, while exposing user-visible entry points
(`dmumps_c → ${LIB_PREFIX}mumps_c`) per target.

Verified end-to-end on 2026-04-29 by hand-linking a JOB=-1/-2
roundtrip.

## 2026-04-29 — B5: continuation-line inline comment loss

The migrator's `_strip_inline_comment()` helper, applied inside
`_segment_fixed_form_statements` when joining continuation lines,
removes `!` comments that previously survived into the joined logical
line and got mangled mid-statement after `reformat_fixed_line`'s
column-66 re-split. Inline comments are intentionally lost from the
joined / reformatted output (the price of correctness); single-line
statements with no transform keep their comments verbatim.

Original failure: `qlr_core.F:1379` had `! L Panel` swallowing
`K, BLR_L)` after continuation joining.

## 2026-04-29 — B6: `FAC_FUTURE_NIV2_MOD` recipe gap

The integer-only module `MUMPS_FUTURE_NIV2`
(`external/MUMPS_5.8.2/src/fac_future_niv2_mod.F`) was being classified
as PRECISION content but copy_files-resident `MUMPS_LOAD` (in
`mumps_common`) `USE`s it. mumps_common builds first → `mumps_load.F`
failed with `Cannot open module file 'mumps_future_niv2.mod'`.

Fixed by adding `FAC_FUTURE_NIV2_MOD` to `copy_files` in
`recipes/mumps.yaml`. (`copy_files` matches uppercased filename
basenames, not module names.)

## UNK-1 — Migrated symbol names

Migrator's `prefix: direct` rule **replaces** the leading precision
letter (D/Z → Q/X for kind16) for **filenames, modules, and
subroutines** but **leaves type names unchanged**.

| Upstream                  | Migrated kind16        |
| ------------------------- | ---------------------- |
| `DMUMPS` subroutine       | `QMUMPS`               |
| `DMUMPS_STRUC_DEF` module | `QMUMPS_STRUC_DEF`     |
| `dmumps_struc.h` header   | `qmumps_struc.h`       |
| `DMUMPS_STRUC` type       | `DMUMPS_STRUC` (unchanged) |
| `A`, `RHS`, `COLSCA`, `RINFO` fields | Promoted to `REAL(KIND=16)` / `COMPLEX(KIND=16)` |

The wrapping module is enough to disambiguate at use sites; cf.
`qmumps_driver.F:197` —
`TYPE (DMUMPS_STRUC), TARGET :: id`. Test code therefore writes:

```fortran
USE QMUMPS_STRUC_DEF                    ! migrated module
TYPE(DMUMPS_STRUC), TARGET :: id        ! type name unchanged
CALL QMUMPS(id)                         ! migrated subroutine
```

## UNK-2 — C bridge struct layout (N/A under header-override)

Original concern was that a hand-rolled `DMUMPS_STRUC_C` would have to
match the migrated Fortran derived type byte-for-byte. The
header-override approach (B2) sidesteps this entirely: the C struct is
declared by upstream `dmumps_c.h` itself, and is INDEPENDENT from the
Fortran `DMUMPS_STRUC` (upstream's `dmumps_c()` decouples them via
field-by-field extraction). No layout match required.
