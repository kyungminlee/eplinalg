# Fuzz / perf coverage survey — 2026-05-17

Follow-up to `doc/fuzz-stride-coverage-audit-20260518.md`. That audit
covered stride coverage only; this one expands the lens to all
categorical BLAS parameters (`trans`, `side`, `uplo`, `diag`) and
re-checks the stride gap.

Same risk shape as before: a code path that nothing exercises will
silently regress (correctness or perf). The loop-direction survey's
`ytrsv` Sub-class A bug is the canonical example.

## TL;DR

- **Stride gap from 2026-05-18 is still wide open.** Only kind10/kind16
  trsv/trmv (the shared `_body.fypp`) were fixed. 20 other L2 routine
  families still hardcode `incx=incy=1` in both fuzz and perf.
- **Multifloats fuzz was *not* included in the 2026-05-18 fix.**
  `target_multifloats/fuzz_m{trsv,trmv}.fypp` are standalone files
  (not body-include shims like kind10/kind16) and still pass `incx=1`.
- **DIAG='U' is never benched** for any triangular routine
  (trmv/trsv/tbmv/tbsv/tpmv/tpsv/trmm/trsm). The unit-diagonal path
  has its own arithmetic shape; a perf regression there is invisible.
- **GEMMTR perf samples only 3 of 9 (ta,tb) pairs.** TT, CN, NC, CT,
  TC, CC are unbenched. **GEMM complex** samples 6 of 9 (skips
  CT/TC/CC).
- **Real L2 triangular fuzz never passes `trans='C'`.** trmv/trsv/
  tbmv/tbsv/tpmv/tpsv real bodies use inline 50/50 N/T; the L3 real
  family (gemm/gemmtr/trmm/trsm) uses the common `rand_trans()`
  helper which emits N/T/C (with `'C'` semantically equal to `'T'`
  for real, but exercising any kernel branch on the literal).
- **Stride helper still fragmented.** `rand_incx`, `rand_inc`,
  `pick_inc`/`pick_inc_pos` coexist; the 2026-05-18 cleanup
  recommendation hasn't landed.

## Fuzz consistency — categorical coverage

Two code paths to audit:
- **Shared bodies** `consistency/fuzz_*_body.fypp` — used by kind10
  (`target_kind10/fuzz_e*.fypp`) and kind16 (`target_kind16/fuzz_q*.fypp`)
  via 2-line `#:include` shims.
- **Multifloats per-target** `consistency/target_multifloats/fuzz_m*.fypp`
  and `fuzz_iw*.fypp` / `fuzz_im*.fypp` — standalone programs, NOT
  shims (each ~65-100 lines with its own randomization logic).

Columns: ✅ randomized, ❌ hardcoded, — n/a.

### Shared bodies (kind10 + kind16 targets)

| Layer | Routine family | trans | side | uplo | diag | strides |
|---|---|---|---|---|---|---|
| L1 | axpy, copyswap, dot, rot, rotm, scal | — | — | — | — | ✅ (±1, ±2) |
| L1 | iamax, norm | — | — | — | — | ✅ pos only |
| L2 | gemv (R/C) | ✅ N/T (R), N/T/C (C) | — | — | — | ❌ |
| L2 | gbmv (R/C) | ✅ N/T (R), N/T/C (C) | — | — | — | ❌ |
| L2 | ger, gerc, geru | — | — | — | — | ❌ |
| L2 | symv (R), hemv | — | — | ✅ | — | ❌ |
| L2 | sbmv, hbmv | — | — | ✅ | — | ❌ |
| L2 | spmv, hpmv | — | — | ✅ | — | ❌ |
| L2 | syr, her, syr2, her2, spr, hpr, spr2, hpr2 | — | — | ✅ | — | ❌ |
| L2 | trmv, trsv (R/C) | ✅ N/T (R), N/T/C (C) | — | ✅ | ✅ | ✅ (fixed 05-18) |
| L2 | tbmv, tbsv, tpmv, tpsv (R/C) | ✅ N/T (R), N/T/C (C) | — | ✅ | ✅ | ❌ |
| L3 | gemm (R/C) | ✅ N/T/C (both R and C — `rand_trans()`) | — | — | — | n/a |
| L3 | gemmtr (R/C) | ✅ ta,tb N/T/C | — | ✅ | — | n/a |
| L3 | symm, hemm | — | ✅ | ✅ | — | n/a |
| L3 | syrk, syr2k (R) | ✅ N/T | — | ✅ | — | n/a |
| L3 | syrk, syr2k (C) | ✅ N/T (no 'C' — invalid for complex syrk per BLAS) | — | ✅ | — | n/a |
| L3 | herk, her2k | ✅ N/C | — | ✅ | — | n/a |
| L3 | trmm, trsm (R/C) | ✅ N/T/C (`rand_trans()`) | ✅ | ✅ | ✅ (30/70 U/N) | n/a |

### Multifloats per-target (separate code, not body-include)

Two parallel families — `fuzz_m*` (real DD, `real64x2`) and `fuzz_w*`
(complex DD, `cmplx64x2`). Both are standalone programs.

| Layer | Routine family | trans | side | uplo | diag | strides |
|---|---|---|---|---|---|---|
| L1 | maxpy, mcopy, mdot, mrot, mrotm, mscal, mswap | — | — | — | — | ✅ (±1, ±2) |
| L1 | masum, mnrm2, mwasum, mwnrm2, imamax, iwamax | — | — | — | — | ✅ pos only |
| L1 | mcabs1, mrotg, mrotmg | — | — | — | — | — |
| L1 | waxpy, wcopy, wdotc, wdotu, wmrot, wmscal, wrotg, wscal, wswap | — | — | — | — | ❌ (use legacy `pick_inc`) |
| L2 | mgemv | ✅ N/T | — | — | — | ❌ |
| L2 | mgbmv | ✅ N/T | — | — | — | ❌ |
| L2 | mger | — | — | — | — | ❌ |
| L2 | msymv | — | — | ✅ | — | ❌ |
| L2 | msbmv | — | — | ✅ | — | ❌ |
| L2 | mspmv | — | — | ✅ | — | ❌ |
| L2 | msyr, msyr2, mspr, mspr2 | — | — | ✅ | — | ❌ |
| L2 | mtrmv, mtrsv | ✅ N/T | — | ✅ | ✅ | ❌ **(missed by 05-18 fix)** |
| L2 | mtbmv, mtbsv, mtpmv, mtpsv | ✅ N/T | — | ✅ | ✅ | ❌ |
| L2 | wgemv | ✅ N/T/C | — | — | — | ❌ |
| L2 | wgbmv | ✅ N/T/C | — | — | — | ❌ |
| L2 | wgerc, wgeru | — | — | — | — | ❌ |
| L2 | whemv | — | — | ✅ | — | ❌ |
| L2 | whbmv | — | — | ✅ | — | ❌ |
| L2 | whpmv | — | — | ✅ | — | ❌ |
| L2 | wher, wher2, whpr, whpr2 | — | — | ✅ | — | ❌ |
| L2 | wtrmv, wtrsv | ✅ N/T/C | — | ✅ | ✅ | ❌ **(missed by 05-18 fix)** |
| L2 | wtbmv, wtbsv, wtpmv, wtpsv | ✅ N/T/C | — | ✅ | ✅ | ❌ |
| L3 | mgemm | ✅ ta,tb (via `rand_trans()` → N/T/C) | — | — | — | n/a |
| L3 | mgemmtr | ✅ ta,tb | — | ✅ | — | n/a |
| L3 | msymm | — | ✅ | ✅ | — | n/a |
| L3 | msyrk, msyr2k | ✅ N/T | — | ✅ | — | n/a |
| L3 | mtrmm, mtrsm | ✅ N/T/C (`rand_trans()`) | ✅ | ✅ | ✅ (30/70) | n/a |
| L3 | wgemm, wgemmtr, wsymm, whemm, wsyrk, wsyr2k, wherk, wher2k, wtrmm, wtrsm | ✅ (mirrors `m*`) | (per routine) | ✅ | (trmm/trsm only) | n/a |

Common helpers in `common/fuzz_util_body.fypp`: `rand_trans()` →
N/T/C ~1/3 each; `rand_trans_complex()` is just an alias; `rand_incx()`
→ {1,2,3,-1,-2} (added 05-18). No `rand_uplo`/`rand_side`/`rand_diag`
helpers exist — each fypp inlines its own 50/50 (or 30/70 for trmm/
trsm diag).

## Perf harness — categorical coverage

`scripts/gen_perf_harnesses.py` and generated harnesses under
`tests/blas_parallel/perf/target_*/`.

| Emitter | Routines | trans | side | uplo | diag | strides |
|---|---|---|---|---|---|---|
| emit_gemv | gemv | ✅ | — | — | — | ❌ |
| emit_ger | ger, gerc, geru | — | — | — | — | ❌ |
| emit_gbmv | gbmv | ✅ | — | — | — | ❌ |
| emit_symv_hemv | symv, hemv | — | — | ✅ | — | ❌ |
| emit_sbmv_hbmv | sbmv, hbmv | — | — | ✅ | — | ❌ |
| emit_spmv_hpmv | spmv, hpmv | — | — | ✅ | — | ❌ |
| emit_syr_her | syr, her | — | — | ✅ | — | ❌ |
| emit_spr_hpr | spr, hpr | — | — | ✅ | — | ❌ |
| emit_trmv_trsv | trmv, trsv | ✅ | — | ✅ | ❌ **N only** | ✅ `BLAS_PERF_INCX` |
| emit_tpmv_tpsv | tpmv, tpsv | ✅ | — | ✅ | ❌ **N only** | ❌ |
| emit_tbmv_tbsv | tbmv, tbsv | ✅ | — | ✅ | ❌ **N only** | ❌ |
| emit_gemm | gemm | ✅ R: {NN,TN,NT,TT}; C: 6 of 9 (skips CT/TC/CC) | — | — | — | n/a |
| emit_symm_hemm | symm, hemm | — | ✅ | ✅ | — | n/a |
| emit_syrk_herk | syrk, herk | ✅ | — | ✅ | — | n/a |
| emit_trmm_trsm | trmm, trsm | ✅ N/T/C | ✅ | ✅ | ❌ **N only** | n/a |
| emit_gemmtr | gemmtr | ⚠ 3 of 9 (ta,tb) pairs (NN, TN, NT) | — | ✅ | — | n/a |

## Gaps, ranked by risk

1. **L2 stride coverage (still open from 2026-05-18 audit).** 20 L2
   routine families × ~28 shared-body fypp variants still pass
   `incx=incy=1`. **Plus ~33 multifloats L2 fuzz files**: 16 in `m*`
   (mgemv, mgbmv, mger, msymv, msbmv, mspmv, msyr/2, mspr/2,
   mtrmv/mtrsv, mtbmv/mtbsv, mtpmv/mtpsv) and 17 in `w*` (wgemv,
   wgbmv, wgerc, wgeru, whemv, whbmv, whpmv, wher/her2, whpr/hpr2,
   wtrmv/wtrsv, wtbmv/wtbsv, wtpmv/wtpsv) — the 05-18 fix never
   touched these because they don't share the body include. Per the
   loop-direction survey, the cost of leaving this open is undetected
   correctness *or* perf cliffs in the strided fallback.
2. **DIAG='U' never benched.** trmv/trsv/tbmv/tbsv/tpmv/tpsv/trmm/trsm
   all hardcode `char diag = 'N'` in the emitted harness. Unit-diag
   takes a different inner-loop shape; a regression there would be
   invisible.
3. **GEMMTR perf samples 3 of 9 (ta,tb) pairs.** TT, CN, NC, CT, TC,
   CC unbenched. Per the loop-direction survey this matters because
   the (ta,tb) combination controls whether the inner walk is along
   rows or columns of A and B.
4. **GER family.** No categorical params (n/a) but also no stride
   randomization in fuzz or perf — gerc/geru complex-conjugate paths
   are exercised only at incx=incy=1.
5. **Real L2 triangular fuzz never exercises `trans='C'`.** trmv/trsv/
   tbmv/tbsv/tpmv/tpsv real bodies inline 50/50 N/T. For real types
   'C' is mathematically identical to 'T', but kernels often branch on
   the literal character; an `'C'` branch left to rot would surface
   only at link time against a real-typed caller. Cheap fix: replace
   the inline if/else with `rand_trans()`.
6. **Stride helper fragmentation.** `rand_incx` (added 05-18 in
   `common/fuzz_util_body.fypp`), per-fypp `rand_inc` (axpy R/C),
   `pick_inc` (~14 files), and `pick_inc_pos` (iamax/norm — for
   routines where negative stride is invalid per BLAS spec) coexist
   with subtly different distributions. Audit recommended unifying to
   `rand_incx`; not done. Also worth promoting `rand_uplo()`,
   `rand_side()`, `rand_diag()` to `common/` (all fypp files inline
   their own 50/50 — or 30/70 in the trmm/trsm diag case).

## Recommended sweep

Mirror the trsv/trmv playbook from 2026-05-18 in one bulk pass so
the audit/survey tables collapse to all-✅:

1. **Plumb `rand_incx()` / `rand_incy()` into all remaining L2 fuzz
   bodies AND the standalone multifloats L2 fuzz files (both `m*`
   real-DD and `w*` complex-DD).** Allocate buffers as
   `1 + (n-1)*abs(inc)`; thread strides through both kernel calls;
   add stride values to FAIL prints. Per the 05-18 workflow,
   sanity-test each by injecting a bug into one `incx ≠ 1` path and
   confirming non-zero failures, then revert.
2. **Extend the 10 affected `emit_*` functions** in
   `gen_perf_harnesses.py` to loop `BLAS_PERF_INCX` (and
   `BLAS_PERF_INCY` where applicable). Emit `/xK` key suffixes for
   non-unit strides. Use `perf_parse_int_list` (already added in
   `perf_common.h`).
3. **Add `BLAS_PERF_DIAG` loop** (default `{N, U}`) to
   `emit_trmv_trsv`, `emit_tpmv_tpsv`, `emit_tbmv_tbsv`,
   `emit_trmm_trsm`. Encode diag into the key (e.g. `ULU` vs `ULN`).
4. **Extend `emit_gemmtr` (ta,tb) sample set** to all 9 pairs. The
   loop-direction survey shows trans choice flips the inner walk
   direction; cheap to cover fully.
5. **Unify stride helper.** Keep `rand_incx()` in
   `common/fuzz_util_body.fypp`. Delete per-fypp `rand_inc()` /
   `pick_inc()` / `pick_inc_pos()`; redirect call sites. Add
   `rand_uplo`/`rand_side`/`rand_diag` while at it.
6. **Switch L2 real triangular bodies to `rand_trans()`** so 'C' is
   covered (one-line per file).
7. **Regenerate** all perf C files via
   `uv run scripts/gen_perf_harnesses.py`.
8. **Rebuild + re-run** fuzz/perf across kind10/kind16/multifloats.
   Catalog any sub-parities as new Addenda in
   `doc/parallel-blas-optimization-findings-20260515.md`, mirroring
   the loop-direction survey's Addendum + Rule format.

Estimated scope: ~28 shared-body fypps (done 2026-05-18 follow-up
session) + ~33 multifloats fypps (16 m* + 17 w*) + ~10 perf emitters
+ ~150 perf C regen. Single sweep; the survey tables collapse on
completion.

## What this survey does *not* cover

- **`alpha` / `beta` scalar values.** Most fuzz bodies use
  `rand_alpha_beta()` (samples ±1, 0, random). Not separately
  audited.
- **Edge sizes.** N=0, N=1, M=0, n=k for banded. `rand_int_log(1,...)`
  floors at 1, so N=0 is not fuzzed. Separate concern (same as 05-18).
- **lda alignment.** `rand_pad()` covers 0-128 padding; whether it
  stresses cache-line/SIMD boundaries specifically is out of scope.
- **Buffer aliasing.** Mixed-stride overlap (e.g. gemv with x and y
  pointing into the same allocation) — BLAS spec says undefined, but
  worth confirming no silent corruption. Out of scope.
