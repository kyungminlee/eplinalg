# Parallel BLAS overlay — optimization findings

Date: 2026-05-15
Branch: `parallel-blas`
Hardware: Intel i3-1315U (Raptor Lake-U, Alder Lake-derived; 2P+4E cores, P-core L1d=48 KB, L2=1.25 MB, L3=10 MB shared)

This document records what was tried, what worked, and what didn't, while tuning the GEMM overlays. It complements `parallel-blas-20260513.md` (the design doc). Read this when you're about to "obviously" optimize a kernel — many of the obvious moves have already been tried and rejected on bench.

The unifying theme: **structural choices that paid off for one (precision, complex-ness) pair often regressed another**. Don't assume.

---

## TL;DR results table

GFLOPs at OMP=1, median speedup over migrated, range across 9 transpose combos × sizes 64/128/256(/512):

| routine | precision  | shape                  | median over migrated | shape of curve   |
|---------|------------|------------------------|----------------------|------------------|
| egemm   | kind10 real    | OpenBLAS-style packed + MR=2/NR=2 tile + TN fast-path | **1.70×** | wins everywhere ≥ s=128 |
| ygemm   | kind10 complex | reference + OMP-over-j, rank-1 paths K-unrolled by 2 | **1.00×** | parity, zero regressions |
| qgemm   | kind16 real    | reference + OMP-over-j (libquadmath-bound)            | ~1.0× | parity, libquadmath dominates |
| xgemm   | kind16 complex | reference + OMP-over-j                                 | ~1.0× | parity |
| mgemm   | mf real    | full GotoBLAS structure, MR=4 register tile via AVX2 FMA | **~28×** | hits 3.86 GF/s asymptote at s≥128 |
| wgemm   | mf complex | same as mgemm with 4-buf SoA pack (rh/rl/ih/il)        | **~19×** | hits 3.9 GF/s asymptote at s≥128 |

---

## Findings by precision

### kind10 real (egemm) — full OpenBLAS structure pays off

The original packed kernel sat at ~1.6 GFLOPs across all transposes, with migrated zgemm reference hitting 2.4 GFLOPs on `TN`/`CN` due to stride-1 access.

What worked:

1. **GotoBLAS three-level cache blocking + copy-and-conquer packing** (MC=64, KC=256, NC=512). Same shape as OpenBLAS — `js → ls → is` loop order, pack B once per (jc, pc), pack A per (ic, pc).
2. **MR=2, NR=2 register-tile outer-product micro-kernel** with 4 scalar `long double` accumulators kept on the x87 register stack across K. The compiler emits 4 independent FMA chains per K iter, masking x87 fmul latency. **This single change lifted asymptote from 1.6 → 2.2 GF/s (+40%)**.
3. **Adaptive MC when K ≤ KC**: grow `MC` so `MC*KC ≈ L2_target` (768 KB target, 4× cap). Helps small-K shapes.
4. **Explicit `TA='T', TB='N'` reference fast-path**: when A is transposed and B is not, the inner k-loop already accesses both A and B stride-1 — packing adds cost without recovering anything. The reference DDOT loop wins at TN/CN by ~5–7% over packed. Kept this fast-path despite "OpenBLAS doesn't have it" because OpenBLAS has an asm kernel that we can't match.

Result: median 1.70× over migrated across 36 (transpose, size) combos; only outlier is CN s=64 at 0.88× (single-run noise).

### kind10 complex (ygemm) — packed structure *regressed*; reverted to reference+OMP

The same MR=2 NR=2 register-tile that worked for egemm **regressed** ygemm from 1.6 → 1.5 GF/s.

**Why**: each `_Complex long double` multiply needs ~6 fp80 temp slots on the x87 stack. Four complex accumulators = 8 fp80 slots = the entire 8-deep x87 register stack, leaving no room for inputs → every input load spills to L1. The win that comes from "scalar accumulators kept in registers across K" simply doesn't exist for complex long double.

Replaced with **reference algorithm + `#pragma omp parallel for` over j**, four explicit orientation branches:

- `NN` rank-1 update (`temp = alpha*B[l,j]; c[i,j] += temp * A[i,l]`)
- `(T|C)N` dot product
- `N(T|C)` rank-1 update with B-row broadcast
- dual-transpose dot product

What further worked: **K-axis manual unroll by 2 on the rank-1 paths** (NN, NT, NC):

```c
for (l = 0; l + 1 < K; l += 2) {
    t0 = alpha * b[j2*ldb + l];
    t1 = alpha * b[j2*ldb + l + 1];
    al0 = &a[l * lda];
    al1 = &a[(l + 1) * lda];
    for (i = 0; i < M; ++i)
        cj[i] += t0 * al0[i] + t1 * al1[i];
}
```

This matches what gfortran emits for Netlib's reference zgemm rank-1 loop — two independent FMA chains per `i` mask the ~10-cycle x87 fmul latency. Recovered NN/NT/NC from **0.78× → 1.08×**.

What did *not* work and was reverted:
- **K-axis unroll on the dot-product paths** (T*N, dual-transpose): single-accumulator form was already well-scheduled by gcc; splitting to `acc0 + acc1` regressed CT/CC to 0.35×.
- **I-axis unroll on rank-1**: gcc already pipelines independent C cells across rows; manual i-unroll buys nothing.

Result: median 1.00×, range 0.97×–1.09×, **zero regressions** across all 36 combos.

### kind16 (qgemm, xgemm) — leave alone

Both already use reference + OMP-over-j. libquadmath function-call overhead per `__multf3` dominates per-op time (~hundreds of cycles), so structural optimizations (blocking, register tiling, K-unroll) **cannot help** — there's no latency to hide and no cache pressure to ease.

The K-axis unroll trick that helped ygemm rank-1 does not transfer: libquadmath calls don't pipeline within the CPU reorder window the way x87 fmuls do. Two `__multf3` calls "in flight" run sequentially.

Bench confirms: q/xgemm sit at 0.95–1.04× of migrated. No change tried; none expected to help.

### multifloats (mgemm, wgemm) — already OpenBLAS-style, one good knob

Both implement the full GotoBLAS structure from day one:
- SoA-packed Ap (column-major MR-row tile) and Bp (4 separate arrays for wgemm: rh/rl/ih/il)
- Templated MR×NR_PAN register-tile SIMD micro-kernel using `simd_dd::dd_mul` + `dd_add` (AVX2+FMA intrinsics from `mgemm_simd_kernel.h`)
- Sub-NC chunking by NR×NR_PAN
- Loop order js → ls → is
- BETA pre-pass

**Only one change shipped**: `MBLAS_SIMD_MR` default 3 → 4. Bench at 80 iters pinned to a P-core:
- s=128: +5.07% (stdev 0.38%)
- s=256: +4.65% (stdev 0.20%)
- s=64: −1.22% (stdev 3.09% — statistically indistinguishable)

MR=4 uses 8 acc + 4 broadcasts + 2 B-loads = 14 ymm of 16 available. The bigger tile amortizes per-tile setup (writeback, register fills) better at s≥128.

---

## Optimizations that did NOT work

These were tested empirically and rejected. **Don't re-try without a different motivating observation.**

### 1. Persistent thread-local scratch buffers (multifloats)
**Hypothesis**: at s=64, three `aligned_alloc` + `free` pairs per call (Ap=128 KB, Bp_hi/lo=256 KB each) dominate per-call overhead.

**Test**: replaced with `static thread_local Scratch { T*Ap; double*Bp_hi; double*Bp_lo; ... }`, sized on first call and reused.

**Result**: s=64 NN 2.56 → 2.55 GF/s. **Zero change**.

**Why**: glibc's tcache + large-bin recycling makes the 128–256 KB alloc/free pairs ~1 µs after warmup. Not the bottleneck.

### 2. OMP region disabled (multifloats)
**Hypothesis**: `GOMP_parallel` setup/teardown adds 5–20 µs/call even at `OMP_NUM_THREADS=1`.

**Test**: defined `MGEMM_DISABLE_OMP` to `#ifdef`-out all four `_OPENMP` guards in `mgemm.cpp`. Compiled, rebuilt, bench.

**Result**: NN s=64 +11% (but NT, NC, CN: 0% to −3%). Mean of 9 combos: **statistical noise**.

**Why**: libgomp's single-thread path is near-trivial. Not a real source of per-call overhead.

### 3. MR=3 tail-row hypothesis (multifloats)
**Hypothesis**: at M=64 with MR=3, one tail row of M%MR=1 goes through the scalar fallback path inside `inner_kernel_simd_t`. Eliminating it (e.g. MR=4 divides 64 evenly) should help s=64.

**Test**: rebuilt with `MBLAS_SIMD_MR=4`. Bench at 80 iters.

**Result**: s=64 actually **regressed** by 1.2% (within noise). s≥128 improved by 5%.

**Why**: bigger tile has more register pressure (8 acc + 4 broadcasts = 12 ymm before B-loads). The win from no-tail-row at s=64 is offset by reduced register headroom and scheduling slack. The tail row was never the bottleneck.

### 4. simd_writeback scratch round-trip avoidance (multifloats)
**Hypothesis**: `_mm256_store_pd(ph_a, ph); ... r.limbs[0] = ph_a[j];` causes a SIMD-store→scalar-load round-trip with store-forwarding stall.

**Test**: replaced with gcc vector-element indexing `((v4df)ph)[j]` to extract lanes directly from the ymm register, no memory round-trip.

**Result**: **0% change** across all sizes (±0.3% noise).

**Why**: Alder Lake's store-forward path handles this case well. The round-trip wasn't the bottleneck — the **scalar `dd_add` on column-strided C** was. Confirmed by stubbing `simd_writeback` entirely to a no-op: that gave +6.7% at s=256, all of which is in the C read-modify-write itself, not the alpha-scale or lane extract.

### 5. Adaptive MC for multifloats
**Hypothesis**: same heuristic that helped egemm — grow MC when K is small so MC*KC ≈ L2.

**Test**: applied to mgemm with `if (K <= KC) target_mc = L2_TARGET / (K * sizeof(T))`. Bench at s=64, 128.

**Result**: s=64: noise (M=64 is already ≤ MC, so the grown MC is capped to M in the ic loop — no effect). s=128: −1 to +4%, noise.

**Why**: mgemm's MR=3 (4) micro-kernel hits asymptote at s=128 already; adaptive MC doesn't move the L2 working set into a faster regime at these sizes.

### 6. MR×NR register tile for kind10 complex (ygemm)
**Hypothesis**: same MR=2 NR=2 tile that worked for egemm should work for ygemm.

**Test**: implemented and benched.

**Result**: 1.6 → 1.5 GF/s, **slight regression**.

**Why**: see kind10 complex section above — x87 stack capacity exhausted by complex multiply temps.

### 7. K-axis unroll on ygemm dot-product paths
**Hypothesis**: K-unroll helped rank-1 (NN/NT/NC); should also help dot-product paths.

**Test**: split `acc` into `acc0, acc1`, accumulate each over l, l+1 strides.

**Result**: TN/CN regressed 1.00× → 0.80×; **CT/CC regressed 0.99× → 0.35–0.56×**.

**Why**: gcc's scalar single-acc loop was already well-scheduled; the split disrupted its register allocation. The conditional `conj_b ? ~bv0 : bv0` inside an unrolled hot loop wasn't constant-folded, putting a runtime branch in the critical path.

---

## Diagnostic methodology that worked

When chasing a perf gap, follow this order:

1. **Measure the actual time/cycle cost first.** Compute `cycles/flop = freq / GFLOPs`. The difference between observed and asymptote gives the absolute per-call overhead in cycles or µs. For mgemm at s=64: 2.5 vs 3.7 GF/s = 0.58 extra cycles/flop × 524 K flops = ~67 µs/call of fixed overhead.

2. **Get IPC from `perf stat`.** Pin to a P-core with `taskset -c 0`. If IPC is identical at small and large sizes, the kernel itself is steady-state at both — the gap is in non-kernel paths or per-call setup. (mgemm: IPC = 1.97 at s=64, 1.95 at s=256 — same.)

3. **Test one hypothesis at a time with a single-variable change**, ideally one that's structurally reversible. Don't bundle.

4. **Rule things out with ablation** (stub-to-no-op), not just optimization. Stubbing `simd_writeback` told us writeback is ~7% of total — useful even though we can't fully eliminate it.

5. **Bench at high iter counts (80+) pinned to a P-core.** Below that, run-to-run variance is ±3–5% and small effects vanish in noise. Without pinning, the kernel may migrate between P-cores and E-cores, contaminating cycle counts.

6. **Suspect noise on `(s=N, single combo)` jumps; trust `(mean over 9 combos at s=N, stdev < 1%)` for structural effects.**

---

## Hardware-specific gotchas on this machine

- **Hybrid CPU**: `perf stat` reports `cpu_atom` (E-core) and `cpu_core` (P-core) separately. Threads can migrate; pin with `taskset -c 0` (or 1, 2, 3 — P-cores 0–3, E-cores 4–7).
- **P-core L2 = 1.25 MB**. Adaptive-MC targets should aim for ~768 KB working set to leave headroom for Bp + code/stack.
- **No AVX-512** on the U-skus — AVX2 + FMA is the ceiling. SIMD-DD kernels use ymm registers (16 available).
- **x87 long-double FPU** is the only path for `long double` and `_Complex long double`. 8-deep register stack. No SIMD. ~10-cycle fmul latency → 2 independent FMA chains saturate it.
- **System `perf` mismatched** with running kernel — use `/usr/lib/linux-tools/6.8.0-111-generic/perf` (any close-enough version works; some events show `<not supported>`).

---

## Open items / things worth trying later

These weren't tested but might be worth a controlled bench:

1. **Specialized small-size kernel for mgemm s=64**: fully unrolled MNK=64 with no blocking, possibly compile-time generated. Target: close the s=64 vs s=128 gap (currently 2.5 vs 3.6 GF/s).
2. **MR=2, NR_PAN=2 shape** for mgemm: same tile area as MR=4, NR_PAN=1 but different shape. Different register pressure profile (8 acc + 2 broadcasts + 4 B-loads = 14 ymm), might suit some access patterns.
3. **Software prefetching** for Ap/Bp in inner kernel: OpenBLAS does this; we don't. May help at s≥256 where data starts to exceed L2.
4. **Bench mgemm/wgemm at OMP=4 or 8** — adaptive MC and persistent buffers may behave differently with cross-thread cache contention.
5. **Egemm + ygemm at OMP>1** — current numbers are OMP=1. Threading is the actual point of the overlay.

---

# Addendum: kind10 L3 (symm/hemm/syrk/herk) round — 2026-05-15

A second optimization pass on the symmetric and Hermitian L3 routines, after the GEMM/TRSM work above. Targets: esymm, ysymm, yhemm, esyrk, ysyrk, yherk. Pattern: diagonal-block scalar kernel plus an off-diagonal trailing `[ey]gemm` call per block-pair, blocked over column panels.

The original baseline had several catastrophic outliers — esyrk U/N s=64 at 0.24×, yhemm L/L s=64 at 0.43×, ysymm L/U at 0.75–0.95× — buried inside otherwise-passable mean speedups (1.0–1.3×). Three structural fixes, applied as a recipe to each file, cleared all of them. Plus one negative result and one "no further closure possible" finding worth recording.

## TL;DR results (OMP=1, 40 iters, P-core)

| routine | mean before | mean after | min before | min after |
|---------|-------------|------------|------------|-----------|
| esymm   | 1.16× | **1.41×** | 0.85 | 0.89 |
| ysymm   | 1.00× | **1.20×** | 0.75 | 1.09 |
| yhemm   | 0.99× | **1.18×** | **0.43** | 0.93 |
| esyrk   | 1.22× | **1.31×** | **0.24** | 0.85 |
| ysyrk   | 1.06× | **1.12×** | 0.82 | 0.99 |
| yherk   | 1.05× | **1.13×** | 0.80 | 0.97 |

## What worked

### A. Stride-1 A access in the diagonal kernel (symm/hemm)
The Netlib reference DSYMM/ZSYMM/ZHEMM walks A column-major — for SIDE='L' UPLO='L', `A(K,I)` with `K>I` reads the stored lower triangle stride-1 in K. Our overlay was reading `A_(i, k)` (row i, column k) — the *symmetric counterpart*, also stored, but **stride-lda in k** at fixed i. Every inner-iteration A load crossed a cache line.

Fix: reverse the i loop direction (backward for UPLO='L', forward for UPLO='U') so the inner k loop reads `A_(k, i)` — same column, varying row, stride-1.

For Hermitian: the conjugate moves from the cj[k] update onto the temp2 accumulation when iteration direction flips (since `conj(A(i,k)) = A(k,i)` for the stored half). Same total ops, same result, stride-1.

Impact on its own: ysymm L/U 0.75–0.95 → 1.12–1.18; yhemm L/L 0.43 → 0.93 at s=64.

### B. Single-block fast path (symm/hemm)
At small M (≤ nb), the framed `for(jc) → for(ic) → diag_kernel + gemm` path runs exactly one ic iteration with no gemm call. The framing — `omp_get_max_threads`, OMP pragma runtime check, function call into the static diag kernel, jc/ic outer loops — adds ~5–10% overhead on top of an irreducible scalar inner loop.

Fix: early fast path inside the entry function that inlines the Netlib loop directly, folds BETA into the diagonal write, and skips all framing.

Modest additional gain (~3–5%) and removes the BETA pre-scale separate pass.

### C. Smaller `nb` cap (all six)
Default `nb=64` (esymm/esyrk) or `64–256` (ysymm/yhemm/ysyrk/yherk). At s=128 with nb=64, the entire problem is a single diagonal block — zero gemm trailing-update work, 100% of the time in the scalar diag kernel. The gemm-driven speedup that justifies this overlay never fires.

Dropped to:
- **esymm: cap 128 → 64**. egemm at K=32 doesn't amortize for real long-double (packing overhead too high), so the floor stays at 64.
- **ysymm/yhemm/ysyrk/yherk: cap 256 → 32**. Complex long-double per-element work is heavy enough that ygemm tolerates K=32; packing fraction is small.
- **esyrk: cap stays 64** (same real long-double reasoning).

Impact: s=128 and s=256 fragment into 2×2 and 4×4 block grids with the trailing `[ey]gemm` calls picking up most of the off-diagonal work.

### D. `restrict` pointers on all A/B/C parameters (largest single per-file win)
C `const T*` doesn't tell gcc the buffers don't alias. The diag inner loop is `cj[k] += temp1 * ai[k]; temp2 += bj[k] * cconj(ai[k]);` — gcc has to assume the `cj[k]` store might have changed `ai[k]` or `bj[k]`, so it reloads them from memory on every iteration. Verified by objdump: `fldt (a.real); fldt (a.imag); ... fstpt (cj.real); fstpt (cj.imag); fldt (a.real); fldt (a.imag);` — three redundant loads per iteration.

Adding `restrict` to `a`, `b`, `c` on both the entry-point signatures (yhemm_, ysymm_, esymm_, esyrk_, ysyrk_, yherk_) and the static helper signatures eliminates the reloads. The BLAS spec already forbids aliasing, so this is a free correctness assertion.

Impact: ysymm L/L 1.04 → 1.09, yhemm L/L 0.93 → 0.96, ysymm L/U all sizes +1–2%. Biggest on the complex routines where the inner loop has to keep both real/imag halves of A live across the C store. esymm (real long-double) was flat — a single scalar write doesn't trigger the same reload chain.

## What didn't work

### Manual `__real__`/`__imag__` real-arithmetic expansion (yhemm)
After noticing that `bj[k] * cconj(ai[k])` compiles to ~25 x87 ops (verified by objdump), tried rewriting the inner loop as explicit scalar `long double` arithmetic:

```c
__real__ cj[k] += t1r * ar - t1i * ai;
__imag__ cj[k] += t1r * ai + t1i * ar;
t2r += br * ar + bi * ai;
t2i += bi * ar - br * ai;
```

**Result**: regressed L/L 64 from 0.93× → 0.86×. Disassembly showed gcc had successfully eliminated the redundant A/B loads, but spilled `t1.i` to the stack (`fldt 0x10(%rsp)`) instead and reloaded it twice per iteration. The x87 stack is 8 slots and the expanded form needs to keep more values live than the implicit form. Net: more loads, slower.

`-fcx-fortran-rules` already gets gcc to inline complex multiplication (no `__mulxc3` libcall), so the manual expansion buys nothing the compiler isn't already doing.

### gcc scheduler flags `-fschedule-insns -fsched-pressure` (yhemm)
Both default-off on x86 due to historic x87 quirks. Tried enabling them per-file via `set_source_files_properties`. Result: small regression (L/L 64: 0.93 → 0.89). Tried `-fmodulo-sched` separately: net negative too. The default scheduler is already doing the best it can on this inner loop.

## Findings worth recording

### x87 codegen ceiling on `_Complex long double` (yhemm L/L)
Even with all fixes applied, yhemm L/L sits at 0.92–0.96× across all sizes — overlay ~2.15 GF/s, migrated ~2.30 GF/s. The 7% gap is a hard floor: gcc and gfortran are producing **the same inlined-cmul + fchs-conjugate inner-loop structure** (verified via objdump on both `.o` files); the difference is x87 register scheduling, which gfortran does slightly better on this specific pattern.

Closing the gap would require:
- **Type-level rewrite**: SoA (`long double` re/im arrays) instead of `_Complex long double`. The Fortran ABI forbids changing the public signature, so this only works as a local repack into scratch buffers inside the kernel — and the ~2 KB unpack/repack per 32×32 diag block likely eats the gain.
- **Inline x87 assembly**: technically feasible but gcc's extended asm with x87 constraints is famously buggy; you'd be chasing gfortran's already-optimal schedule with a more fragile encoding.

Not pursued. The remaining 7% is below the noise floor of any realistic OMP>1 win.

### The `nb` choice splits real vs complex
Real long-double (`egemm`, `esyrk`): keep nb ≥ 64. egemm at K=32 doesn't amortize its packing overhead — the scalar inner kernel has enough ILP on x87 that the gemm-via-pack isn't a win until K is large.

Complex long-double (`ygemm`, `ysyrk`, etc.): nb=32 is fine. Each complex op is ~8× heavier than scalar long-double, so the packing fraction of total work shrinks and gemm pays off at K=32.

### Aliasing pessimism is real and `restrict` is free
This is the most important takeaway from the L3 round. Three lines of code (`restrict` on a, b, c at the entry-point + helper) closed several percentage points of margin across complex routines. The disassembly proof made this obvious — but it would have been easy to miss had I not been looking at the actual generated code.

Whenever a BLAS-shape kernel has the inner-loop pattern `c[i] += stuff(a[i], b[i])`, `restrict` on those three pointers is almost always a free win on x87 long-double — because the compiler can otherwise assume the C write may have changed A or B.

## Additions to "diagnostic methodology"

7. **When the overlay loses on a single (uplo, trans, size) cell, check whether the framed path is hitting a *single-block regime* — i.e. nb is large enough that the gemm trailing-update never fires.** That's almost always where the catastrophic outliers live (esyrk U/N 0.24×, yhemm L/L 0.43×, ytrmm L/U/N/N 0.45×). The fix is rarely "improve the diag kernel" — it's "shrink nb so the trailing gemm picks up the slack."

8. **For x87 long-double kernels, always check the objdump for redundant `fldt` after a `fstpt` to A/B addresses.** Two reloads of the same memory location per inner iteration is a tell that the compiler is being conservative about aliasing. `restrict` fixes it.

9. **Negative results (this round): two compiler flag experiments and one source-level manual expansion** all regressed. Save effort: prove the disassembly is suboptimal *before* trying to fix it, and verify the fix changes the disassembly in the expected direction *before* benching.

10. **Compare migrated to the theoretical x87 ceiling before "fixing" the overlay.** Real long-double dot product has three regimes:
    - Single-accumulator: ~1.3 GF/s (7-cycle fadd dep chain)
    - 2 independent accumulators: ~2.7 GF/s
    - Fully pipelined fmul+fadd, no chain: ~9.4 GF/s (rarely achievable)

    Migrated DSYRK at TR='T' s=64 hits 2.36 GF/s — that's right at the 2-acc ceiling. **If migrated is at a known x87 ceiling, you can't beat it; you can only match it by replicating the same dep-chain count.** Compute migrated's % of peak first; if it's ≥85% of some ceiling, identify which ceiling and replicate.

## K-axis 2-unroll: works on real dot products, breaks on complex ones

Counterpart to the documented ygemm dot-product regression (§7 above). The pattern:

```c
T s0 = 0.0L, s1 = 0.0L;
int l = 0;
for (; l + 1 < K; l += 2) {
    s0 += Ai[l]     * Aj[l];
    s1 += Ai[l + 1] * Aj[l + 1];
}
T s = s0 + s1;
for (; l < K; ++l) s += Ai[l] * Aj[l];
```

**Real long-double (esyrk TR='T')**: cleanly doubles throughput on a serial fadd-chain kernel. esyrk L/T 64: 0.85× → 0.98×; every other L/T,U/T cell pinned to ≥0.98×. Matches gfortran's reference DSYRK rate.

**Complex long-double dot product paths (ygemm TN/CN, CT/CC; ysyrk/yherk T/C)**: regressed under the same pattern, or already at parity so the unroll buys nothing. Two reasons:
1. Complex multiplication is itself ~4× heavier than real, so the single-chain bottleneck is less of the runtime (chain latency hidden by cmul work).
2. ygemm had a `conj_b ? ~bv0 : bv0` conditional inside the loop body that wasn't constant-folded; the K-unroll put a runtime branch in the critical path. Even without that, naïve K-unroll on complex would help less.

**Rule of thumb:** K-axis 2-unroll is the right structural move on **real scalar long-double dot-product reductions**. For complex types, leave the inner loop single-accumulator — gcc's scheduler does better when there's only one cmul chain to lay out on the x87 stack.
