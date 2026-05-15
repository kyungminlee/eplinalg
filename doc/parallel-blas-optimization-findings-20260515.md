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
