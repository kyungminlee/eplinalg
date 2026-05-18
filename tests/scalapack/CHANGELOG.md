# tests/scalapack — CHANGELOG

Resolved items, reverse-chronological. Open work lives in `TODO.md`.

## 2026-05-05 — p?larz / p?larzc PBxTRNV M-vs-L fix

Lifted upstream commit `c408393` ("fix PBxTRNV N-arg in p?larz /
p?larzc") into `pdlarz.f`, `pzlarz.f`, and `pzlarzc.f` overrides:
the row-V SIDE='L' and column-V SIDE='R' branches passed `M` /
`N` as the global vector length to PBxTRNV, but only the trailing
`L` entries of V are physically stored — same pattern as the
P?LARZB fix landed earlier. Changed to `L` in 4 sites per file.

Not observed in our test matrix because our companion MPV/NQV
redefinition constrains the consuming GEMV to the meaningful L
local rows, but the OOB read in PBxTRNV is unsound on its own
and aborts on NPROW=3 upstream. Documented in
`doc/upstream-bugs/README.md`.

## 2026-05-05 — pzunmrz SIDE='L' resolved (complex path)

Found and fixed a complex-specific algorithmic bug in `PZLARZ` and
`PZLARZC` SIDE='L' branches that caused `pzunmrz` SIDE='L' to fail
with O(1) residuals (~1.3-1.8x) on every grid configuration.

`PZLARZ` accumulates `w = v^H * sub(C)` via
`ZGEMV('Conjugate transpose', ...)` followed by `ZAXPY` of the
'1'-row of C, and applies the rank-1 update via `ZGERC`. But
`ZGEMV('C', C2, V)` produces `C2^H * V = conj(V^H * C2)`, not
`V^H * C2`, and `ZGERC` then conjugates the second operand again —
so the imaginary part of `v^H * C` is sign-flipped twice
inconsistently. LAPACK's serial `ZLARZ` handles this by calling
`ZLACGV` twice (before adding C1, then again after) and using
`ZGERU` for the rank-1 update; `PZLARZ`/`PZLARZC` translated the
ZGEMV trick but skipped the `ZLACGV` cleanup AND used `ZGERC`.
Real path `PDLARZ` works because `DGEMV`/`DGER` have no
conj/no-conj distinction.

Override applied to all 5 SIDE='L' MPV×NQC2 sites in each of
`pzlarz.f` and `pzlarzc.f`: insert `CALL ZLACGV(NQC2, accum, 1)`
after the `ZGEMV` and replace `ZGERC(MPV, NQC2, ...)` with
`ZGERU(MPV, NQC2, ...)`. SIDE='R' branches are unchanged
(they correctly use `ZGERC` because `H` from the right is
`C := C * (I - tau v v^H) = C - tau (Cv) v^H` where the second
operand really is conjugated).

**Test impact.** `test_pzunmrz` now exercises SIDE in {L, R} ×
TRANS in {N, C} (4/4 PASS on kind16, residuals ~1e-33 on 1, 2,
and 4 ranks). `test_pztzrzf` continues to pass (uses `PZLATRZ` →
`PZLARZ` SIDE='R', branch unchanged). Documented in
`doc/upstream-bugs/README.md`. Not covered by any existing upstream
`fix-*` branch in `../scalapack-bugfix/scalapack`.

## 2026-05-05 — pdormrz SIDE='L' resolved (real path); pzunmrz SIDE='L' narrowed

Carried in two upstream bug-fix branches from
`scalapack-bugfix/scalapack` and added a companion fix for the
conjugate variant:

- `fix-larz-buffer-sizing` (commits `0a58017` MPV / `c429b7b` NQV):
  `P?LARZ` derived the receiver-buffer size for the SIDE='L' / SIDE='R'
  PBxTRNV transpose from V's distribution; the actual write goes to
  sub(C2)'s distribution. Under PDORMR3's IV iteration the alignment
  restriction is violated, the buffer is undersized by 1-3 elements,
  and the unpacked write into Y overruns into adjacent aux WORK —
  observed as glibc "corrupted size vs. prev_size" aborts on
  NPROW > 1 grids with non-uniform row distribution.
- `fix-larz-daxpy-stride` (commit `36abce7`): the six (real) /
  five (complex) `?AXPY` calls in `P?LARZ` that target `WORK` /
  `WORK(IPW)` used `INCY = MAX(1, NQC2)` (a leading dimension) where
  `WORK` is a contiguous vector — correct stride is `1`. Looks like
  copy-paste of the legitimate LD on the surrounding `?LASET` /
  `?GSUM2D` matrix calls.
- Override `pzlarzc.f` (the `Q**H` variant called from `PZUNMR3` when
  `TRANS='C'`) carries the same MPV/NQV + AXPY-stride fixes — the
  byte-identical formulas exist there too, but the upstream bug-fix
  branches did not patch the conjugate variant.

Wired via `recipes/scalapack/source_overrides/p[dz]larz.f`,
`pzlarzc.f`, and `prefer_source: PDLARZ, PZLARZ, PZLARZC` pins.
Documented in `doc/upstream-bugs/README.md`.

**Test impact.** `test_pdormrz` now exercises SIDE in {L, R} ×
TRANS in {N, T} (4/4 PASS on kind16 / 2×2 grid; previously SIDE='R'
only). `test_pzunmrz` SIDE='R' continues to pass; SIDE='L' was
parked as a complex-only residual (~1.3-1.8x), then resolved
later the same day — see the entry above.

## 2026-05-03 — PDDBTRS / PZDBTRS workspace under-allocation (`*trsv` oracle)

Extended the LWMIN override (originally for `*dbtrs`) to every
`pddbtrs` / `pzdbtrs` call site, including the oracle call in
`test_p[dz]dbtrsv` — trusting the LQUERY result there leaks the same
upstream under-allocation and corrupts the heap (surfaces at
`MPI_Finalize` → libfabric `dlclose` → `free`).

## 2026-05-02 — `*trsv` (single-half banded/tridiag triangular solves)

The 8 smoke tests at
`tests/scalapack/linear_solve/test_p[dz]{pt,dt,pb,db}trsv.f90` were
upgraded from INFO=0 checks to validation tests that compose the two
halves and compare against the full `*trs` oracle on the same
distributed factor. Per-family composition:

* **PT** (real LDL^T): `pdpttrsv('L')` → mid-step → `pdpttrsv('U')`
  ≡ `pdpttrs`. The mid-step mirrors `pdpttrs.f:707-716` exactly: ranks
  0..NPCOL-2 divide entries 1..nb-1 by `d_loc` (the modified LDL^T
  diagonal D'), then divide entry nb by the reduced-system pivot
  stored at `af(nb+1)`; the last rank divides all nb entries by
  `d_loc` (no boundary). The original "fails by 60-200x" smoke test
  was missing both the D'^-1 step and the AF reduced-system pivot.
* **PT** (complex Hermitian LDL^H): `pzpttrsv('L','N')` → mid-step
  (real D' divide + complex AF pivot divide, mirroring
  `pzpttrs.f:744-758`) → `pzpttrsv('L','C')` ≡ `pzpttrs(UPLO='L')`.
  Note backward half uses `('L','C')` (apply L^H of stored L), not
  `('U','N')` — pzpttrf stores the L view, and `pzpttrs(UPLO='L')`
  internally uses the same pattern.
* **DT/DB** (LU): `trsv('L','N')` ∘ `trsv('U','N')` ≡ `*trs` (no
  scaling). For DB, the `*trsv` and `*trs` paths must use distinct AF
  copies — even though AF is documented read-only, in practice the
  `*trs` divide-and-conquer path perturbs internal scratch in AF on
  certain ranks, so the test snapshots AF after the factor and runs
  the two paths sequentially with restore in between. Workspace must
  be the bumped formula `nrhs*max(bwl,bwu) + (nrhs-1)*(bwl+bwu)`.
* **PB** (Cholesky, UPLO='U'): `trsv('U','T'/'C')` then `trsv('U','N')`
  ≡ `pbtrs` — apply U^-T (or U^-H) before U^-1.

All 8 PASS to ~target precision on kind16 / 2×2 grid.

## 2026-05-02 — pdtrsen LQUERY heap corruption

Root cause was an upstream LQUERY-contract bug: `pdtrsen.f:499-538`
writes `IWORK(1:N)` (the SELECT→integer copy plus the IGAMX2D broadcast)
before the LQUERY return at line 619-622, so a caller passing
`IWORK(1)` for the workspace query writes past the end of its buffer
and corrupts the heap. The follow-up "real" call ran cleanly
(LIWMIN ≥ N), which is why the eigenvalue spectrum was correct even as
the heap was already damaged.

Fix: `target_pdtrsen` in
`tests/scalapack/common/target_scalapack_body.fypp` allocates a local
`iwork_t(max(1,n))` in the LQUERY branch and forwards it to upstream in
place of the caller's IWORK; `iwork_t(1)` (= LIWMIN) is copied back to
`iwork(1)` afterwards. Mirrors the `WORK_T(3)` pattern for
`target_pdsyevx`'s upstream `WORK(1:3)` early-write.

New driver `tests/scalapack/eigenvalue/test_pdtrsen.f90` exercises the
fix over sizes [32, 64, 96]; all PASS without the previous `free()`
abort. Documented in `doc/upstream-bugs/README.md`.

## 2026-05-01 — pdposvx / pzposvx LWMIN bugs

Two upstream LWMIN bugs in `pdposvx.f` / `pzposvx.f`. Both fixed via
`recipes/scalapack/source_overrides/p[dz]posvx.f` plus
`prefer_source: PDPOSVX, PZPOSVX` in `recipes/scalapack.yaml`.

1. `LWMIN = 3*DESCA(LLD_)` is too small for the internal PDPOCON /
   PZPOCON contract (param-10 abort). Documented contract is
   `MAX(PDPOCON_LWMIN, PDPORFS_LWMIN)`.
2. Companion bug in `pzposvx.f`: `LRWMIN = MAX(2*NQ, NP)` doesn't cover
   PZLANHE('1', ...) which needs `2*NQMOD + NPMOD + LDW`. PZLANHE then
   writes past RWORK and corrupts the heap, surfacing as a `free()`
   SIGSEGV during cleanup.

Tests: `tests/scalapack/linear_solve/test_pdposvx.f90` /
`test_pzposvx.f90` both 3/3 PASS at ~33-digit accuracy on kind16 /
2×2 grid. Documented in `doc/upstream-bugs/README.md`.

## 2026-05-01 — pdgecon / pdpocon condition-number estimators

Implemented the κ-true reference helper and Hager-Higham bound check.
`tests/scalapack/common/cond_helpers.f90` exposes
`true_kappa1_general` / `true_kappa1_general_z` (LU + dgetri / zgetri
inverse) and `true_kappa1_posdef` / `true_kappa1_posdef_z` (Cholesky +
dpotri / zpotri).
`tests/scalapack/auxiliary/test_pdgecon.f90` and `test_pdpocon.f90`
factor the matrix on the distributed side, call the estimator, and
assert `rel_err_scalar(kappa_est, kappa_true) <= 0.7` — covers the
worst-case 1/3 LACON ratio with margin for fp noise.

Results (kind16, 2×2 grid):
- `test_pdgecon`: 3/3 PASS, κ_est ≈ 0.78 · κ_true (~22% rel-err)
- `test_pdpocon`: 3/3 PASS, κ_est ≈ 0.83 · κ_true (~17% rel-err)

## 2026-04-30 — Banded test sources (workspace / IPIV size fixes)

13 banded-family tests went from crash/FAIL to PASS on impi after three
test-side fixes:

1. **PDGBTRF / PZGBTRF IPIV under-allocation** — docs claim IPIV size
   `>= DESCA(NB)`, but the divide-and-conquer merge step at
   `pdgbtrf.f:1045` writes `IPIV(LN+1 .. LN+BM+BMN)` with
   `LN <= NB-BW` and `BM+BMN` up to `2*(BW+BWU)`. Tests now allocate
   `nb + 2*(bwl+bwu)`. Affects `test_p[dz]gbtrf`, `test_p[dz]gbtrs`,
   `test_p[dz]gbsv`.
2. **PDDBTRS / PZDBTRS workspace under-allocation** — upstream
   `WORK_SIZE_MIN = max(bwl,bwu)*nrhs` is too small; the local-update
   `QLAMOV` at `pqdbtrsv.f:1500` writes `BWU` rows × `NRHS` cols at
   stride `MAX_BW+BWL` starting from `WORK(1+MAX_BW-BWU)`. Tests now
   override `lwork` to `nrhs*max(bwl,bwu) + (nrhs-1)*(bwl+bwu)`.
   Affects `test_p[dz]dbtrf`, `test_p[dz]dbtrs`, `test_p[dz]dbtrsv`.
3. **test_pdpttrf invalid reference comparison** — original test
   compared `d_loc, e_loc` against LAPACK `dpttrf`. PDPTTRF's modified
   d/e are NOT LDL^T factors; multi-rank cyclic reduction stores
   reduced-system fillin in `af`, leaving d/e in a representation only
   PDPTTRS understands. Test now does factor+solve and compares solve
   result, mirroring `test_pzpttrf` (which was already correct).

## 2026-04-30 — Banded solver families delivered

Tridiagonal piece (`pddtsv` / `pdptsv`) plus all 22 banded drivers
(`p[dz]gb{sv,trf,trs}`, `p[dz]pb{sv,trf,trs,trsv}`,
`p[dz]db{sv,trf,trs,trsv}`) landed in phases 45/46/55/56. Each driver
inlines its own packing — no shared `banded_pack` module materialized;
the per-driver loops are 3 lines and vary in `bw` vs `bwl/bwu` vs
`uplo`, so consolidation has no payoff. Diagonal row in PDGBTRF
storage: `BWL+2*BWU+1` (verified in `test_pdgbsv.f90`); PDDBTRF /
PDPBTRF use the standard LAPACK packing (diagonal at `BWU+1` / `BW+1`).

## 2026-04-30 — pdsyevx / pzheevx fixed

Root cause was upstream PJLAENV's strict S/D/C/Z precision-letter gate
(`external/scalapack-2.2.3/SRC/pjlaenv.f:201-204`) that returns the
function result *uninitialised* when the routine name doesn't start
with one of those letters. After migration, callers pass names like
'PESYTTRD' / 'PXHETTRD' / 'PMSYTTRD' — the gate trips, ANB picks up a
stack-resident garbage integer, NSYTRD_LWOPT overflows to INT32_MIN,
the caller allocates a tiny WORK and the next call scribbles past it.

Fixed via `PJLAENV_EP` (`recipes/scalapack/extras/pjlaenv_ep.f`) +
`extra_renames PJLAENV → PJLAENV_EP` in `scalapack.yaml`. The same
upstream bug affects LAPACK's ILAENV / ILAENV2STAGE / IPARAM2STAGE;
those got the parallel `ILAENV_EP` / `ILAENV2STAGE_EP` /
`IPARAM2STAGE_EP` treatment in `lapack.yaml`.

A second (independent) bug surfaced after PJLAENV_EP unblocked the
workspace query: `target_pdsyevx` and `target_pzheevx` allocated
`WORK_T(1)` (or `RWORK_T(1)`) for the LQUERY call, but upstream
P{D,Z}SYEVX writes WORK(1:3) on rank 0 (ABSTOL/VL/VU broadcast) before
reporting the optimal size. AddressSanitizer caught it; bumped to
(R)WORK_T(3) in the wrapper.

Test drivers: `tests/scalapack/eigenvalue/test_pdsyevx.f90` and
`test_pzheevx.f90` (RANGE='A', JOBZ='N'). Both PASS on kind10 /
kind16 / multifloats to ~18 digits. pdsygvx / pzhegvx (the generalized
eigensolver siblings) sit on the same code path and are unblocked by
the same fix; drivers TBD.

## 2026-04-30 — pdlanhs / pxlanhs fixed

Upstream ScaLAPACK 2.2.3 bug in `pdlanhs.f` / `pzlanhs.f`: the NPROW=1
path failed to advance `II` after the first column block, so the
inner-loop bound `MIN(II+LL-JJ+1, IIA+NP-1)` collapsed to row 2 for
every column past the first MB. M-norm passed by luck (max element
typically lay in the kept range); 1/F/I norms underestimated by 10-20%.
Fixed via `source_overrides` in `recipes/scalapack.yaml` — the patched
upstream bodies live in `recipes/scalapack/source_overrides/` and go
through the normal migration pipeline.

Tests: `tests/scalapack/auxiliary/test_pdlanhs.f90` /
`test_pzlanhs.f90` both PASS to full target precision on all 3
targets.

## 2026-04-30 — pxhetrd / pxheev fixed

Originally these returned garbage (rel err ~3-5) on the Hermitian side
while their real-symmetric counterparts (pdsytrd / pdsyev) passed
cleanly. Root cause was the precision-independent PBLAS dispatcher
files (PB_C* helpers) shipping with stale `(double*)` pointer-cast
strides on KIND targets — fixed by the std/extension archive split,
which routes the alias-widened dispatcher through the migrated archive
while the standard double-precision body lives separately. Both
`test_pzhetrd` (matrix-element residual) and `test_pzheev` (eigenvalue
+ residual on JOBZ='V') now pass to full target precision (32+ digits
on kind16).

## 2026-04-30 — Real-MPI exercise via `linux-impi` preset (kind16)

`cmake --preset=linux-impi` invokes Intel MPI 2021.18 with
`I_MPI_ADJUST_REDUCE=1` and `I_MPI_ADJUST_ALLREDUCE=1` injected via
`MPIEXEC_PREFLAGS` to skip the impi shortpaths that mishandle
BLACS-registered user MPI ops over `REAL(KIND=16)`. ScaLAPACK tests
went from `0% real signal` (singleton MPI) to **170/170 (100%)** on
the kind16 target with a real 2×2 grid.

The 13 ScaLAPACK tests that originally failed under impi were resolved
by 4 distinct fixes (commits `cfcf023`, `67c57d2`, `e469670`):

* **norms (pdlanhs/pzlanhs)** — upstream IAROW double-advance bug;
  patched in source_overrides. See `doc/upstream-bugs/README.md`
  "p?lanhs.f IAROW double-advance".
* **equiv (pzgeequ)** — upstream wrong-axis reduction bug; patched in
  source_overrides. See `doc/upstream-bugs/README.md`
  "p?geequ.f column-scale reduction wrong axis". (pdgeequ was passing
  only by coincidence; the same fix applies.)
* **QR-fac (pdggrqf/pzggrqf)** — wrapper TAUB used the row axis where
  TAUB is column-distributed; fixed `loccols(descb)` for TAUB and
  bumped the test pre-allocation to `locn_b`.
* **refine (rfs ×6) + expert (svx ×2)** — the impi preset was missing
  `I_MPI_ADJUST_ALLREDUCE=1`; BLACS *GAMX2D / IGAMX2D paths through
  `MPIR_Reduce_local` crashed in the intranode shortpath on user
  MPI_Op + REAL(KIND=16) payload.
