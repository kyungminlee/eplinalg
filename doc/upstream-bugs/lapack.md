# Upstream bugs: LAPACK 3.12.1

*Last catalogued: 2026-05-11. See [`README.md`](README.md) for the cross-library
index, audit methodology, bug summary table, and how fixes are carried.*

This file collects LAPACK 3.12.1 bugs in the vendored
`external/lapack-3.12.1/` source. ScaLAPACK and MUMPS bugs are
catalogued in [`scalapack.md`](scalapack.md) and
[`mumps.md`](mumps.md) respectively. Bugs that span both LAPACK and
ScaLAPACK (e.g. the XERBLA-string sweep) live in the index file.

## 2026-05-11 LAPACK residual-divergence audit

After the symmetric-fix sweep below cleared the asymmetric-patch
backlog, an audit of the remaining 122 raw S↔D / C↔Z divergences
surfaced 13 more upstream LAPACK bugs (and one ScaLAPACK / one
migrator-side fix). All patches live under `recipes/lapack/patches/`
unless noted otherwise.

Numerical or correctness-affecting bugs:

| File | Bug | Sibling state |
|---|---|---|
| `clahef.f` | IPIV-undo loop boundary off-by-one in **both** branches: upper uses `J.LE.N` where `slasyf` / `dlasyf` / `zlahef` use `J.LT.N`; lower uses `J.GE.1` where the family uses `J.GT.1`. One extra iteration of the swap-back loop, redundant `CSWAP` at edge | clahef alone — fixed via `clahef.f.patch` |
| `slaqr5.f` | Line 604 typo: `T3 = T2*VT(3)` introduces a stray `VT(2)` factor (since `T2 = T1*VT(2)`). Every other `T3 = ...` site in slaqr5 (5 of them) and dlaqr5 all use `T3 = T1*VT(3)`. Affects the bulge-start tolerance test at line 608, so the bug can change bulge-start decisions in the single-precision QR iteration | slaqr5 alone — fixed via `slaqr5.f.patch` |
| `cgesvj.f` | Lines 512-513: `BIG = ONE / SFMIN` is active while `BIG = SLAMCH('Overflow')` is commented out. The other three siblings (sgesvj, dgesvj, zgesvj) all use the LAMCH form (the canonical LAPACK overflow constant). At edge values the two forms differ | cgesvj alone — fixed via `cgesvj.f.patch` |
| `sgeesx.f` | Line 648: `IWORK(1) = SDIM*(N-SDIM)` returns 0 when `SDIM=0` or `SDIM=N`. dgeesx wraps in `MAX(1, ...)` per LAPACK API convention (optimal-workspace report must be ≥1) | sgeesx alone — fixed via `sgeesx.f.patch` |
| `zgemlq.f` | Lines 261 and 286: `WORK(1) = LW` where the three siblings (cgemlq, sgemlq, dgemlq) all return `LWMIN`. `LWMIN = MAX(1, LW)`; zgemlq reports 0 instead of 1 for empty input | zgemlq alone — fixed via `zgemlq.f.patch` |
| `stgex2.f` | Lines 289 and 291: `LWORK.LT.MAX(N*M, M*M*2)` and `WORK(1) = REAL(MAX(N*M, M*M*2))` missing the LAPACK-conventional `MAX(1, ...)` floor. dtgex2 has the guard | stgex2 alone — fixed via `stgex2.f.patch` |
| `cheevx.f` | Line 360: `LWKOPT = (NB+1)*N` missing the `MAX(1, ...)` guard that zheevx has at the same site. Returns LWKOPT=0 when N=0 | cheevx alone — fixed via `cheevx.f.patch` |

EXTERNAL-list-vs-CALL-site bugs (declared routine missing, or wrong name):

| File | Bug | Sibling state |
|---|---|---|
| `clahef_rk.f` | Line 302: EXTERNAL declares `CGEMM` but the routine actually `CALLs CGEMMTR` at lines 760 and 1190. slasyf_rk / dlasyf_rk / zlahef_rk all correctly list `GEMMTR` | clahef_rk alone — fixed via `clahef_rk.f.patch` |
| `dopmtr.f` | Line 180: EXTERNAL lists only `DLARF` but the routine `CALLs DLARF1L` at line 264 (undeclared) and `CALL DLARF` at line 320 wrapped in manual `AII` save/restore (the inline expansion of what `DLARF1F` does internally). sopmtr / cunmr2 / zunmr2 all use the modern `LARF1F/LARF1L` wrappers | dopmtr alone — fixed via `dopmtr.f.patch` (also switches the line-320 call to `DLARF1F` and drops the manual `AII` tweak) |
| `dormr2.f` | Lines 188 + 265-269: declares `DLARF` and does manual `AII` save/restore around the call. sormr2 / cunmr2 / zunmr2 all use `DLARF1L` (the modern wrapper handles the `AII` trick internally) | dormr2 alone — fixed via `dormr2.f.patch` |
| `sggev3.f` | Lines 262-265: EXTERNAL omits `XERBLA`, but `CALL XERBLA('SGGEV3 ', -INFO)` at line 359. dggev3 declares it | sggev3 alone — fixed via `sggev3.f.patch` |
| `sgges.f` | Lines 325-327: EXTERNAL omits `XERBLA`, but `CALL XERBLA('SGGES ', -INFO)` at line 420. dgges declares it correctly | sgges alone — fixed via `sgges.f.patch` |
| `slaqr2.f` | Line 315: EXTERNAL declares `SLARF1L`, but the routine actually `CALLs SLARF1F` at lines 604, 606, 608. `SLARF1L` is dead, `SLARF1F` is undeclared. dlaqr2 declares the right one; slaqr3 / dlaqr3 were correct upstream | slaqr2 alone — fixed via `slaqr2.f.patch` |

### Severity

All 13 are real correctness or interface bugs. Most have visible
impact only at edge cases (N=0, SDIM=0, empty-matrix workspace
queries) or with vendor LAPACK that rounds workspace differently
from Netlib. The `clahef` boundary fixes, `slaqr5` typo, and `cgesvj`
overflow-constant difference are the only ones that could affect
numerics in a non-edge-case path.

### Upstream report

Not yet filed.

### Comparer-side fixes landed alongside

`_canonicalize_for_compare` extended with 5 additional W-class
normalizations (per `doc/archive/lapack-residual-divergence-categorization.md`):
`ELSE IF` ↔ `ELSEIF`, `GO TO` ↔ `GOTO`, `DOUBLE PRECISION` ↔ `REAL`,
`COMPLEX*16` ↔ `COMPLEX`, `''` / `' '` / `'  '` CHARACTER-literal
collapse, XERBLA trailing-space strip, `::` whitespace strip,
`)\s*THEN` collapse. Plus label canonicalization (`DO 100` ↔ `DO 110`
re-numbering) moved into the diverge comparator.

### Migrator-side fix landed alongside

`_strip_roundup_lwork` (`src/migrator/fortran_migrator.py:3300`)
extended to recognise the F90 attribute-list form
``REAL, EXTERNAL :: SROUNDUP_LWORK``. The existing F77 form
``REAL ... + EXTERNAL ...`` on separate lines was already handled;
the modern form (slaqz0 / slaqz3 / slaqz4) leaked its orphan
declaration into the migrated output.

### B? candidates after deep-dive (2026-05-11 follow-up)

After per-pair upstream-source review the seven B? candidates resolve as:

**Reclassified to W (cosmetic / documented-design, not bugs):**

- `slasq2.f:279` — `IEEE = .FALSE.` is upstream-documented at line 276
  (`[11/15/2008]` comment) as a deliberate single-precision-only
  workaround for type-16 test matrices. Not a bug; the limitation is
  precision-specific by design.
- `sorcsd.f` — passing `DUMMY(1)` for A/TAU during LWORK=-1 workspace
  query is equivalent to passing the actual arrays (`U1`, `THETA`):
  LAPACK does not read those slots in query mode. Cosmetic idiom split.
- `sgedmdq.f90` — original "missing `.OR.WNTVCQ` branch" claim was a
  mis-read; both halves use identical `IF (WNTVEC .OR. WNTVCF)` test.
  Only divergence is literal `MIN(M,N)` vs the local `MINMN` (same
  value). Cosmetic.

**Real upstream inconsistency, no migrator impact:**

- `ssytri2.f:166` queries `ILAENV` for `'SSYTRF'`'s blocksize where
  `dsytri2.f:165` queries for `'DSYTRI2'`. Across the six-routine
  family: ssytri2 → SSYTRF, dsytri2 → self, csytri2/zsytri2 → self,
  chetri2/zhetri2 → HETRF. Genuinely inconsistent. The returned
  blocksize feeds only the MINSIZE workspace-hint computation; no
  numerical effect. Worth a Netlib report. Not patched.

**Genuine B? remaining — domain-expertise gated, in the GEDMD/GEJSV
family (newer LAPACK, known asymmetric maintenance):**

- `sgejsv.f:1711` passes `'L'` to `SGESVJ` where `dgejsv.f:1711` passes
  `'G'` (JOBA: Lower-triangular hint vs General). Operand is `U` from
  prior operations; requires tracing whether U is L-structured at that
  call site.
- `cgedmd.f90` line 704/916 GEJSV JOBR `'N'` vs `'R'` in zgedmd; line
  751 OFL conservative `SLAMCH('O')*SLAMCH('P')` (cgedmd) vs raw
  `DLAMCH('O')` (zgedmd). C-half looks more carefully written.
- `sgedmd.f90` line 716 passes `-1` to GESVDQ's LIWORK (consistent
  with query mode); dgedmd passes actual `LIWORK`. Line 724 sgedmd
  has `LWRSVQ = INT(RDUMMY(1))` where dgedmd has
  `LWRSVQ = MAX(MWRSVQ, INT(RDUMMY(1)))` — D has a safety floor S
  lacks.

> **Reporting note.** These three items, plus the `ssytri2` ILAENV
> inconsistency above, should be filed **as issues** on the Netlib
> LAPACK tracker (https://github.com/Reference-LAPACK/lapack/issues),
> *not* as PRs. We cannot determine the correct resolution from
> reading the source alone — picking one branch over the other for
> a PR would be a guess. The maintainers / GEDMD-GEJSV authors have
> the algorithmic context to decide which half is canonical.

Net of the deep-dive: zero new migrator patches; one upstream-only
report candidate (ssytri2 ILAENV typo).

See `doc/archive/lapack-residual-divergence-categorization.md` for the full
110-pair classification and the deep-dive table.


## 2026-05-12 B? deep-dive (LAPACK portion)

The B? items remaining in the `2026-05-11` audit were re-examined; the
LAPACK-side findings are below. See [`scalapack.md`](scalapack.md) for
the ScaLAPACK portion (`pslaqr3` / `pdlaqr3` LWK8 divergence).

### Promoted to B: `cgedmdq.f90` — three workspace-query bugs

**Symptom.** `CGEDMDQ` is the Q-projected variant of the DMD driver.
Three independent bugs in its workspace-query path, all absent from
`ZGEDMDQ`:

1. **Line 638: LQUERY missing LZWORK==-1.** The query predicate
   recognizes only `LWORK==-1` or `LIWORK==-1`, not `LZWORK==-1`.
   A caller asking only for the complex-workspace size by passing
   `LZWORK=-1` and positive LWORK/LIWORK hits the validation block
   instead of the query block, and returns INFO=-31 (or similar
   negative INFO) instead of populating ZWORK(1)/(2). Same shape as
   the `pdsyevd.f` and `pzheevd.f` LIWORK gaps documented above —
   this is now the **fourth occurrence of this exact bug shape** in
   the LAPACK/ScaLAPACK family.

2. **Lines 700-701: missing ZWORK(1)/(2) init in the void-input
   quick-return.** When `N==0` or `N==1` and `LQUERY` is true, the
   routine sets `IWORK(1)=1`, `WORK(1)=2`, `WORK(2)=2` and returns.
   It does NOT set `ZWORK(1)/(2)`, so the caller reads stale memory
   for the complex-workspace size. ZGEDMDQ at the analogous line
   (699-702) sets all five (`IWORK`, `ZWORK(1)/(2)`, `WORK(1)/(2)`).

3. **Line 725: inner CGEDMD workspace query passes wrong arg
   shape.** During the LQUERY path, CGEDMDQ calls inner CGEDMD to
   discover its workspace requirements. The c-half call:

   ```fortran
   CALL CGEDMD( ..., S, LDS, ZWORK, LZWORK, WORK, -1, IWORK,
                LIWORK, INFO1 )
   ```
   passes the user-supplied LZWORK and LIWORK rather than `-1`. So
   CGEDMD only enters query mode for LWORK and only populates
   `WORK(1)`. The subsequent reads `MLWDMD = INT(ZWORK(1))` and
   `IMINWR = MAX(IMINWR, IWORK(1))` see stale/garbage values, so
   CGEDMDQ's reported optimal LZWORK and LIWORK are wrong.

   ZGEDMDQ at line 726 passes `-1` to all three workspace args:

   ```fortran
   CALL ZGEDMD( ..., S, LDS, ZWORK, -1, WORK, -1, IWORK,
                -1, INFO1 )
   ```

**Affected files.**
- `external/lapack-3.12.1/SRC/cgedmdq.f90` line 638 (LQUERY); lines
  700-701 (void-input quick-return missing ZWORK init); line 725
  (inner CGEDMD query call).

**Fix.** Patch all three from the Z-half reference. See
`recipes/lapack/patches/cgedmdq.f90.patch`.

**Severity.** Workspace-query API contract violation. Callers
following the LAPACK convention (either `LZWORK=-1` for a complex-
workspace-only query, OR full `LWORK=LZWORK=LIWORK=-1`) get wrong
results. Migrated `qgedmdq` is unaffected (Z is canonical), but the
upstream `cgedmdq.f90` is broken for direct single-complex users.

**Why upstream's tests miss it.** Same as the pzheevd / pcheevd /
pdsyevd LIWORK family: workspace queries by individual workspace
arg (rather than the bundled all-args query) are uncommon in
upstream test drivers. The ZGEDMDQ branch is exercised; CGEDMDQ
likely went through a single mechanical S→C find-replace pass that
preserved the (broken) shape from an earlier draft.

**Upstream report.** File as a PR at
https://github.com/Reference-LAPACK/lapack — the fix is mechanical
(copy from ZGEDMDQ).

### Stays B?: `cgedmd` vs `zgedmd` — two independent inconsistencies

**Symptom.** Two distinct cross-half asymmetries in the GEDMD
driver, each with a 3-vs-1 vote shape (three halves agree, one is
the outlier — but the *outlier is a different file* in each case):

| Asymmetry | sgedmd | dgedmd | cgedmd | zgedmd |
|---|---|---|---|---|
| Inner `?GEJSV` JOBR arg (positional 4 in the call) | `'N'` | `'N'` | `'N'` | **`'R'`** |
| `OFL` overflow threshold formula | `SLAMCH('O')` | `DLAMCH('O')` | **`SLAMCH('O')*SLAMCH('P')`** | `DLAMCH('O')` |

`zgedmd` is the outlier on JOBR; `cgedmd` is the outlier on OFL.

**Why each matters.** JOBR `'R'` tells GEJSV to flag and discard
singular values smaller than the overflow threshold; `'N'` keeps
them all. Different filtering strictness. The OFL formula sets the
threshold used in `IF (SCALE .GE. (OFL / ROOTSC))` at lines 774 and
848 — a scaling-safety guard. `LAMCH('O')*LAMCH('P')` is much
smaller than `LAMCH('O')` (multiplying by precision epsilon),
making cgedmd's guard fire on a much larger range of inputs.

Both differences are real algorithmic choices. Each could be (a) a
fix that landed on one half and wasn't backported, or (b) a typo
that escaped review. Without DMD-domain expertise we cannot tell
which half is correct in either case.

**Affected files.**
- `external/lapack-3.12.1/SRC/cgedmd.f90:703, 916` (JOBR `'N'`).
- `external/lapack-3.12.1/SRC/zgedmd.f90:704, 916` (JOBR `'R'`).
- `external/lapack-3.12.1/SRC/cgedmd.f90:751` (OFL conservative).
- `external/lapack-3.12.1/SRC/zgedmd.f90:751` (OFL raw).
- `external/lapack-3.12.1/SRC/sgedmd.f90:772` (OFL raw, JOBR `'N'`).
- `external/lapack-3.12.1/SRC/dgedmd.f90:772` (OFL raw, JOBR `'N'`).

**Severity.** Numerical accuracy / robustness, not memory safety
or API contract. JOBR `'R'` vs `'N'` changes which singular values
GEJSV reports — for ill-conditioned inputs the two halves can
return different DMD spectra. OFL formula difference changes the
threshold at which cgedmd's scaling-safety guard fires; for inputs
near the overflow boundary, cgedmd takes the safe-scaling branch
while sgedmd/dgedmd/zgedmd take the fast-scaling branch. Whether
the safe branch produces a *more* or *less* accurate result is
algorithm-dependent.

**Why upstream's tests miss them.** GEDMD landed in LAPACK 3.10
(2021) and the GEDMD test suite focuses on well-conditioned DMD
spectra. Neither asymmetry causes a test failure on the standard
matrices. The 3-vs-1 vote shape with the outlier on a different
file in each case strongly suggests these are uncoordinated
hand-edits made during separate review passes rather than a
deliberate per-precision design choice.

**Upstream report.** File as **two separate issues** (not PRs) on
the Netlib LAPACK tracker — one per asymmetry. Authors of the
GEDMD/GEJSV family (cited in the file headers as Z. Drmac et al.)
have the context to decide which form is canonical for each.


## LAPACK 3.12.1: `?orbdb3.f` `?ROT` stride uses LDX11 where it should be LDX21

**Symptom.** Migrated extended-precision `?orbdb3` runs cleanly (the
override carries the fix), but the kind4 / kind8 baseline ctest cycle
— which links the unmodified `external/lapack-3.12.1/SRC/` archive —
trips one of two heap-corruption manifestations. `lapack_test_zunbdb3`
SIGSEGVs in `arena_for_chunk → __libc_free` during the wrapper's
auto-deallocate. `lapack_test_dorbdb3` doesn't crash but reports
`digits ~= 1` against the quad reference (it presents as a precision
shortfall — the actual mechanism is corruption of an adjacent buffer
on the heap, presented as wrong values in `theta` / `phi`). The
`s` / `c` precision halves carry the byte-identical bug; would
surface under the kind4 baseline once that column is wired.

**Root cause.** Inside the `IF( I .GT. 1 )` branch of `?orbdb3.f`'s
main loop, the plane-rotation call passes `LDX11` as the `INCY` for
`X21`:

```fortran
DO I = 1, M-P
   IF( I .GT. 1 ) THEN
      CALL DROT( Q-I+1, X11(I-1,I), LDX11, X21(I,I), LDX11, C, S )
   END IF
   ...
END DO
```

`LDX11` is correct for the `X` operand (which lives in `X11`). But
`Y` is `X21` — its leading dimension is `LDX21`. The orbdb3 regime
requires `M-P <= min(P, Q, M-Q)`, so `M-P < P` and the two leading
dimensions disagree (typical: `LDX11 = P`, `LDX21 = M-P`). With
`INCY = LDX11 > LDX21`, the rotation strides `Q-I+1` elements through
`X21` at stride `LDX11`, walking off the end of `X21`'s allocated
buffer and overwriting whatever heap chunk follows.

Layout decides whether the corruption surfaces as a SIGSEGV (when
the trampled bytes are chunk metadata) or as silently wrong numerics
(when they're data words of an adjacent allocation). Both
manifestations are heap-state-dependent — the same input shape can
flip between modes across builds.

**Why M-norm-style accidental passes don't apply.** Unlike the
ScaLAPACK lanhs-style "missed elements happen to be zero" stories
below, every test invocation that reaches `I > 1` triggers the
out-of-bounds write. The bug is data-independent; what varies is
only whether the heap corruption surfaces immediately or silently.

**Affected files.**

* `external/lapack-3.12.1/SRC/dorbdb3.f` (used by our migrated D-half).
* `external/lapack-3.12.1/SRC/zorbdb3.f` (used by our migrated Z-half).
* `external/lapack-3.12.1/SRC/sorbdb3.f` (S-half — same bug, exercised
  only under the kind4 baseline column).
* `external/lapack-3.12.1/SRC/corbdb3.f` (C-half — same status).

**Fix.** One-character change: `LDX11` → `LDX21` on the `?ROT` `INCY`
argument inside the `I .GT. 1` branch. Carried in
`recipes/lapack/patches/{d,s,z,c}orbdb3.f` and wired via
`recipes/lapack.yaml`'s `patches:` list. Every migrated target
(kind10, kind16, multifloats) builds against the patched form. The
`prefer_source` pin isn't needed here — LAPACK's canonical-rank picker
prefers D over the other halves by default, and the override matches
that choice.

**Standard-precision archive still buggy.** The kind4 / kind8 baseline
path (`migrator stage --target kind{4,8}`) skips migration and stages
upstream `external/lapack-3.12.1/SRC/` verbatim into `_reflapack_src/`,
so the std `lapack` archive that the baseline links against still has
the typo. Two of the residual failures in
`doc/archive/kind48-baseline-status-20260506.md` (`zunbdb3`, `dorbdb3`)
trace to this
gap and are parked under `tests/lapack/TODO.md` until the staging path
overlays `recipes/lapack/patches/` on top of the upstream
copy. Confirmed both clear when the override is overlaid.

**Why upstream's tests miss it.** Reference LAPACK's `?orbdb3` test
shapes happen to use `LDX11 = LDX21` (the bug is silent when the two
strides match). The mismatch surfaces only when a caller invokes
`?orbdb3` with `LDX11 /= LDX21`, which the upstream test driver
doesn't exercise but a real CSD-driver call path does.

**Upstream report.** Not yet filed.

---

## LAPACK 3.12.1: `dgejsv.f` workspace size argument off by N to inner DGESVJ

**Symptom.** Silent stack/workspace overflow. In the
`RSVEC && !LSVEC && !ALMORT` branch of DGEJSV (the "two-more-QR-
factorizations" path used when neither left singular vectors nor an
"almost-orthogonal" shortcut applies), DGEJSV calls DGESVJ with the
workspace base offset by N elements but passes the *full* `LWORK` as
the workspace-size argument. DGESVJ may then write up to N entries
past the end of `WORK(LWORK)`. The float sibling `sgejsv.f` correctly
passes `LWORK-N`.

**Root cause.** `dgejsv.f:1174-1175`:

```fortran
CALL DGESVJ( 'Lower', 'U','N', NR, NR, V,LDV, SVA, NR, U,
$            LDU, WORK(N+1), LWORK, INFO )
```

The 13th argument is `WORK(N+1)` (workspace base, N elements in)
but the 14th argument (workspace size) is `LWORK`, not `LWORK-N`.
DGESVJ trusts its 14th arg as the upper bound on writes. Compare
`sgejsv.f:1174-1175`:

```fortran
CALL SGESVJ( 'Lower', 'U','N', NR, NR, V,LDV, SVA, NR, U,
$            LDU, WORK(N+1), LWORK-N, INFO )
```

The same routine on the float side passes `LWORK-N`; the bug is in
the D copy only.

**Affected files.**
- `external/lapack-3.12.1/SRC/dgejsv.f` (line 1175, 14th arg of the
  inner DGESVJ call).

**Fix.** Single-token change `LWORK → LWORK-N` on the 14th argument.
Carried in `recipes/lapack/patches/dgejsv.f`. Wired via
`recipes/lapack.yaml`'s `patches:` list. (No `prefer_source`
pin needed: D-half is canonical for the migrated archive by default,
and the S-half value is already correct, so converge folds.)

**Why upstream's tests miss it.** The DGEJSV test driver
(`TESTING/EIG/dchksvj.f`) calls DGEJSV with workspace generously
sized — typically allocated as `LWORK = N*N + 5*N` or larger. The
N-element overrun lands inside the spare margin and produces no
detectable corruption on the test hosts. The bug surfaces when a
caller follows the documented minimum-LWORK formula exactly: the
overrun then hits whatever lives immediately after the WORK array on
the stack.

**Upstream report.** Not yet filed.

---

## LAPACK 3.12.1: `?trsyl3.f` D/C/Z halves miss LIWORK / LDSWORK argument-validation

**Symptom.** `STRSYL3` reports `INFO=-14` when the caller passes a
too-small `LIWORK` and `INFO=-16` when `LDSWORK` is too small (when
not querying). The D/C/Z halves omit those guards: a too-small
`LIWORK` or `LDSWORK` slips past argument validation without an
`XERBLA` callback. Functionally not catastrophic — the routine has a
downstream "use unblocked code" fallback (`MIN(NBA,NBB).EQ.1 .OR.
LDSWORK.LT.MAX(NBA,NBB) .OR. LIWORK.LT.IWORK(1)`) that drops to
`?TRSYL` when workspace is insufficient — but the missing guard
masks user errors, exactly the same class of bug as the
`pssyevd.f` LQUERY gap below.

**Root cause.** `dtrsyl3.f:286-288` ends the validation block at
`LDC` without checking `LIWORK` or `LDSWORK`. Same omission in
`ctrsyl3.f` and `ztrsyl3.f`. `strsyl3.f:286-290` correctly extends:

```fortran
ELSE IF( .NOT.LQUERY .AND. LIWORK.LT.IWORK(1) ) THEN
   INFO = -14
ELSE IF( .NOT.LQUERY .AND. LDSWORK.LT.MAX( NBA, NBB ) ) THEN
   INFO = -16
END IF
```

The `IWORK(1)` reference is well-defined here because the LQUERY-prep
block earlier in the routine writes `IWORK(1) = NBA + NBB + 2`
unconditionally before the test block runs. (Note: `ZTRSYL3` /
`CTRSYL3` have no `IWORK` argument — only the `LDSWORK` guard
applies, and `LDSWORK` is at parameter position 14, not 16.)

**Affected files.**
- `external/lapack-3.12.1/SRC/dtrsyl3.f` (after line 287, before
  `END IF` at line 288).
- `external/lapack-3.12.1/SRC/ztrsyl3.f` (after line 254, before
  `END IF` at line 256).
- `external/lapack-3.12.1/SRC/ctrsyl3.f` (analogous; no override
  carried — C half is not canonical for the migrated archive).

**Fix.** Insert the two `ELSE IF` branches into D's validation block;
insert just the `LDSWORK` branch (at parameter position 14) into Z's.
Carried in `recipes/lapack/patches/dtrsyl3.f` and
`ztrsyl3.f`.

**Why upstream's tests miss it.** Reference LAPACK's `dchktz` /
`zchktz` test drivers always size `IWORK` and `SWORK` by querying
first, so the too-small-workspace path on the validation side is
never exercised.

**Upstream report.** Not yet filed.

---

## LAPACK 3.12.1: stale / missing EXTERNAL declarations in 5+ routines

**Symptom.** None at runtime — Fortran `EXTERNAL` is advisory and the
linker resolves the actual call sites correctly via the global symbol
table. The typos surface only as divergences against the matching
S/C precision sibling in the migrator's convergence reports. Same
flavor as the `zupmtr.f` `ZLARF1` / `ZLARF1L` typo (already patched
at the top of this file) and the ScaLAPACK `pzunmbr.f`
`PCHK1MAT` / `PCHK2MAT` typo.

**Affected files / fixes.** All of these are single-line edits to the
`EXTERNAL` declaration block; the body is correct in each case.

| File | Issue | Fix |
|------|-------|-----|
| `dlaqp2.f:174` | declares `DLARF`, body calls `DLARF1F` (line 222). slaqp2 correctly declares `SLARF1F`. | replace `DLARF` with `DLARF1F` |
| `zlarf1f.f:185` | missing `ZAXPY`. Body calls ZAXPY at lines 287, 293. clarf1f also missing CAXPY (not patched — C-half is non-canonical). | add `ZAXPY` |
| `zhetrf_aa.f:168` | declares stale `ZGEMV` — body never calls ZGEMV (only ZGEMM). chetrf_aa correctly omits CGEMV. | drop `ZGEMV` |
| `zlahef_aa.f:175` | declares stale `ZGEMM` — body never calls ZGEMM (only ZGEMV). clahef_aa correctly omits CGEMM. | drop `ZGEMM` |
| `zgbrfsx.f:499` | missing `ILATRANS`. Body calls ILATRANS at lines 511 and 561 (`TRANS_TYPE = ILATRANS(TRANS)` and the `IF(TRANS_TYPE.EQ.-1)` test). The other three rfsx siblings (`sgbrfsx`, `dgbrfsx`, `cgbrfsx`) all declare ILATRANS in their EXTERNAL list. | add `ILATRANS` |

**Carried in:** `recipes/lapack/patches/{dlaqp2,zlarf1f,zhetrf_aa,zlahef_aa,zgbrfsx}.f`,
wired in `recipes/lapack.yaml`.

**Other instances surfaced by the same audit, NOT patched** because
the migrator picks D/Z as canonical and the sibling S/C bug never
rides into the migrated archive:

- `clarf1f.f:154`: declares `CGER` (dead — body calls `CGERC`); also
  missing `CAXPY`. Real upstream typo, but only the standalone single-
  precision archive is affected, not the migrated extended-precision
  one.
- `slaqr2.f:313`: declares `SLARF1L` (dead — body calls `SLARF1F`).
- `sgelqt.f:143`: declares stale `SGEQRT2` and `SGEQRT3` — body only
  calls `SGELQT3` and `SLARFB`.
- `ssysv_aa.f:227`: `XERBLA('SSYSV_AA', ...)` missing the trailing
  space that D/C/Z all use (`'?SYSV_AA '`). Cosmetic.
- `sgges.f:325-327`: declares `EXTERNAL SGEQRF, ..., STGSEN` —
  missing `XERBLA`. Body calls XERBLA at line 420. dgges correctly
  appends `XERBLA` to the EXTERNAL list.
- `sormr2.f:131` / `sopmtr.f`: declares `SLARF1L` (dead — body calls
  `SLARF` after manual `AII = A(...); A(...) = ONE`). dormr2 / dopmtr
  use `DLARF` directly, with the same AII trick. The S halves carry
  a stale `SLARF1L` that doesn't match the body.

**Why upstream's tests miss them.** Fortran's `EXTERNAL` is a
*declaration*, not an *invocation*. The linker resolves calls by
symbol name regardless of what `EXTERNAL` advertises, and the
upstream test drivers don't check declarations against call sites.
None of these has any observable runtime effect.

**Upstream report.** Not yet filed.

---

## LAPACK 3.12.1: INFO=-N parameter-position drift (mechanical sweep 2026-05-08)

A mechanical sweep parsed every `IF (test) THEN; INFO = -N` arg-
validation block in LAPACK SRC, identified the single signature arg
referenced in the test, and verified the literal `N` matches the
actual position of that arg in the SUBROUTINE signature. **128
mismatches out of 11 352 sites; after filtering the multi-arg-test
false positives and the LAPACK `INFO = -100*PARAM - INDEX` element-
validation convention (slasq2 et al. with `INFO = -201/202/203` for
`Z(1)/Z(2)/Z(3)` failures), 9 real bugs remain.**

### Patched in migrated archive (D/Z canonical)

| File | Test | Old INFO | Fixed INFO | Position of tested arg |
|------|------|---------:|-----------:|----------------------:|
| `dorbdb4.f` | `LWORK.LT.LWORKMIN` | -14 (WORK) | -15 | 15 (LWORK) |
| `zunbdb4.f` | `LWORK.LT.LWORKMIN` | -14 (WORK) | -15 | 15 (LWORK) |
| `dggsvd3.f` | `LWORK.LT.1` | -24 (INFO itself!) | -22 | 22 (LWORK) |
| `zggsvd3.f` | `LWORK.LT.1` | -24 (INFO itself!) | -22 | 22 (LWORK) |
| `dorcsd.f` | `LWORK.LT.LWORKMIN` | -22 (LDU2) | -28 | 28 (LWORK) |
| `zuncsd.f` | `LWORK.LT.LWORKMIN` | -22 (LDU2) | -28 | 28 (LWORK) |
| `zuncsd.f` | `LRWORK.LT.LRWORKMIN` | -24 (LDV1T) | -30 | 30 (LRWORK) |
| `zlaqz0.f` | `LWORK.LT.LWORKREQ` | -19 (RWORK) | -18 | 18 (LWORK) |
| `zlaqz2.f` | `LWORK.LT.LWORKREQ` | -26 (RWORK) | -25 | 25 (LWORK) |

### Documented but not patched (S/C non-canonical, same bugs)

- `sorbdb4.f` -14 (should be -15)
- `cunbdb4.f` -14 (should be -15)
- `sggsvd3.f`, `cggsvd3.f` -24 (should be -22)
- `sorcsd.f`, `cuncsd.f` -22 LWORK (should be -28); cuncsd -24 LRWORK (should be -30)
- `slaqz0.f`, `claqz0.f` -19 (should be -18)
- `claqz2.f` -26 (should be -25)

### Notable: dggsvd3 / zggsvd3 `INFO = -24` references INFO itself

The most striking bug in this batch: `?ggsvd3.f` on illegal LWORK
returns `INFO = -24` — but in those routines, position 24 *is* INFO.
So XERBLA reports "argument INFO had an illegal value", which is
nonsensical (the user can't pass an illegal INFO since it's an OUT
parameter).

### Severity

All diagnostic-only. The routines correctly reject the bad input;
they just print the wrong parameter number in the XERBLA message.
Same severity class as the PCHK?MAT param-position bugs and XERBLA-
string typos.

### Why upstream's tests miss them

Test drivers always pass valid LWORK (or query-then-allocate); the
INFO=-N path for too-small workspace never fires.

### `zrscl.f` doc-header references wrong routine

While auditing, also caught `zrscl.f:1` `\brief \b ZDRSCL` (and line
9, "Download ZDRSCL + dependencies") — the doc header refers to a
different routine, ZDRSCL, which zrscl.f *calls* internally (line
139). The actual routine in zrscl.f is ZRSCL. Auto-generated HTML
docs at netlib.org would label ZRSCL's docs page as "ZDRSCL" or fail
to find ZRSCL entirely. Fixed in
`recipes/lapack/patches/zrscl.f`.

**Upstream report.** Not yet filed.

---
## LAPACK 3.12.1: 8 more EXTERNAL-vs-CALL drifts (mechanical sweep 2026-05-08)

A second mechanical sweep — this one matching every CALL target in
each subprogram against the routine's EXTERNAL declaration — found
462 routines with both stale (declared, never called) and missing
(called, not declared) entries. After filtering keyword false
positives, 8 of these affect the migrated archive's D/Z-canonical
files:

| File | Bug |
|------|-----|
| `dlarf1l.f` | EXTERNAL declares `DGEMV, DGER`. Body also calls `DAXPY` (3×) and `DSCAL` (2×). Fix: add `DAXPY, DSCAL`. |
| `zlarf1l.f` | EXTERNAL declares `ZGEMV, ZGERC, ZSCAL`. Body also calls `ZAXPY` (2×). Fix: add `ZAXPY`. |
| `dlatrs3.f` | EXTERNAL omits `DGEMM`. Body calls `DGEMM` 2×. Fix: add `DGEMM`. |
| `zlatrs3.f` | EXTERNAL omits `ZGEMM`. Body calls `ZGEMM` 3×. Fix: add `ZGEMM`. |
| `zgelss.f` | EXTERNAL omits `ZUNMQR`. Body calls `ZUNMQR` 2×. Fix: add `ZUNMQR`. |
| `dorbdb4.f` | EXTERNAL declares `DLARF` (dead). Body calls `DLARF1F` 4×, never `DLARF`. Fix: replace `DLARF` with `DLARF1F`. |
| `dorgr2.f` | EXTERNAL declares `DLARF` (dead). Body calls `DLARF1L`. Fix: replace `DLARF` with `DLARF1L`. |
| `zrscl.f` | EXTERNAL omits `ZSCAL`. Body calls `ZSCAL` 4×. Fix: add `ZSCAL`. |

Carried as source overrides in `recipes/lapack/patches/`,
wired in `recipes/lapack.yaml`.

**Sibling halves with similar drift, NOT patched** (S/C non-canonical;
documenting for upstream report only):

- `slarf1l.f` missing `SAXPY`, `SSCAL`
- `clarf1l.f` missing `CAXPY`
- `slatrs3.f` missing `SGEMM`
- `clatrs3.f` missing `CGEMM`
- `cgelss.f` missing `CUNMQR`
- `sorbdb4.f` declares dead `SLARF` (calls `SLARF1F`)
- `sorgr2.f` declares dead `SLARF` (calls `SLARF1L`)
- `crscl.f` missing `CSCAL`
- `clarf1f.f` missing `CAXPY`, `CGERC` (already documented separately)
- `dopmtr.f` declares dead `DLARF`, calls `DLARF1L`

Plus four entries in `external/lapack-3.12.1/SRC/VARIANTS/lu/REC/` —
`{c,d,s,z}getrf.f` all missing `?GEMM`. The VARIANTS subdirectory is
not in the migrator's `source_dir` scan path, so these are out of
scope for the migrated archive but worth filing upstream.

**Severity.** All advisory — Fortran linker resolves CALL targets
by symbol name regardless of EXTERNAL. The fixes correct the
declaration to match the body so static analyzers, pretty-printers,
and downstream tooling see consistent metadata.

**Why upstream's tests miss them.** Reference test drivers compile
and link the routines unchanged from upstream — the linker resolves
the symbols, the tests pass. EXTERNAL-vs-CALL consistency is
checkable only by reading source.

**Upstream report.** Not yet filed.

---

## LAPACK 3.12.1: `dlaswlq.f` argument-validation accepts NB=0 where slaswlq rejects it

**Symptom.** `DLASWLQ` fails to reject `NB=0` as an illegal block-size
argument. The downstream block loop divides by `NB` and iterates as
`I=1, N-K, NB`; with `NB=0` this is a divide-by-zero (or infinite
loop, depending on compiler arithmetic).

**Root cause.** `dlaswlq.f:220-221`:

```fortran
ELSE IF( NB.LT.0 ) THEN
   INFO = -4
```

The float sibling `slaswlq.f:221` correctly tests `NB.LE.0`. NB=0
isn't a valid block size — block algorithms require positive blocks.
The single-precision validation is right; the double-precision
validation lost the `=` somewhere during a sync.

**Affected files.**
- `external/lapack-3.12.1/SRC/dlaswlq.f` (line 220).

**Fix.** `NB.LT.0 → NB.LE.0`. Carried in
`recipes/lapack/patches/dlaswlq.f`, wired in
`recipes/lapack.yaml`.

**Why upstream's tests miss it.** Reference test drivers always pass
NB > 0; the NB=0 path is only reachable via a malformed user call.

**Upstream report.** Not yet filed.

---

## LAPACK 3.12.1: `sggev3.f` workspace-query gates on ILVL where real call gates on ILV

**Symptom.** Workspace under-estimation when the caller requests
**right eigenvectors only** (`JOBVL='N', JOBVR='V'` → `ILVL=false,
ILVR=true, ILV=true`). The `LWORK=-1` query path returns a workspace
size computed for `SLAQZ0('E', ...)` (eigenvalues only), but the real
computation later issues `SLAQZ0('S', ...)` (full Schur), which needs
*more* workspace. A caller who follows the standard "query, allocate
exactly, then call" idiom runs out of workspace inside SLAQZ0 — likely
SIGSEGV or a workspace overrun, depending on what's adjacent on the
stack.

**Root cause.** `sggev3.f:339-352` (workspace-query block):

```fortran
IF( ILVL ) THEN                                      ! ← wrong predicate
   CALL SORGQR( N, N, N, VL, LDVL, WORK, WORK, -1, IERR )
   LWKOPT = MAX( LWKOPT, 3*N+INT( WORK( 1 ) ) )
   CALL SLAQZ0( 'S', JOBVL, JOBVR, N, 1, N, A, LDA, B, LDB,
$                ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR,
$                WORK, -1, 0, IERR )
   LWKOPT = MAX( LWKOPT, 2*N+INT( WORK( 1 ) ) )
ELSE
   CALL SLAQZ0( 'E', ..., WORK, -1, 0, IERR )
   LWKOPT = MAX( LWKOPT, 2*N+INT( WORK( 1 ) ) )
END IF
```

vs `sggev3.f:471-477` (real-call CHTEMP setup):

```fortran
IF( ILV ) THEN
   CHTEMP = 'S'
ELSE
   CHTEMP = 'E'
END IF
CALL SLAQZ0( CHTEMP, JOBVL, JOBVR, N, ILO, IHI, ..., WORK(IWRK),
$            LWORK+1-IWRK, 0, IERR )
```

`dggev3.f`, `cggev3.f`, and `zggev3.f` all gate the query block on
`IF(ILV)` (or split it into separate `IF(ILVL)` and `IF(ILV)` blocks
the way `cggev3` does), matching the real call. Only `sggev3.f` uses
`ILVL` where it should use `ILV`.

**Affected files.**
- `external/lapack-3.12.1/SRC/sggev3.f` (workspace-query block ~line
  339).

**Fix.** Restructure the query block to mirror `dggev3.f`:

```fortran
IF( ILVL ) THEN
   CALL SORGQR( ..., -1, IERR )
   LWKOPT = MAX( LWKOPT, 3*N+INT( WORK( 1 ) ) )
END IF
IF( ILV ) THEN
   CALL SLAQZ0( 'S', ..., -1, 0, IERR )
ELSE
   CALL SLAQZ0( 'E', ..., -1, 0, IERR )
END IF
LWKOPT = MAX( LWKOPT, 2*N+INT( WORK( 1 ) ) )
```

**NOT patched** as a source override because S-half is non-canonical
for the migrator. The migrated extended-precision archive uses the
D-derived `QGGEV3`, which is correct. Standalone single-precision
LAPACK builds inherit the bug.

**Why upstream's tests miss it.** The `dchkgg`/`schkgg` test drivers
allocate workspace generously above the queried minimum, so the query
under-estimate doesn't manifest as a fault — the real call still has
enough room. The bug is invisible unless a caller follows the
documented minimum-workspace contract literally.

**Upstream report.** Not yet filed.

---

## LAPACK 3.12.1: `slarf1l.f` missing `LASTV.GT.0` guards around DSCAL paths

**Symptom.** Out-of-bounds memory access in `SLARF1L` when the
input matrix has degenerated to `LASTV=0` (no significant rows/columns
left after the leading-zero scan). `dlarf1l.f` wraps the inner
"apply Householder reflector" block in `IF(LASTV.GT.0) THEN`, so
`LASTV=0` falls through to the `RETURN`. `slarf1l.f` lacks that
outer guard: with `LASTV=0`, the code reaches
`IF(LASTV.EQ.FIRSTV) THEN; CALL SSCAL( LASTC, ONE-TAU, C(LASTV,1), LDC )`
where `C(LASTV,1) = C(0,1)` is an out-of-bounds Fortran array access.

**Root cause.** Compare `dlarf1l.f:192` and `slarf1l.f:193`:

```fortran
* dlarf1l.f (correct)
IF( LASTV.GT.0 ) THEN
   IF( LASTV.EQ.FIRSTV ) THEN
      CALL DSCAL(LASTC, ONE - TAU, C(FIRSTV, 1), LDC)
   ELSE
      ...
   END IF
END IF

* slarf1l.f (missing outer guard)
IF( LASTV.EQ.FIRSTV ) THEN
   CALL SSCAL( LASTC, ONE - TAU, C(LASTV, 1), LDC )
ELSE
   ...
END IF
```

Same omission appears at `slarf1l.f:223` (the `Form C * H` branch):
the inner block lacks the `IF(LASTV.GT.0)` wrapper that
`dlarf1l.f:222` has. Additionally, the corner-case SSCAL uses
`C(LASTV, 1)` where `dlarf1l` uses `C(FIRSTV, 1)`; in the branch
where `LASTV.EQ.FIRSTV` they're the same element so it doesn't
matter, but the slarf1l form looks the wrong way around once the
outer `LASTV.GT.0` guard is absent.

**Affected files.**
- `external/lapack-3.12.1/SRC/slarf1l.f` (lines 193, 197, 223, 227).

**Fix.** Add `IF(LASTV.GT.0) THEN ... END IF` wrappers around both
inner blocks; change `C(LASTV, 1)` and `C(1, LASTV)` to
`C(FIRSTV, 1)` and `C` (consistent with dlarf1l).

**NOT patched** as a source override because S-half is non-canonical
for the migrator. The migrated extended-precision archive uses the
D-derived `QLARF1L`, which has the guards. Standalone single-
precision LAPACK builds inherit the OOB read.

**Why upstream's tests miss it.** `LASTV=0` requires the input matrix
to be effectively zero in the trailing rows/columns scanned by
`ILASLR`/`ILASLC`. The reference LAPACK test drivers seed C with
random data, so the leading-zero scan rarely returns 0. The bug
surfaces when SLARF1L is fed degenerate input from a calling driver
whose iteration has reduced C to zeros.

**Upstream report.** Not yet filed.

---

## LAPACK 3.12.1: `slarrf.f` shift-safety expansion uses TWO×EPS where DLARRF uses FOUR×EPS

**Symptom.** Numerical correctness drift in single-precision
eigenvalue refinement (`?STEMR` / RRR algorithm). The "make sure
shift bounds are properly outside the cluster" expansion uses a
narrower margin in the S half than in the D half:

```fortran
* slarrf.f:277-278
LSIGMA = LSIGMA - ABS(LSIGMA) * TWO * EPS
RSIGMA = RSIGMA + ABS(RSIGMA) * TWO * EPS

* dlarrf.f:277-278
LSIGMA = LSIGMA - ABS(LSIGMA) * FOUR * EPS
RSIGMA = RSIGMA + ABS(RSIGMA) * FOUR * EPS
```

`TWO*EPS` is the natural lower bound for safe shift placement in
single precision; `FOUR*EPS` (the D-half value) reflects a later
robustness fix that doubled the margin. The S half didn't get the
fix. With the narrower margin, SLARRF can occasionally place a
shift right at a cluster boundary instead of safely outside it,
causing the inner SSTEMR loop to stall or report convergence on a
shifted-but-overlapping interval.

**Affected files.**
- `external/lapack-3.12.1/SRC/slarrf.f` (lines 277-278).

**Fix.** Change `TWO * EPS` to `FOUR * EPS` on both lines. **NOT
patched** — S-half non-canonical.

**Upstream report.** Not yet filed.

---

## LAPACK 3.12.1: `slasq2.f` hard-codes IEEE=.FALSE. instead of querying ILAENV

**Symptom.** Loss of an IEEE-arithmetic optimization path in
`SLASQ2` (the dqds eigenvalue iteration used by SBDSQR). The IEEE
fast path runs without the explicit zero-test that the
non-IEEE-safe path performs each iteration; on IEEE-compliant
hardware (which is essentially all modern systems) the fast path
should be selected.

`dlasq2.f:275`:
```fortran
IEEE = ( ILAENV( 10, 'DLASQ2', 'N', 1, 2, 3, 4 ).EQ.1 )
```

`slasq2.f:274-279`:
```fortran
*     IEEE = ( ILAENV( 10, 'SLASQ2', 'N', 1, 2, 3, 4 ).EQ.1 )    ! commented out
...
IEEE = .FALSE.    ! hard-coded
```

The S half has the proper ILAENV query commented out and forces
`IEEE = .FALSE.`. Performance regression only — the result is the
same, just slower (every divisor is checked for zero). Most likely
a leftover from an upstream debugging session; the comment placement
suggests it was meant to be temporary.

**Affected files.**
- `external/lapack-3.12.1/SRC/slasq2.f` (lines 274-279).

**Fix.** Uncomment the ILAENV call and remove the `IEEE = .FALSE.`
override. **NOT patched** — S-half non-canonical and the migrator
already replaces ILAENV with ILAENV_EP at all call sites in the
migrated archive, so the D-derived QLASQ2 picks the proper path
regardless.

**Upstream report.** Not yet filed.

---

## LAPACK 3.12.1: `sgelss.f` uses ILAENV-estimated workspace where DGELSS queries directly

**Symptom.** Workspace size estimate returned by `SGELSS` for the
`M < N` path uses `ILAENV(1, 'SGELQF', ...)` to look up the optimal
block size and then computes `MAXWRK = M + M*NB`. `DGELSS`
instead queries DGELQF directly with `LWORK=-1` and uses the
returned value:

```fortran
* sgelss.f:324
MAXWRK = M + M*ILAENV( 1, 'SGELQF', ' ', M, N, -1, -1 )

* dgelss.f:308-310, 328
CALL DGELQF( ..., DUM, -1, INFO )
LWORK_DGELQF = INT( DUM(1) )
...
MAXWRK = M + LWORK_DGELQF
```

The query form is more accurate when the optimal LWORK includes
overhead beyond `M*NB` (e.g. when DGELQF needs space for triangular
factor T or for a workspace lookahead block). For modest M,N the
two formulas coincide, but for large problems the ILAENV-estimated
form can under-state the required workspace by a small constant.

**Severity.** Borderline. The reference DGELQF/SGELQF
implementations don't currently use significantly more than `M*NB`
overhead, so the under-estimate doesn't typically bite. But it's
inconsistent with the D-half pattern and could break if a vendor
LAPACK substitutes a different SGELQF whose workspace overhead is
higher.

**Affected files.**
- `external/lapack-3.12.1/SRC/sgelss.f` (line 324).

**Fix.** Replace the `M*ILAENV(...)` formula with a proper SGELQF
query, mirroring `dgelss.f:308-310, 328`. **NOT patched** — S-half
non-canonical.

**Upstream report.** Not yet filed.

---

## LAPACK 3.12.1: `sgejsv.f` JOBA='L' should be 'G' in LSVEC && RSVEC branch

**Symptom.** Numerical-correctness drift in the float SVD path
when both left and right singular vectors are requested and the
`U` matrix has just been *populated* (not zeroed) in its strict
upper triangle. The ELSE branch immediately above the call uses
`SLASET('U', NR-1, NR-1, ZERO, ZERO, U(1,2), LDU)` to zero the
strict upper, but the corresponding IF branch *fills* it with
`U(p,q) = -SIGN(TEMP1, U(q,p))`. After either branch, the same
SGESVJ call processes U:

```fortran
CALL SGESVJ( 'L', 'U', 'V', NR, NR, U, LDU, SVA, ...)   ! sgejsv:1711
```

`'L'` tells SGESVJ that `U` is lower triangular and the strict
upper is zero. After the IF branch, that's not true — SGESVJ
silently drops the upper-triangular data the IF branch just
deposited. DGEJSV correctly uses `'G'` (general):

```fortran
CALL DGESVJ( 'G', 'U', 'V', NR, NR, U, LDU, SVA, ...)   ! dgejsv:1711
```

Mathematically: `'G'` is correct for both branches; `'L'` is only
correct for the ELSE (zeroed) branch. The IF branch produces a
biased SVD result.

**Root cause.** `sgejsv.f:1711`. Likely lost during a sync between
the D and S copies: D got the `'L' → 'G'` fix, S didn't.

**Affected files.**
- `external/lapack-3.12.1/SRC/sgejsv.f` (line 1711, 1st arg of
  SGESVJ).

**Fix.** Single-token `'L' → 'G'`. **NOT patched** as a source
override because S-half is non-canonical for the migrator. The
extended-precision archive uses the D-derived QGEJSV, which is
correct. Standalone single-precision LAPACK builds inherit the bug.

**Why upstream's tests miss it.** `dchksvj` shape coverage probably
doesn't trigger the IF branch frequently enough for the biased
result to fall outside the test's tolerance bound, or the tolerance
is loose enough to absorb the bias.

**Upstream report.** Not yet filed.

---

