# Upstream bugs in vendored Netlib sources

*Last catalogued: 2026-05-06*

This document catalogues bugs found in the vendored upstream sources
(`external/lapack-3.12.1/`, `external/scalapack-2.2.3/`, etc.) that
the migrator works around without editing `external/`. Each entry
records the symptom, root cause, and the in-tree workaround. Entries
that have been reported to upstream link the tracking issue.

## How fixes are carried

Recipes accept a ``source_overrides`` field (see ``recipes/README.md``)
that maps an upstream filename to a replacement source written in
upstream shape (``DOUBLE PRECISION`` types, ``pd*``/``pz*`` symbol
names, ``dgemm`` call sites, …). The replacement goes through the
normal migration pipeline, so a single override produces correctly
renamed/promoted output for every target. The standard-precision
archive built from the unmodified ``external/`` tree is unaffected
— only the migrated extended-precision archive carries the fix.

When the convergence picker would otherwise pick the un-fixed C/S
half over the patched D/Z half, the recipe's ``prefer_source`` field
pins the correct canonical (the rank picker doesn't recognize
ScaLAPACK's ``pd*``/``pz*`` two-letter prefix — first character is
always ``P`` — so it sorts alphabetically by file name).

---

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
`recipes/lapack/source_overrides/{d,s,z,c}orbdb3.f` and wired via
`recipes/lapack.yaml`'s `source_overrides:` map. Every migrated target
(kind10, kind16, multifloats) builds against the patched form. The
`prefer_source` pin isn't needed here — LAPACK's canonical-rank picker
prefers D over the other halves by default, and the override matches
that choice.

**Standard-precision archive still buggy.** The kind4 / kind8 baseline
path (`migrator stage --target kind{4,8}`) skips migration and stages
upstream `external/lapack-3.12.1/SRC/` verbatim into `_reflapack_src/`,
so the std `lapack` archive that the baseline links against still has
the typo. Two of the residual failures in
`tests/KIND48_BASELINE_STATUS.md` (`zunbdb3`, `dorbdb3`) trace to this
gap and are parked under `tests/lapack/TODO.md` until the staging path
overlays `recipes/lapack/source_overrides/` on top of the upstream
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
Carried in `recipes/lapack/source_overrides/dgejsv.f`. Wired via
`recipes/lapack.yaml`'s `source_overrides:` map. (No `prefer_source`
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
Carried in `recipes/lapack/source_overrides/dtrsyl3.f` and
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

**Carried in:** `recipes/lapack/source_overrides/{dlaqp2,zlarf1f,zhetrf_aa,zlahef_aa,zgbrfsx}.f`,
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

## ScaLAPACK 2.2.3: `p{d,z}atrmv_.c` ALPHA hardcoded to one in (UPLO=L, TRANS=T/C)

**Symptom.** `PDATRMV` / `PZATRMV` (the absolute-value triangular
matrix-vector product `y := |alpha|*|A|*|x| + |beta|*|y|`) returns
results disagreeing with a serial reference whenever `UPLO='L'` and
`TRANS='T'` (resp. `TRANS='T'` / `'C'` for the complex variant) are
combined with a non-unit `ALPHA`. Every other (UPLO, TRANS)
combination matches. The differential-precision suite catches the
disagreement at `n=80` (10 block-rows of `MB=8`) with `alpha=0.6`;
the residual scales as `(1 - alpha) * |L_below^T * x_below|`.

**Root cause.** `pdatrmv_.c:573` and `pzatrmv_.c:577` invoke the
local off-diagonal contribution via `?agemv_` with the literal `one`
in the ALPHA slot:

```c
dagemv_( TRANS, &Amp0, &Anq0, one,
         Aptr, &Ald, ... );
```

Every other (UPLO, TRANS) branch in the same routine forwards the
caller's `ALPHA` (cast to the BLAS char-pointer ABI). The diagonal
block correctly receives `ALPHA` via `PB_Cptrm`, so the diagonal is
scaled and the off-diagonal isn't — produces `y_got = y_ref + (1 -
alpha) * |L_below^T * x_below|`.

**Affected files.**

* `external/scalapack-2.2.3/PBLAS/SRC/pdatrmv_.c` (line 573 — D-half).
* `external/scalapack-2.2.3/PBLAS/SRC/pzatrmv_.c` (line 577 — Z-half;
  same bug fires for both `TRANS='T'` and `TRANS='C'`).
* `external/scalapack-2.2.3/PBLAS/SRC/p{s,c}atrmv_.c` carry the
  byte-identical bug. Not exercised by us (we don't migrate single
  precision), so no override is wired for them.

**Fix.** One-token change: pass `((char *) ALPHA)` instead of `one`
to `?agemv_`. Carried in
`recipes/pblas/source_overrides/p{d,z}atrmv_.c` and wired via
`recipes/pblas.yaml`'s `source_overrides:` block. The migrator's
PBLAS pipeline applies the patched form for every extended-precision
target.

**Standard-precision archive still buggy.** Same caveat as the LAPACK
`?orbdb3` entry above: the std `pblas` archive built directly from
`external/` carries the upstream form. Standard-precision callers
of `PDATRMV` / `PZATRMV` see the same numerical shortfall.

**Why upstream's tests miss it.** The PBLAS test driver's input
shapes happen to combine (UPLO='L', TRANS='T') only with `ALPHA=1`
(or sets that exercise the asymmetry doesn't surface in). Any caller
exercising the documented contract with `ALPHA != 1` reproduces the
bug immediately.

**Upstream report.** Not yet filed.

---

## ScaLAPACK 2.2.3: `p?lanhs.f` NPROW=1 underestimate (1/F/I norms)

**Symptom.** Migrated `pqlanhs` / `pxlanhs` (and the upstream halves
they came from) return 1-norm, Frobenius-norm, and infinity-norm
values 10–20% smaller than the reference for upper-Hessenberg matrices
of size n ≥ 32 with MB ≥ 8. The max-element norm (`'M'`) appears to
pass — but only by luck on random matrices. The error is independent
of process count: it reproduces with `mpirun -np 1` as well as 2×2 grids.

**Root cause.** The NPROW=1 first-block code fails to advance the
local row counter `II` after processing the first block of columns.
The inner-loop bound `MIN(II + LL - JJ + 1, IIA + NP - 1)` is
supposed to stop at the local row of column `LL`'s subdiagonal —
which depends on `II` tracking the local row corresponding to the
top of the current column block. The structure is:

```fortran
IF( NPROW.EQ.1 ) THEN
   IF( MYCOL.EQ.IACOL ) THEN
      DO LL = JJ, JJ+JB-1
         ...inner loop bounded by II+LL-JJ+1...
      END DO
      JJ = JJ + JB             ! JJ advances
   END IF
   IACOL = MOD( IACOL+1, NPCOL )
   ! II is *not* advanced here — bug.

   DO J = JN+1, JA+N-1, NB
      ...inner loop using stale II...
      JJ = JJ + JB
      II = II + JB             ! main loop advances II every iteration
   END DO
END IF
```

When control enters the main loop for the second column block:

* `JJ` has been advanced to `JN + 1` (correct)
* `II` is still `IIA` (should be `IIA + JB`)

So for the first column of block 2 (`LL = JN + 1`):

```
II + LL - JJ + 1 = IIA + (JN+1) - (JN+1) + 1 = IIA + 1
```

That's row 2. The correct subdiagonal for column `JN + 1` is at
row `JN + 2` (e.g., row 10 for `JN = 8`). The inner loop reads only
2 rows where it should read 10 — eight elements per column are
silently dropped. After the main loop's first iteration `II` finally
advances, but it's still `JB` short of where it should be, and that
gap propagates for the rest of the matrix.

The NPROW>1 first-block code already does the analogous
`IF MYROW.EQ.IAROW THEN II = II + JB` advance after its own
first-block code, so that path is correct.

**Why M-norm passes by luck.** `MAX(|A(i,j)|)` doesn't care if you
skip elements as long as the actual maximum sits inside the rows you
do read. For random uniform entries the max element typically lands
in the upper-left, which is always in the kept range. M-norm tests
silently agreeing with the reference is not evidence that the code is
correct — it's evidence that the test inputs aren't adversarial.

The 1/F/I-norms are sums (or sums of squares); dropped elements
directly reduce the result by the missing fraction.

**Affected files.** All four precision halves carry the identical
buggy NPROW=1 path:

* `external/scalapack-2.2.3/SRC/pdlanhs.f` (used by our migrated D-half)
* `external/scalapack-2.2.3/SRC/pzlanhs.f` (used by our migrated Z-half)
* `external/scalapack-2.2.3/SRC/pslanhs.f` (S-half — bug present but
  not exercised by our extended-precision targets, since migration
  picks D as canonical for the real family)
* `external/scalapack-2.2.3/SRC/pclanhs.f` (C-half — same status)

The single-precision halves remain buggy in our standard-precision
archive (`libblas`, `liblapack`, …) since we link those directly
from `external/`. Standard-precision callers see the upstream
behavior. Only the migrated extended-precision archives carry the
fix.

**Fix.** Add the missing `II = II + JB` after the first-block code
in each of the four norm branches (M / 1 / I / F). Mirrors the
NPROW>1 branch's existing update.

```fortran
IF( NPROW.EQ.1 ) THEN
   IF( MYCOL.EQ.IACOL ) THEN
      DO LL = JJ, JJ+JB-1
         ...
      END DO
      JJ = JJ + JB
   END IF
   II = II + JB                 ! ← add this line
   IACOL = MOD( IACOL+1, NPCOL )
   ...
END IF
```

After the fix, all four norms agree with the reference to full
target precision: ~33 digits on KIND=16, ~19 on KIND=10, ~32 on
multifloats double-double.

**Workaround in tree.**

* `recipes/scalapack/source_overrides/pdlanhs.f`
* `recipes/scalapack/source_overrides/pzlanhs.f`

Wired via `recipes/scalapack.yaml`'s `source_overrides:` map.
`PDLANHS` and `PZLANHS` are pinned in `prefer_source:` so the
patched D/Z halves win convergence over the un-fixed C/S siblings.

**Why upstream Netlib's test suite never caught it.** Their tests
generate random matrices and lean heavily on M-norm coverage; the
1/F/I-norm paths on truly upper-Hessenberg input simply aren't
exercised. The disagreement only surfaces when you generate a
genuinely Hessenberg matrix (zero below the subdiagonal) and compare
the sum-norms to a serial reference — exactly what
`tests/scalapack/auxiliary/test_p[dz]lanhs.f90` do.

**Test drivers.**

* `tests/scalapack/auxiliary/test_pdlanhs.f90` — real Hessenberg, all four norms.
* `tests/scalapack/auxiliary/test_pzlanhs.f90` — complex Hermitian-Hessenberg, all four norms.

Both PASS to full target precision on all three targets after the fix.

**Upstream report.** Not yet filed.

---

## ScaLAPACK 2.2.3: `p?lanhs.f` IAROW double-advance (NPROW>1)

**Symptom.** A sibling bug to the NPROW=1 one above, hiding in the
NPROW>1 branch of the same routines. With a 2×2 (or any NPROW>1)
grid, the 1/F/M norms of an upper-Hessenberg matrix come out 20–50%
under-reported; I-norm sometimes matches by chance because every row
sum is contributed exactly once (just by the wrong rank). With
NPROW=1 the bug is silent because `MOD(_, 1) == 0` collapses the
miscompute.

**Root cause.** The first-block + per-block-iteration code maintains
a tracked owner row `IAROW` that should advance one row per block.
The upstream pattern is:

```fortran
INXTROW = MOD( IAROW + 1, NPROW )    ! before the first-block code
...first-block work...
IAROW = INXTROW                      ! step IAROW forward by one block
IAROW = MOD( IAROW + 1, NPROW )      ! WRONG: advances IAROW *again*
```

`INXTROW` is computed once before the first-block code and never
recomputed inside the main loop, so the second assignment is a typo
for `INXTROW = MOD( INXTROW + 1, NPROW )`. As written, `IAROW` skips
one row owner per iteration; for NPROW=2 this leaves `IAROW` stuck
on row 0 forever, so half the block rows are read by the "this row
owns the diagonal element" branch on the wrong rank. The 1/F norms
sum the dropped contributions to zero and underestimate; the M-norm
similarly misses entries.

The pattern repeats eight times in each of `pdlanhs.f` / `pzlanhs.f`:
once per norm (M, 1, I, F) × {first-block code, main loop body}.

**Affected files.** Same set as the NPROW=1 bug above — all four
precision halves (`p[sdcz]lanhs.f`) carry the identical pattern.

**Fix.** Replace the duplicate `IAROW = MOD(IAROW+1, NPROW)` with
`INXTROW = MOD(INXTROW+1, NPROW)` so `INXTROW` advances each
iteration and IAROW correctly tracks the next-block owner:

```fortran
IAROW   = INXTROW
INXTROW = MOD( INXTROW + 1, NPROW )
```

**Workaround in tree.** Same pair as the NPROW=1 fix:

* `recipes/scalapack/source_overrides/pdlanhs.f`
* `recipes/scalapack/source_overrides/pzlanhs.f`

Both the NPROW=1 `II = II + JB` patch and this IAROW fix live in the
same override file.

**Why upstream's tests miss it.** The upstream LIN driver feeds
random general matrices and matches against `DLANGE`, not against a
hand-computed reference for genuinely upper-Hessenberg input. Since
`PDLANHS` returns a value that's a slight underestimate on random
matrices, the driver's "result is finite and not too far from
DLANGE's value" assertion accepts it. Only when the matrix is
zeroed below the first subdiagonal (so the inner-loop bound
`II + LL - JJ + 1` actually constrains anything) does the bug
surface, which is what our `test_p[dz]lanhs.f90` drivers do.

**Upstream report.** Not yet filed.

---

## ScaLAPACK 2.2.3: `p?geequ.f` column-scale reduction wrong axis

**Symptom.** `PDGEEQU`/`PZGEEQU` return correct row scaling and
`AMAX` on multi-MYCOL grids, but `COLCND` and the `INFO` zero-column
detection are computed per grid column instead of globally. For a
2×2 grid, ranks `(_,0)` see a `COLCND` derived from columns 0..NB-1,
2*NB..3*NB-1, ...; ranks `(_,1)` see `COLCND` from the other set.
With our `tests/scalapack/auxiliary/test_p[dz]geequ.f90` random
matrices the disagreement is typically ~5% on `COLCND` and on
individual `C(j)` values that depend on the global RCMAX. The
real-precision test passes by chance on the seeds we use; the
complex test fails consistently.

**Root cause.** The C scale-factor computation has two reductions:

```fortran
! Each rank computes max over its local column-block:
DO J = JJA, JJA+NQ-1
   DO I = IIA, IIA+MP-1
      C(J) = MAX( C(J), CABS1( A(IOFFA+I) ) * R(I) )
   END DO
END DO
! Combine per-column maxes across grid rows (same MYCOL):
CALL DGAMX2D( ICTXT, 'Columnwise', COLCTOP, 1, NQ, C(JJA), 1, ...)

! Now compute scalar RCMAX/RCMIN over the local NQ columns:
DO J = JJA, JJA+NQ-1
   RCMAX = MAX( RCMAX, C(J) )
   RCMIN = MIN( RCMIN, C(J) )
END DO
! Combine across processes — this is the bug:
CALL DGAMX2D( ICTXT, 'Columnwise', COLCTOP, 1, 1, RCMAX, 1, ...)
CALL DGAMN2D( ICTXT, 'Columnwise', COLCTOP, 1, 1, RCMIN, 1, ...)
```

After the first `'Columnwise'` reduction on `C(JJA..)`, every rank in
a given grid column already holds the max-folded `C` values for its
column block. The second reduction (on the scalar `RCMAX`) is also
`'Columnwise'` — combining across `MYROW`-with-same-`MYCOL`, ranks
that already agree. So it's a no-op, and ranks holding *different*
column blocks (different `MYCOL`) keep their disjoint local
extrema, leading each grid column to compute `COLCND` from its own
block max/min instead of the global. The `IGAMX2D` for `INFO` (zero
column detection) has the same axis mistake.

The mirror axis would be `'Rowwise'` with `ROWCTOP`, which combines
across `MYCOL`-with-same-`MYROW` — exactly what's needed to merge
the disjoint per-grid-column blocks into a single global value.

The R scale-factor reductions are correct: after a `'Rowwise'` C(R)
reduce across `MYCOL`, the scalar RCMAX reduce uses `'Columnwise'` to
combine across `MYROW`, which is the right mirror.

**Affected files.**

* `external/scalapack-2.2.3/SRC/pdgeequ.f`
* `external/scalapack-2.2.3/SRC/pzgeequ.f`
* `external/scalapack-2.2.3/SRC/psgeequ.f` (same bug, untested by us)
* `external/scalapack-2.2.3/SRC/pcgeequ.f` (same bug, untested by us)

Lines (in `pdgeequ.f`): 332, 334, 346 — the three calls to
`DGAMX2D / DGAMN2D / IGAMX2D` after the C(JJA) accumulation block.
The same three lines exist in `pzgeequ.f` (with `DGAMX2D` /
`DGAMN2D` / `IGAMX2D` — note the `D`, not `Z`, for the real-typed
extrema in the complex routine).

**Fix.** Change `'Columnwise'` to `'Rowwise'` and `COLCTOP` to
`ROWCTOP` on those three calls. The patched override in
`recipes/scalapack/source_overrides/p[dz]geequ.f` carries the fix
plus an inline comment block explaining the mirror.

**Why upstream's tests miss it.** The driver compares row/column
scalings against an LU-based condition-number proxy, not against
sequential `DGEEQU` directly. The block-local `COLCND` happens to
fall in the same "no scaling needed" bucket for random matrices,
so the no-op path through the rest of the driver agrees. A
seed-by-seed numerical comparison surfaces the disagreement
immediately, which is what `test_p[dz]geequ.f90` do.

**Upstream report.** Not yet filed.

## ScaLAPACK 2.2.3: `p?posvx.f` LWMIN too small (PDPOCON / PZPOCON aborts)

**Symptom.** Calling migrated `pqposvx` / `pxposvx` (or the upstream
`PDPOSVX` / `PZPOSVX`) in single-thread or any nontrivial grid
configuration aborts inside the internal `PDPOCON` / `PZPOCON` call
with `On entry to P[QX]POCON parameter number 10 had an illegal value`.
The outer `*POSVX` workspace query returns `WORK(1) = 3*DESCA(LLD_)`,
which the caller dutifully allocates — but `PDPOCON` then enforces a
much larger `LWMIN` and aborts.

**Root cause.** `pdposvx.f:430` and `pzposvx.f:429` set
``LWMIN = 3*DESCA( LLD_ )`` where the documented contract (line 311 of
either file's prologue) is

```
LWORK = MAX( PDPOCON( LWORK ), PDPORFS( LWORK ) )
```

`PDPOCON`'s actual `LWMIN` is

```
2*NPMOD + 2*NQMOD +
MAX( 2, MAX( NB*ICEIL(NPROW-1, NPCOL),
             NQMOD + NB*ICEIL(NPCOL-1, NPROW) ) )
```

and `PDPORFS`'s is `3*NPMOD`. For any nontrivial grid (e.g. 2×2 with
`NB ≥ 8`), the documented formula exceeds `3*DESCA(LLD_)`. Same shape
for `PZPOSVX` with `LWMIN_PZPOCON = 2*NPMOD + MAX(2, …)` (no `2*NQMOD`
term — that part lives in `LRWMIN`) and `LWMIN_PZPORFS = 2*NPMOD`. The
`LRWMIN = MAX( 2*NQ, NP )` formula in `pzposvx.f:430` is also too small
versus `LRWMIN_PZPOCON = 2*NQMOD` and `LRWMIN_PZPORFS = NPMOD`.

**Files affected.**
- `external/scalapack-2.2.3/SRC/pdposvx.f` (line 430).
- `external/scalapack-2.2.3/SRC/pzposvx.f` (lines 429–430).

**Fix.** Patched overrides in
`recipes/scalapack/source_overrides/p[dz]posvx.f` recompute
`NPMOD`/`NQMOD` (PDPOCON's unadjusted NUMROCs), then set
`LWMIN = MAX( PDPOCON_LWMIN, PDPORFS_LWMIN )` (and `LRWMIN` accordingly
for the complex variant). `recipes/scalapack.yaml` declares the
overrides plus the matching `prefer_source: PDPOSVX, PZPOSVX` pin to
keep the canonical-rank picker from selecting the un-fixed S/C halves.

**Companion bug — `pzposvx.f` LRWMIN too small.** The same file's
`LRWMIN = MAX( 2*NQ, NP )` doesn't cover `PZLANHE('1', ...)`'s
documented `RWORK` requirement of `2*Nq0 + Np0 + LDW`, which `PZPOSVX`
calls at line 610 to compute ANORM. The shortfall doesn't abort
immediately — `PZLANHE` writes past the end of the caller's RWORK and
corrupts the heap, surfacing as a `free()` SIGSEGV during cleanup or a
few iterations later. Fix: `LRWMIN >= 2*NQMOD + NPMOD + LDW` (LDW=0
on square grids; we use `NPMOD` as a safe upper bound).

**Why upstream's tests miss it.** The shipped `*posvx` test drivers
allocate `WORK` / `RWORK` from precomputed bounds that exceed the
flawed `LWMIN` / `LRWMIN`, so the in-tree workspace queries are never
the binding constraint. A driver that trusts the queries (the natural
call pattern) reproduces the LWORK abort immediately and the LRWORK
heap corruption a few iterations later.

**Upstream report.** Not yet filed.

---

## ScaLAPACK 2.2.3: `pdsyevx.f` / `pzheevx.f` LQUERY-path WORK/RWORK(1:3) early-write

**Symptom.** Calling `PDSYEVX` (or `PZHEEVX`) with `LWORK = -1`
(workspace query) and a 1-element `WORK` (the documented contract for
the query) corrupts the heap, surfacing as a `free(): invalid
pointer` or "corrupted size vs. prev_size" abort during cleanup of
any unrelated buffer. The query *also* returns the correct optimal
`LWORK` (and `LRWORK` for the complex variant) in `WORK(1)` /
`RWORK(1)`, so a follow-up properly sized call runs to completion and
produces correct eigenvalues — the heap damage shows up only at
deallocate time.

`PZHEEVX` carries the same shape on `RWORK(1:3)` (rank 0 broadcasts
`ABSTOL` / `VL` / `VU` into the real-typed `RWORK` before reporting
`LRWORK`); a 1-element `RWORK` is similarly trampled.

**Root cause.** Inside `PDSYEVX`'s LQUERY branch, before computing
the optimal `LWORK` to return in `WORK(1)`, rank 0 broadcasts the
real-typed `ABSTOL`, `VL`, and `VU` into `WORK(1:3)` so the other
ranks see the same values. The broadcast happens unconditionally
(both LQUERY and full-call paths execute it), so on the LQUERY path
`WORK(2)` and `WORK(3)` get written before the LQUERY return at the
end of the routine. A caller passing `WORK(1)` for the query —
exactly the documented two-pass contract (`WORK(1)` for the query,
then size to the returned `LWORK`) — has `WORK(2..3)` written past
the end of its 1-element buffer.

`PZHEEVX` does the same broadcast, but into `RWORK(1:3)` since
`ABSTOL` / `VL` / `VU` are real-typed (the complex variant's `WORK`
is `complex`).

**Files affected.**
- `external/scalapack-2.2.3/SRC/pdsyevx.f` (rank-0 broadcast block
  ahead of the LQUERY return).
- `external/scalapack-2.2.3/SRC/pzheevx.f` — same shape on `RWORK`.
- `external/scalapack-2.2.3/SRC/{ps,pc}{sy,he}evx.f` carry the
  byte-identical pattern; not exercised by us (we don't migrate
  single precision).

**Fix.** Wrapper-side, in `tests/scalapack/common/target_scalapack_body.fypp`:
- `target_pdsyevx` / `target_pssyevx` / `target_p{q,e,m}syevx` allocate
  a local `work_t(3)` for the LQUERY branch and forward it to upstream
  instead of the caller's `WORK`; `work_t(1)` (= `LWMIN`) is copied
  back to `work(1)` after the query returns.
- `target_pzheevx` / `target_pcheevx` / `target_p{x,y,w}heevx`
  additionally allocate `rwork_t(3)` for the LQUERY branch.

Same pattern as the `target_pdtrsen` `iwork_t(max(1,n))` fix
documented in the next section. Both fixes live wrapper-side rather
than as `recipes/scalapack/source_overrides/` entries because they
adjust caller workspace expectations rather than the algorithm; they
don't need to ship into migrated builds the way an algorithm fix does.

**Why upstream's tests miss it.** ScaLAPACK's own test driver allocates
`WORK` / `RWORK` from a precomputed `LWMIN` upper bound *before*
the query, so the early-write fits comfortably. A caller that follows
the documented two-pass contract trips the heap corruption.

**Caveat — alternative interpretation.** `pdsyevx.f:247` documents
`WORK` as `dimension max(3, LWORK)`, and `pzheevx.f:282` says the
same about `RWORK`. Read strictly, the contract requires *every*
caller (LQUERY included) to provide at least 3 elements — the
ABSTOL/VL/VU broadcast then sits inside the documented buffer and
isn't an OOB write at all. Under that reading the wrapper-side
`work_t(3)` / `rwork_t(3)` allocation is plain contract compliance,
not a fix. The reason this entry stays in the catalogue: standard
LAPACK / ScaLAPACK LQUERY convention everywhere else is "WORK(1) is
the only element touched on the query," so a caller naturally
assuming that convention (as our wrappers initially did) gets bitten.
Whether the routines should be tightened to honour the universal
LQUERY convention or the doc clarified to flag the deviation is a
project-policy call; either way the surprise is real and the
wrapper-side allocation is the load-bearing fix.

**Upstream report.** Not yet filed.

---

## ScaLAPACK 2.2.3: `pdtrsen.f` IWORK(1:N) early-write during LQUERY

**Symptom.** Calling `PDTRSEN` with `LIWORK = -1` (workspace query)
and a 1-element `IWORK` corrupts the heap, surfacing as a
`free(): invalid pointer` abort on the next `deallocate` of any
unrelated buffer. The query *also* returns the correct `LIWMIN` in
`IWORK(1)`, so a follow-up full call with a properly sized `IWORK`
runs to completion and produces a numerically correct eigenvalue
spectrum — the heap damage shows up only at cleanup.

**Root cause.** The `LQUERY` gate at `pdtrsen.f:475` enters the
parameter-validation block on either `INFO=0 OR LQUERY`. Inside that
block, the `SELECT(1:N) → IWORK(1:N)` integer-conversion loop at
lines 499-525 writes `IWORK(K)` (and `IWORK(K+1)` for 2-by-2 block
boundaries) for `K = 1..N`. Line 537 then broadcasts `IWORK(1:N)` via
`IGAMX2D(...,IWORK,N,...)`. Only after all of this is the `LQUERY`
return reached at line 619-622, which sets `IWORK(1) = LIWMIN`. The
caller passing `IWORK(1)` for the query (the documented contract) has
`IWORK(2..N)` written past the end of its 1-element buffer.

**Files affected.**
- `external/scalapack-2.2.3/SRC/pdtrsen.f` (lines 499-538).

**Fix.** Wrapper-side, mirroring the `target_pdsyevx` `WORK_T(3)` bump
for upstream `PDSYEVX`'s `WORK(1:3)` early-write. `target_pdtrsen` in
`tests/scalapack/common/target_scalapack_body.fypp` now allocates a
local `iwork_t(max(1,n))` for the LQUERY branch and forwards it to
upstream instead of the caller's `IWORK`; `iwork_t(1)` (= LIWMIN) is
copied back to `iwork(1)` after the query returns.

**Why upstream's tests miss it.** Most callers (including ScaLAPACK's
own test harness) allocate `IWORK` to a precomputed `LIWMIN` upper
bound *before* the query, so the early-write fits. A caller that
follows the documented two-pass query contract (`IWORK(1)` for the
query, then size to `LIWMIN`) trips the heap corruption.

**Upstream report.** Not yet filed.

---

## ScaLAPACK 2.2.3: `pdlarzb.f` PBDTRAN N-arg vs LV-buffer mismatch

**Symptom.** `PDORMRZ` / `PZUNMRZ` (apply Q from RZ-factorization) on a
multi-rank grid silently produces a numerically corrupt C: gathered
result disagrees with LAPACK `dormrz`/`zunmrz` by a factor of ~1.2-1.3
regardless of SIDE/TRANS. Stderr shows
"On entry to PBDTRAN parameter number 11 had an illegal value" but
PXERBLA returns rather than aborting, so the routine continues with
uninitialised workspace and returns garbage. Single-rank execution
also trips the issue.

**Root cause.** `pdlarzb.f:374-377` allocates the PBDTRAN destination
buffer:

```fortran
LV = MAX( 1, MPC20 )    ! MPC20 = local rows of L-row distribution
WORK( IPV ) is MPC20 x K
```

then calls

```fortran
CALL PBDTRAN( ICTXT, 'Rowwise', 'Transpose', K, M+ICOFFV,
$             DESCV( NB_ ), WORK( IPW ), LW, ZERO,
$             WORK( IPV ), LV, IVROW, IVCOL, ICROW2, -1, WORK( IPT ) )
```

But `pbdtran.f:457-459` (the ROWFORM/'Rowwise' branch) checks
`LDC ≥ Np` where `Np = NUMROC(N, NB, MYROW, ICROW, NPROW)` is the
local-row count of an `N`-row distribution, and `N = M+ICOFFV` is
PDORMRZ's outer-call M (= the order of Q = `K + L` for the canonical
shape with `mA <= nA`). Since `mC = K + L > L` for any non-trivial K,
`LV < Np` and PBDTRAN refuses the call.

**Files affected.**
- `external/scalapack-2.2.3/SRC/pdlarzb.f` (lines 374-395).
- `external/scalapack-2.2.3/PBLAS/SRC/PBBLAS/pbdtran.f` (lines 457-459).
- Same shape for the complex variant `pzlarzb.f` / `pbztran.f` driving
  `pzunmrz`.

**Fix (PARTIAL).** Patched overrides in
`recipes/scalapack/source_overrides/p[dz]larzb.f`. The repair is a
one-line change to the SIDE='L' branch: PBDTRAN/PBZTRAN's N argument
goes from `M+ICOFFV` to `L+ICOFFV` (and the matching `MQV0` line is
adjusted likewise). Per `pdlarzb`'s comment "WORK(IPW) is K x MQV0 = [
. V(IOFFV) ]'" and the row-stored RZ convention, only the trailing
K×L portion of V holds meaningful Householder data, so the transpose
only needs to operate on those L columns. The `MPC20`-sized IPV
buffer (= local rows of L) then satisfies PBDTRAN's
`LDC >= NUMROC(L+ICOFFV, NBV, MYROW, ICROW2, NPROW)` requirement
under the alignment hypothesis `IROFFC2 = ICOFFV` (= `NB_V = MB_C`).
`recipes/scalapack.yaml` declares the overrides plus matching
`prefer_source: PDLARZB, PZLARZB` to keep the canonical-rank picker
from choosing the un-fixed S/C halves.

This fix alone is insufficient. A second, independent upstream bug in
`pdormrz.f` / `pzunmrz.f` (post-loop condition copy-paste — see the
section below) also applies; both fixes are required for SIDE='R' to
pass. SIDE='L' continues to fail by ~factor of 2 even after both
fixes — the residual bug appears to live in the PDORMR3/PDLARZ chain
for SIDE='L' and is still under investigation.

**Why upstream's tests miss it.** ScaLAPACK ships no `pdormrz`
differential test in its standard test suite. The routine is exercised
only as a building block of higher-level GSVD/GELSY drivers, where the
specific shape that triggers PBDTRAN's check apparently isn't hit.

**Test drivers.**
- `tests/scalapack/factorization/test_pdormrz.f90` — real, SIDE='R'
  only (TRANS='N','T'), mA=32, mC=48, nC=nA=64.
- `tests/scalapack/factorization/test_pzunmrz.f90` — complex
  counterpart with TRANS in {'N','C'}.

Both pass to ~target precision on kind16 / 2×2 grid after the fixes.
SIDE='L' is currently disabled in the test drivers pending the
remaining PDORMR3/PDLARZ investigation (see `tests/scalapack/TODO.md`).

**Upstream report.** Not yet filed.

---

## ScaLAPACK 2.2.3: `p?ormrz.f` post-loop condition copy-paste

**Symptom.** Calling `PDORMRZ` (or `PZUNMRZ`) with SIDE/TRANS
combinations `(LEFT,NOTRAN)` or `(RIGHT,TRANS)` returns C with the
leading partial block of reflectors (rows/cols `IA:I2-1`) unapplied.
The remaining K reflectors are applied correctly, so the disagreement
with the LAPACK serial reference scales with how many reflectors fell
into the leading partial block — for the canonical-shape case
`K = mA`, that's the entire first `MIN(K, MB)` reflectors.

**Root cause.** `pdormrz.f:458-459` reads

```fortran
IF( ( LEFT .AND. .NOT.NOTRAN ) .OR.
$    ( .NOT.LEFT .AND. NOTRAN ) ) THEN
   IB = I2 - IA
   ...
   CALL PDORMR3( SIDE, TRANS, MI, NI, IB, L, A, IA, JA, ... )
END IF
```

This is a *byte-for-byte copy* of the pre-loop condition at
`pdormrz.f:416-417`. Because the same condition gates both, both fire
together for `(LEFT,!NOTRAN)` and `(!LEFT,NOTRAN)`, and *neither*
fires for the complementary cases `(LEFT,NOTRAN)` and
`(!LEFT,!NOTRAN)`. The complementary cases need the post-loop
PDORMR3 call to handle the leading partial block left untouched by
the main DO loop's asymmetric I1/I2 bounds.

The intended structure mirrors `pdormrq.f:454`:

```fortran
IF( ( RIGHT .AND. TRAN ) .OR.
$    ( LEFT .AND. NOTRAN ) ) THEN
```

— exactly the complement of the pre-loop. Same shape applies in
`pzunmrz.f:459-460`.

**Affected files.**
- `external/scalapack-2.2.3/SRC/pdormrz.f` (lines 458-459).
- `external/scalapack-2.2.3/SRC/pzunmrz.f` (lines 459-460).
- Same shape in `psormrz.f` / `pcunmrz.f` (untested by us).

**Fix.** Patched overrides in
`recipes/scalapack/source_overrides/p[dz]ormrz.f` change the post-loop
condition to `(LEFT .AND. NOTRAN) .OR. (.NOT.LEFT .AND. .NOT.NOTRAN)`.
Wired via `recipes/scalapack.yaml` plus matching `prefer_source:
PDORMRZ, PZUNMRZ` pins.

**Upstream report.** Not yet filed.


## ScaLAPACK 2.2.3: `p?larz.f` MPV/NQV undersizing (heap overrun)

**Symptom.** On NPROW > 1 grids with non-uniform A row distribution
(e.g. NPROW = 3, 5, 6 with K/NB blocks not divisible by NPROW),
`PDORMRZ` / `PZUNMRZ` (and `PCUNMRZ` / `PSORMRZ`) crash with glibc
"corrupted size vs. prev_size" inside `free()` during cleanup, or
produce silently wrong results when the overrun lands somewhere
benign. Reproduces with the migrator's standard 3×2 / 1×3 grids.

**Root cause.** `pdlarz.f:288` (and the s/c/z byte-identical copies)
defines the receiver-buffer size for the SIDE='L' branch as

```fortran
MPV = NUMROC( L+IROFFV, DESCV( MB_ ), MYROW, IVROW, NPROW )
```

— derived from V's row distribution (`IVROW` / `IROFFV` / `MB_V`).
But the row-V SIDE='L' `PBDTRNV` calls a few lines later transpose V
into a column at sub(C2)'s procrow `ICROW2` and write the local row
count of sub(C2). When `PDORMR3` iterates `IV` across blocks of A,
the `PDLARZ` alignment restriction documented in its header is
violated: V's row distribution no longer matches sub(C2)'s. `MPV`
underestimates the receiver write by 1-3 elements per call, and the
unpacked write into Y overruns into the adjacent aux WORK region,
landing on glibc heap metadata.

The same shape appears in the SIDE='R' branch at `pdlarz.f:292` for
`NQV` (col-V branch, exposed by direct callers passing `INCV=1`):

```fortran
NQV = NUMROC( L+ICOFFV, DESCV( NB_ ), MYCOL, IVCOL, NPCOL )
```

**Affected files.**
- `external/scalapack-2.2.3/SRC/pdlarz.f`, `pclarz.f`, `pslarz.f`,
  `pzlarz.f` (lines ~288 and ~292).
- `external/scalapack-2.2.3/SRC/pzlarzc.f` (and `pclarzc.f`) — the
  Q**H variant called from `PZUNMR3` when `TRANS='C'` / `'H'`. Same
  byte-identical formulas.

**Fix.** Override redefines `MPV` / `NQV` in the relevant branches
from sub(C2)'s distribution (`L+IROFFC2` / `DESCC(MB_)` / `ICROW2` for
SIDE='L'; `L+ICOFFC2` / `DESCC(NB_)` / `ICCOL2` for SIDE='R'). The
upstream `LWMIN` already sizes the leading WORK region to `MPC0`
(full sub(C) local rows), which dominates the new `MPV`, so no
`LWMIN` change is needed. Wired via
`recipes/scalapack/source_overrides/p[dz]larz.f` and
`recipes/scalapack/source_overrides/pzlarzc.f`, with
`prefer_source: PDLARZ, PZLARZ, PZLARZC` pins.


## ScaLAPACK 2.2.3: `p?larz.f` ZAXPY stride uses LDA where vector wants 1

**Symptom.** Companion bug to the `MPV`/`NQV` undersizing above: even
on aligned-grid calls where the buffer is sized correctly, the
SIDE='L' update path inside `P?LARZ` skips elements in `WORK` and
can read/write past the buffer when `NQC2 > 1`.

**Root cause.** `pdlarz.f` makes ten `DAXPY` calls into `WORK` /
`WORK(IPW)` using `INCY = MAX(1, NQC2)`, e.g.

```fortran
CALL DAXPY( NQC2, ONE, C( IOFFC1 ), LDC,
$           WORK( IPW ), MAX( 1, NQC2 ) )
```

`WORK` here is a contiguous local vector, so the correct stride is
`1`. `MAX(1, NQC2)` is a legitimate leading dimension elsewhere
nearby (for matrix-shaped `DLASET` / `DGSUM2D`); the bug looks like a
copy-paste of that LDA into the AXPY's stride slot. Compare LAPACK's
`dlarz.f`:

```fortran
CALL DAXPY( N, -TAU, WORK, 1, C, LDC )
```

**Affected files.**
- `external/scalapack-2.2.3/SRC/p[dscz]larz.f` (six pairs in the
  d/s variants, five pairs in the c/z and `pzlarzc` variants).
- `external/scalapack-2.2.3/SRC/p[cz]larzc.f` (same shape).

**Fix.** Override changes `INCY` from `MAX(1, NQC2)` to `1` in every
AXPY into `WORK`. The `MAX(1, NQC2)` expression is retained where it
is the legitimate LD of the surrounding `LASET` / `GSUM2D` matrix
calls. Wired via the same `source_overrides` entries as the
`MPV`/`NQV` fix.


## ScaLAPACK 2.2.3: `pzlarz.f` / `pzlarzc.f` SIDE='L' missing ZLACGV; ZGERC should be ZGERU

**Symptom.** `PZUNMRZ` SIDE='L' returns numerically wrong results,
producing residuals on the order of 1.3–1.8× the input magnitude
(both `TRANS='N'` via `PZLARZ` and `TRANS='C'` via `PZLARZC`). The
real-precision analogue (`PDORMRZ` SIDE='L') is unaffected.
Reproduces bit-identically on 1, 2, and 4 ranks — pure local
arithmetic, not a distribution issue. Surfaced by
`tests/scalapack/factorization/test_pzunmrz.f90` after the
`MPV`/`NQV` and AXPY-stride fixes had been lifted.

**Root cause.** Each SIDE='L' branch of `PZLARZ` / `PZLARZC`
accumulates `w = v^H * sub(C)` via

```fortran
CALL ZGEMV( 'Conjugate transpose', MPV, NQC2, ONE,
$           C( IOFFC2 ), LDC, V, 1, ZERO, WORK( IPW ), 1 )
IF ( MYROW.EQ.ICROW1 )
$   CALL ZAXPY( NQC2, ONE, C( IOFFC1 ), LDC, WORK( IPW ), 1 )
...
CALL ZAXPY( NQC2, -TAULOC( 1 ), WORK( IPW ), 1, C( IOFFC1 ), LDC )
CALL ZGERC( MPV, NQC2, -TAULOC( 1 ), WORK, 1,
$           WORK( IPW ), 1, C( IOFFC2 ), LDC )
```

But `ZGEMV('C', C2, V)` computes `C2^H * V`, which equals
`conj(V^H * C2)` — not `V^H * C2`. The subsequent `ZAXPY` of `C1`
(unconjugated) into `WORK(IPW)` therefore mixes
`conj(V^H*C2) + C1` rather than producing `w = V^H*C2 + C1`. The
final `ZGERC` (which computes `A := alpha * x * conj(y^T) + A`) then
conjugates `WORK(IPW)` again, leaving the rank-1 update with the
imaginary part of `C2` sign-flipped relative to the math.

LAPACK's serial `zlarz.f` handles this correctly by calling
`ZLACGV(N, WORK, 1)` twice — once after copying the '1'-row of C
(to bring it to `conj(C1)`), then again after the GEMV (so
`WORK = w`, not `conj(w)`) — and using `ZGERU` (no conjugation) for
the rank-1 update:

```fortran
CALL ZCOPY( N, C, LDC, WORK, 1 )
CALL ZLACGV( N, WORK, 1 )
CALL ZGEMV( 'Conjugate transpose', L, N, ONE, C(M-L+1,1), LDC,
$           V, INCV, ONE, WORK, 1 )
CALL ZLACGV( N, WORK, 1 )
CALL ZAXPY( N, -TAU, WORK, 1, C, LDC )
CALL ZGERU( L, N, -TAU, V, INCV, WORK, 1, C(M-L+1,1), LDC )
```

ScaLAPACK's `PZLARZ` translation evidently inlined the GEMV trick
(swapping the `ZCOPY+ZLACGV` for `ZGEMV` with `beta=ZERO`) but
dropped the post-GEMV `ZLACGV` and used `ZGERC` instead of `ZGERU`.
The bug never bites the real path (`PDLARZ` uses `DGEMV`/`DGER`
which have no conj/no-conj distinction) and never bites SIDE='R'
(where the rank-1 update wants `ZGERC` legitimately, since
`H` from the right is `C := C - tau (Cv) v^H`).

**Affected files.**
- `external/scalapack-2.2.3/SRC/pzlarz.f` — 5 SIDE='L' MPV×NQC2
  sites.
- `external/scalapack-2.2.3/SRC/pzlarzc.f` — same 5 sites in the
  Q**H variant called from `PZUNMR3` when `TRANS='C'`.
- `external/scalapack-2.2.3/SRC/pclarz.f` and `pclarzc.f` carry the
  byte-identical buggy formulas in single-complex; the same fix
  would apply but is not wired here (we don't migrate single-complex
  RZ).

**Fix.** Override inserts `CALL ZLACGV(NQC2, accum, 1)` after each
SIDE='L' `ZGEMV('Conjugate transpose', ...)` (where `accum` is
`WORK(IPW)` or `WORK` depending on the branch) and replaces each
`ZGERC(MPV, NQC2, -TAULOC(1), x, 1, accum, 1, C(IOFFC2), LDC)` with
the corresponding `ZGERU(...)`. The SIDE='R' `ZGERC` calls
(`MPC2×NQV`) are left untouched. `ZLACGV` and `ZGERU` are added to
the `EXTERNAL` declaration in both files. Wired via the same
`source_overrides` entries and `prefer_source` pins as the
`MPV`/`NQV` and AXPY-stride fixes.

**Not in `../scalapack-bugfix/scalapack`.** The four prior `fix-*`
branches (`fix-larz-buffer-sizing`, `fix-larz-daxpy-stride`,
`fix-larzb-pbtran-buffer`, `fix-ormrz-post-loop-guard`) addressed
distribution and stride bugs that hit the real path first; this
algebraic bug is complex-only and was masked by the other failures
until they were lifted. Worth submitting upstream as a new
`fix-larz-side-l-conj` topic branch.


## ScaLAPACK 2.2.3: `p?larz.f` / `p?larzc.f` PBxTRNV reads M elements where only L are stored

**Symptom.** `P?LARZ` / `P?LARZC` SIDE='L' rowwise (and SIDE='R'
columnwise) calls `PBxTRNV` with the global vector length set to
`M` (resp. `N`), but only the trailing `L` entries of `V` are
physically stored — the leading `M-L` (resp. `N-L`) entries are
implicit zeros from the RZ Householder convention. `PBxTRNV` reads
`M` (or `N`) elements from `V(IOFFV)`, picking up garbage past the
stored tail. The transpose buffer ends up corrupt; the subsequent
`?GEMV`/`?GER` reads it. On grids with non-uniform row distribution
the OOB read can also touch heap-managed memory and produce
"corrupted size vs. prev_size" aborts.

**Why we didn't observe it under our overrides.** The companion
`MPV` / `NQV` redefinition (sub(C2)'s distribution rather than V's)
already constrains the subsequent `?GEMV` to the meaningful `L`
local rows, so the garbage in the trailing `M-L` slots of the
transpose buffer is never consumed. The specific NPROW=3 heap-abort
case that surfaced the bug upstream is not in our test matrix
(ctest runs at 4 ranks → 2×2 grid).

**Affected files.**
- `external/scalapack-2.2.3/SRC/pdlarz.f` (and `pslarz.f`).
- `external/scalapack-2.2.3/SRC/pzlarz.f` (and `pclarz.f`).
- `external/scalapack-2.2.3/SRC/pzlarzc.f` (and `pclarzc.f`).

**Fix.** Override changes the third argument of `PBxTRNV` from `M`
(rowwise SIDE='L') / `N` (columnwise SIDE='R') to `L` in all four
sites of each file. Same `N`-vs-`L` pattern as the `P?LARZB`
PBxTRAN fix landed earlier (commit `75714cb` upstream). Wired via
the same `source_overrides` entries.

---

## ScaLAPACK 2.2.3: `pzungql.f` / `pzunml2.f` post-loop `PB_TOPGET` should be `PB_TOPSET`

**Symptom.** `PZUNGQL` / `PZUNML2` (the Z-half QL-form orthogonal-
generator and apply routines) leak modified BLACS broadcast topology
state on return: a caller that observes `'Broadcast', 'Rowwise'` /
`'Broadcast', 'Columnwise'` after the call sees the algorithm's
internal `'I-ring'` / `' '` settings instead of the caller's original
topology. Subsequent BLACS operations from that context inherit the
wrong topology until the caller manually re-sets it. The C-half
(`PCUNGQL` / `PCUNML2`) does the correct restore.

**Root cause.** Save/restore copy-paste error. The standard
ScaLAPACK pattern is `PB_TOPGET` at entry to save, `PB_TOPSET` to
install algorithm-specific settings, then `PB_TOPSET` at exit to
restore the saved values. In `pzungql.f`:

```fortran
*  L247-248 — save original
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
*  L249-250 — install algorithm's needs
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise', 'I-ring' )
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', ' ' )
      ...
*  L292-293 — should be PB_TOPSET to restore
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
```

The exit-side calls re-read into the same locals (no-op) instead of
writing the saved values back. `pzunml2.f` carries the same pattern
on lines 394-395.

**Affected files.**
- `external/scalapack-2.2.3/SRC/pzungql.f` (lines 292-293).
- `external/scalapack-2.2.3/SRC/pzunml2.f` (lines 394-395).
- `external/scalapack-2.2.3/SRC/pcungql.f` — *not* affected (uses
  `PB_TOPSET` correctly).
- `external/scalapack-2.2.3/SRC/pcunml2.f` — *not* affected.

The S/D halves (`psorgql.f` / `psormql.f` / `pdorgql.f` /
`pdormql.f`) — checking is left to a follow-up; the immediate
migrator concern is the Z half because that's what the differential
suite exercises.

**Fix in-tree.** No `source_overrides` body — the C-half is already
correct. Recipe pins the C-half as canonical via
`prefer_source: PCUNGQL / PCUNML2` so the convergence picker takes
the C body and the migrator generates correct `Q`/`X`/`E`/`Y` clones
with `PB_TOPSET` at the restore site. Standard-precision archive
built from unmodified `external/` still has the bug.

**Why upstream's tests miss it.** The leak is silent on tests that
don't observe BLACS topology between calls — most ScaLAPACK test
drivers run a single algorithm per context and tear the grid down,
so the state-leak never surfaces as a test-detectable behavior.
Reproducing requires calling `PZUNGQL` (or `PZUNML2`), then
inspecting `PB_TOPGET` results, then calling another routine that
relies on default broadcast topology — a usage pattern the upstream
test cycle doesn't exercise.

**Upstream report.** Not yet filed.

---

## ScaLAPACK 2.2.3: `pzheevd.f` PCHK2MAT parameter index 11 should be 12

**Symptom.** Calling `PZHEEVD` with an illegal `DESCZ` (e.g. invalid
`MB_`/`NB_`, mismatched `CTXT_`) produces a `PXERBLA` diagnostic that
reports parameter number 11 — but `DESCZ` sits at position 12 in
PZHEEVD's signature (`JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ,
DESCZ, ...`); position 11 is `JZ`. The error path *does* fire, just
mislabels which argument the caller got wrong, surfacing as a
confusing diagnostic during the parameter-validation stage.
`PCHEEVD` is correct (passes 12). The bug doesn't affect numerics.

**Root cause.** Copy-paste during the Z-half adaptation of
`PCHEEVD`. `pzheevd.f:299-300`:

```fortran
CALL PCHK2MAT( N, 3, N, 3, IA, JA, DESCA, 7, N, 3, N, 3, IZ,
$              JZ, DESCZ, 11, 2, IDUM1, IDUM2, INFO )
```

The 16th argument of `PCHK2MAT` is `DESCBPOS0` — the parameter-number
to forward to `PXERBLA` if `DESCB` (= `DESCZ` here) is invalid.
`PCHEEVD` (line 299) passes `12`, which is the actual position of
`DESCZ`. `PZHEEVD` passes `11`.

**Affected files.**
- `external/scalapack-2.2.3/SRC/pzheevd.f` (line 300, 16th arg of
  `PCHK2MAT`).

**Fix.** Single-token change `11 → 12` on the 16th argument of the
`PCHK2MAT` call. Carried in
`recipes/scalapack/source_overrides/pzheevd.f`. Wired via
`recipes/scalapack.yaml`'s `source_overrides:` map plus a matching
`prefer_source: PZHEEVD` pin so the patched Z half wins convergence
over the un-fixed `PCHEEVD` sibling (which carries the *correct* 12,
but the migrator's canonical-rank picker doesn't know that and would
otherwise sort `pcheevd.f` alphabetically first).

**Why upstream's tests miss it.** The shipped test driver feeds
valid descriptors and never exercises the parameter-error path on
`DESCZ` specifically; the misreported parameter number is invisible
unless a caller deliberately passes an invalid `DESCZ` and inspects
PXERBLA output.

**Upstream report.** Not yet filed.

---

## ScaLAPACK 2.2.3: `pzunmbr.f` EXTERNAL declares PCHK1MAT but body calls PCHK2MAT

**Symptom.** None at runtime — Fortran `EXTERNAL` is advisory and the
linker resolves `PCHK2MAT` correctly via the global symbol table. The
typo surfaces only as a divergence against `pcunmbr.f` in the
convergence reports (the C half lists `PCHK2MAT` correctly).

**Root cause.** `pzunmbr.f:302`:

```fortran
EXTERNAL           BLACS_GRIDINFO, CHK1MAT, PCHK1MAT, PXERBLA,
$                   PZUNMLQ, PZUNMQR
```

The body calls `PCHK2MAT` four times (lines 511, 515, 521, 525) and
never calls `PCHK1MAT`. `PCHK1MAT` in the EXTERNAL list is dead — and
it's the wrong name: the C half's `pcunmbr.f:302` correctly lists
`PCHK2MAT`. Pure copy-paste error in upstream's Z-half edit.

**Affected files.**
- `external/scalapack-2.2.3/SRC/pzunmbr.f` (line 302).

**Fix.** Replace `PCHK1MAT` with `PCHK2MAT` in the EXTERNAL list.
Carried in `recipes/scalapack/source_overrides/pzunmbr.f`. Wired via
`recipes/scalapack.yaml`'s `source_overrides:` map plus a matching
`prefer_source: PZUNMBR` pin so the patched Z half wins convergence
over the (already-correct) C half (canonical-rank picker would
otherwise sort `pcunmbr.f` first).

**Why upstream's tests miss it.** The bug has no observable effect.
Only a textual diff between the two precision halves surfaces it.

**Upstream report.** Not yet filed.

---

## ScaLAPACK 2.2.3: `pssyevd.f` LQUERY misses LIWORK=-1

**Symptom.** Calling `PSSYEVD` with `LWORK >= 1` and `LIWORK = -1`
(querying just the integer workspace) is rejected as if the caller
supplied real values: the routine treats `LWORK >= 1` as a real call
and proceeds, but then fails to populate `IWORK(1)` with `LIWMIN`
because the LQUERY branch never fires. `PDSYEVD` (and the migrated
extended-precision `P{Q,E,M}SYEVD` halves derived from it) handle
the asymmetric query correctly — `LQUERY` triggers when *either*
`LWORK == -1` or `LIWORK == -1`.

**Root cause.** `pssyevd.f` defines

```fortran
LQUERY = ( LWORK.EQ.-1 )
```

`pdsyevd.f` (and the migrated D-derived halves) correctly defines

```fortran
LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
```

A caller asking only for `LIWMIN` via `LIWORK = -1` therefore gets
treated as a real call on the S half.

**Affected files.**
- `external/scalapack-2.2.3/SRC/pssyevd.f` (LQUERY definition).
- `external/scalapack-2.2.3/SRC/pcheevd.f` not affected (its LQUERY
  guard is correct — only the real S-precision half drifts).

**Fix.** None wired in-tree. The migrator picks `PDSYEVD` as canonical
for the real family (D-half over S-half by precision-rank policy), so
every migrated extended-precision build (`PQSYEVD`/`PESYEVD`/
`PMSYEVD`) gets the correct LQUERY guard automatically. The bug
remains in the standard-precision archive built directly from
`external/`. Documented here for completeness; no `source_override`
needed.

**Why upstream's tests miss it.** Their test driver always passes
both `LWORK = -1` and `LIWORK = -1` together, which the broken guard
catches via the `LWORK == -1` clause; the `LIWORK == -1` -only path
isn't exercised.

**Upstream report.** Not yet filed.

---

## ScaLAPACK 2.2.3: `pslaed3.f` clobbers user INFO and skips bounds guard

**Symptom.** Two distinct upstream-bug shapes in the S half of the
divide-and-conquer eigenvalue sub-step `PSLAED3`, both absent in
`PDLAED3`:

1. The S half writes `INFO = 0` inside an internal helper loop —
   clobbering the user-visible output `INFO` with a stale zero, so a
   non-zero error set earlier in the routine is silently overwritten
   before return. The D half correctly uses a local `IINFO` for the
   internal loop's status and leaves the public `INFO` untouched.

2. The S half's main loop indexes `D(I+J)` without first checking
   `I + J <= N`, so for large enough `I` the loop reads past the end
   of the eigenvalue array `D(1:N)`. The D half wraps the loop body
   in `IF (I+J <= N) THEN ... ENDIF`. The OOB read on the S path is
   typically a small step into adjacent stack/heap and doesn't
   crash; it either returns a slightly wrong updated eigenvalue or
   nothing observable depending on what sat past `D(N)`.

**Affected files.**
- `external/scalapack-2.2.3/SRC/pslaed3.f` — both bugs.
- `external/scalapack-2.2.3/SRC/pdlaed3.f` — clean (canonical for our
  migrated builds).

**Fix.** None wired. The migrator picks `PDLAED3` as canonical for
the real family, so every migrated extended-precision derivative
inherits the correct `IINFO=0` and `IF(I+J.LE.N)` guard. The S-half
bugs persist in the standard-precision archive built from
`external/` unchanged. Documented here for completeness.

**Why upstream's tests miss it.** The eigenvalue test drivers report
correctness via residuals computed against shipped reference
matrices, not against an alternate decomposition; the OOB read at
`D(I+J)` for `I+J > N` happens to fall in zero-initialized stack
padding on the test build hosts and produces no detectable numerical
shift. A driver running on a host where the trailing memory holds
non-zero data reproduces the disagreement.

**Upstream report.** Not yet filed.

---

## MUMPS 5.8.2: bad-input cases SIGSEGV instead of returning `INFOG(1) < 0`

**Symptom.** The MUMPS user guide documents that input-validation
failures populate `id%INFOG(1)` with a negative error code (e.g.
`-16` for an invalid `N`, `-6` for IRN/JCN out of range) and return
control to the caller. In practice, MUMPS 5.8.2 SIGSEGVs on every
invalid-input case the differential-precision suite tried, before
any diagnostic is written. Concretely:

- `id%N = -1` (or `id%N = 0`) at JOB=1 — SIGSEGV inside analysis,
  no `INFOG(1)` set.
- `id%NNZ = -1` (negative count) — SIGSEGV during pattern setup.
- `id%IRN(1) = N + 5` (out-of-range row index, high) — SIGSEGV
  inside the input-reformatting path.
- `id%IRN(1) = 0` (out-of-range row index, low — zero is invalid
  under 1-based indexing) — same crash signature.
- `id%JCN(1) = N + 3` / `id%JCN(1) = 0` — symmetric to IRN.
- `id%NNZ` mismatched against the actual `IRN/JCN/A` sizes
  (e.g. `NNZ = 0` with non-empty arrays, or `NNZ` larger than
  what the arrays hold) — SIGSEGV during reformatting.

The crashes are unconditional: the analysis / pattern-setup code
doesn't bounds-check the user's IRN/JCN/N/NNZ before indexing into
internal work arrays sized from those same values, so a wrong value
either dereferences past the end of a fresh allocation or hits an
arena boundary.

**Root cause (inferred — not from upstream).** `*MUMPS_ANA_DRIVER`
and the IRN/JCN reformatters trust the assembled-format inputs the
user sets in the `id` struct. There is a `MUMPS_ANA_F_CHECK` (and
similar) helper internally, but it isn't invoked early enough — by
the time it would run, the analysis code has already indexed past
end-of-buffer using N or IRN values it never validated. So the
manual's promised `INFOG(1) < 0` return is a documented-but-unwired
contract: the validation pass exists in name but doesn't gate the
indexing it's supposed to protect.

**Affected files.** Symptom is observed at the `?MUMPS` driver
entry point; the actual crashing site varies per case (analysis
driver, IRN/JCN reformatter, pattern-set work-array sizing). All
under `external/MUMPS_5.8.2/src/`. Not narrowed to a specific
file because the workaround sidesteps the issue at the caller
rather than patching MUMPS internals.

**Workaround (in-tree, not an override).** Per
`tests/mumps/TODO.md` D1, the differential-precision wrapper at
`tests/mumps/common/target_mumps_body.fypp` exports
`check_dmumps_input` / `check_zmumps_input` that mirror the
documented validation contract before any `?MUMPS` call. They
return integer codes (`MIC_BAD_N`, `MIC_BAD_NNZ`, `MIC_BAD_IRN`,
`MIC_BAD_JCN`, `MIC_SIZE_MISMATCH`, `MIC_OK`). Tests at
`tests/mumps/fortran/test_{d,z}mumps_errors.f90` exercise each
class plus a final valid-input pass that reaches `MIC_OK` and
factors via `JOB=6` to confirm the wrapper isn't over-rejecting.
This is a *test-side* contract simulator — it lets us assert the
documented behavior without crashing the test harness, but does
not modify MUMPS itself. A real caller still has to pre-validate
their inputs externally to avoid the SIGSEGV.

**Why upstream's tests miss it.** MUMPS's own test drivers
(`*examples/`) feed valid matrices only. The validation contract
documented in the user guide is asserted by the prose, not by any
shipped test that passes a deliberately bad `id` struct and
checks `INFOG(1)`. So the gap between "documented" and "wired" is
invisible to the upstream test cycle.

**Upstream report.** Not yet filed.
