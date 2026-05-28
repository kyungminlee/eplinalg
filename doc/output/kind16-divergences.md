# Notes — KIND=16 Divergence Analysis

This document provides a detailed technical analysis of convergence divergences
for the `kind16` target. For a current cross-target summary, run
`python -m migrator diverge <recipe>`.

---

## Remaining BLAS convergence divergences

After type-migration, three BLAS routine pairs disagree between the
single- and double-precision halves. These are not migrator bugs —
they reflect genuine source-level differences between the S/C and D/Z
files in the upstream LAPACK 3.12.1 reference BLAS. The migrator keeps
the D/Z version on disk (our canonical output) and reports the
mismatch.

### 1. `sdot.f` vs `ddot.f` → `qdot.f`

Trivial statement-order difference in the initializer.

```fortran
! sdot.f
STEMP = 0.0e0
SDOT  = 0.0e0

! ddot.f
DDOT  = 0.0d0
DTEMP = 0.0d0
```

The two statements are reordered between halves. Semantically
identical; only statement order differs.

### 2. `scasum.f` vs `dzasum.f` → `qxasum.f`

The complex-absolute-value computation is inlined in the single
version but delegated to the helper `DCABS1` in the double version.

```fortran
! scasum.f (inline)
STEMP = STEMP + ABS(REAL(CX(I))) + ABS(AIMAG(CX(I)))

! dzasum.f (helper call)
STEMP = STEMP + DCABS1(ZX(I))
```

`DCABS1(z) = |Re z| + |Im z|`, so the two formulations compute the
same quantity. Upstream has simply never been harmonized; the single
side was written before a matching `SCABS1` helper was considered.

### 3. `srotmg.f` vs `drotmg.f` → `qrotmg.f`

Two differences:

**(a) Declaration order.** The `DATA`-initialized constants are
grouped differently in the local variable list:

```fortran
! srotmg.f
REAL GAM,GAMSQ,ONE,RGAMSQ,SFLAG,SH11,SH12,SH21,SH22,SP1,
     SP2,SQ1,SQ2,STEMP,SU,TWO,ZERO

! drotmg.f
DOUBLE PRECISION DFLAG,DH11,DH12,DH21,DH22,DP1,DP2,DQ1,DQ2,
                 DTEMP,DU,GAM,GAMSQ,RGAMSQ,TWO,ZERO,ONE
```

Declaration ordering has no semantic effect.

**(b) Constant literal precision.** The two halves write the shared
overflow/underflow constants at different precisions:

| constant  | srotmg.f       | drotmg.f        | value   |
|-----------|----------------|-----------------|---------|
| `GAM`     | `4096.E0`      | `4096.D0`       | 2^12    |
| `GAMSQ`   | `1.67772E7`    | `16777216.D0`   | 2^24    |
| `RGAMSQ`  | `5.96046E-8`   | `5.9604645D-8`  | 2^−24   |

`GAMSQ` and `RGAMSQ` differ beyond the 7-digit single-precision
mantissa: `16777216 = 2^24` exactly, whereas `1.67772E7` drops the
trailing `16`; similarly `5.9604645E-8 = 2^-24` while `5.96046E-8`
truncates earlier. After migration to `KIND=16`, the literals in the
S source still carry only single-precision accuracy, so the migrated
`QROTMG` produced from `srotmg.f` has slightly coarser constants than
the one produced from `drotmg.f`. The D-sourced version is retained.


## LAPACK convergence divergences

LAPACK has **431 divergent pairs** out of 1018 migrated routines.
All are genuine upstream source-level differences between the
S/C and D/Z halves of LAPACK 3.12.1 — not migrator bugs. The D/Z
version is retained on disk as canonical for each pair.

The convergence test auto-filters precision-prefix local-variable
drift (S↔D, C↔Z swaps in unclassified local names) via
`_filter_precision_drift`, so these 431 are structural differences
only.

### 1. `QROUNDUP_LWORK` asymmetry (235 files)

`SROUNDUP_LWORK(N)` (in `INSTALL/sroundup_lwork.f`) converts an
INTEGER workspace size to a REAL such that `INT(result) >= N`. It
exists because LAPACK returns the optimal workspace size in
`WORK(1)`, which has the same floating-point type as the matrix
data. A 32-bit REAL cannot represent every INTEGER exactly:

```
REAL(16777217)           = 16777216.0    ! truncates down — unsafe
REAL(16777217)*(1+eps)   = 16777218.0    ! rounds up — safe
```

If a caller casts `WORK(1)` back to INTEGER and allocates that many
elements, a truncated-down value yields a buffer too small.
`SROUNDUP_LWORK` multiplies by `1+eps` to force rounding up.

**Asymmetry.** The S- and C-precision routines call
`WORK(1) = SROUNDUP_LWORK(LWKOPT)`. The D- and Z-precision routines
simply write `WORK(1) = LWKOPT`. Reason: 64-bit REAL has a 53-bit
mantissa, so every INTEGER up to `2^53 ≈ 9×10^15` round-trips through
REAL exactly. No workspace will ever approach that bound, so the
double-precision side has no bug to fix and skips the helper call
entirely. Only `dgesdd.f` and `zgesdd.f` use `DROUNDUP_LWORK`
(defined alongside, and equivalent under D precision).

This is a deliberate choice by upstream LAPACK 3.12: the release
notes describe `SROUNDUP_LWORK` as fixing "a subtle bug with
returning LWORK as a Float" that only manifests in single precision,
and the D/Z halves were intentionally left bare. After migration to
KIND=16 the migrated `QROUNDUP_LWORK` is harmless (128-bit REAL has
a 113-bit mantissa, so the round-up branch never triggers), but the
two halves remain textually distinct and are flagged as divergences.
The D/Z version is retained.

### 2. Explicit `REAL(int_expr, KIND=16)` casts (83 files)

The S/C halves wrap INTEGER expressions in explicit `REAL(...)`
casts before combining with a REAL operand; the D/Z halves rely
on implicit INTEGER→DOUBLE promotion.  Examples:

```fortran
! cgebal (S/C):  SCALE(L) = REAL(I, KIND=16)
! zgebal (D/Z):  SCALE(L) = I

! cgbrfs (S/C):  SAFE1 = REAL(NZ, KIND=16) * SAFMIN
! zgbrfs (D/Z):  SAFE1 = NZ * SAFMIN

! cbdsqr (S/C):  THRESH = MAX(TOL*SMINOA, REAL(MAXITR, KIND=16) * ...)
! zbdsqr (D/Z):  THRESH = MAX(TOL*SMINOA, MAXITR * ...)
```

Post-migration to KIND=16 the two are semantically identical but
textually distinct. The migrator does not strip these — doing so
would require type inference to prove the cast is a no-op, and the
cast is not formally wrong.

### 3. `ILAPREC` precision-character argument (10 files)

`ILAPREC(c)` (`SRC/ilaprec.f`) is a character→integer dispatcher
that maps a single-letter selector to one of the BLAST XBLAS
intermediate-precision enum constants:

```
'S' → BLAS_PREC_SINGLE     = 211
'D' → BLAS_PREC_DOUBLE     = 212
'I' → BLAS_PREC_INDIGENOUS = 213
'X'/'E' → BLAS_PREC_EXTRA  = 214
```

The value is passed to the XBLAS `BLAS_xGEMV_X` / `BLAS_xGEMV2_X`
calls inside the extended iterative-refinement drivers
(`[sdcz]la_*rfsx_extended.f`) to choose the residual precision.

**Asymmetry.** The single/complex drivers select `'D'` (one rung
above their base type):

```fortran
! sgerfsx.f  (and sgbrfsx / sporfsx / ssyrfsx,
!             cgerfsx / cgbrfsx / cporfsx / csyrfsx / cherfsx)
PREC_TYPE = ILAPREC( 'D' )         ! BLAS_PREC_DOUBLE = 212
```

The double/complex-double drivers select `'E'`:

```fortran
! dgerfsx.f  (and dgbrfsx / dporfsx / dsyrfsx,
!             zgerfsx / zgbrfsx / zporfsx / zsyrfsx / zherfsx)
PREC_TYPE = ILAPREC( 'E' )         ! BLAS_PREC_EXTRA = 214
```

Upstream picks whichever level is one rung above the base precision
of the driver — single uses double, double uses extra. The two
halves therefore request different XBLAS residual-precision modes
and the character literals (`'D'` vs `'E'`) survive into the
migrated output unchanged.

Diverging pairs: `sgbrfsx/dgbrfsx`, `sgerfsx/dgerfsx`,
`sporfsx/dporfsx`, `ssyrfsx/dsyrfsx`, `sgbsvxx/dgbsvxx`,
`cgbrfsx/zgbrfsx`, `cgerfsx/zgerfsx`, `cporfsx/zporfsx`,
`csyrfsx/zsyrfsx`, `cherfsx/zherfsx`.

These files also exhibit declaration-ordering drift (LOGICAL LSAME,
RCEQU positions differ between halves) and DZTHRESH/RTHRESH
parameter grouping differences.

After migration to `KIND=16` both halves are physically quad
precision; the character literal divergence is genuine upstream
source drift that the migrator cannot (and should not) paper over.
The D/Z version is retained on disk.

### 4. Other source-level drift (~103 files)

These are miscellaneous places where the S/C and D/Z halves of a
routine were maintained independently and drifted apart in
token-level but semantically-benign ways. They group into roughly
six patterns:

#### 4a. DO-loop / GO TO label renumbering (10 files)

The loop-label values were chosen differently between halves even
though the control structure is identical:

```fortran
! strsyl.f                  ! dtrsyl.f
DO 70 L = 1, N              DO 60 L = 1, N
  IF( L.LT.LNEXT ) GO TO 70   IF( L.LT.LNEXT ) GO TO 60
```

Files: `clatbs/zlatbs` (+58), `cptts2/zptts2` (+54),
`clatrs/zlatrs` (+52), `clatps/zlatps` (+48), `slatrs/dlatrs`
(+32), `slatbs/dlatbs` (+30), `slatps/dlatps` (+28),
`clapmt/zlapmt` (+22), `slapmt/dlapmt` (+22), `cpttrf/zpttrf`
(+16).

#### 4b. Unused `PARAMETER(ONE=…)` declaration on one side only (25 files)

One half declares a constant that the code path never references:

```fortran
! cgehd2.f declares ONE; zgehd2.f does not (or vice versa):
COMPLEX ONE
PARAMETER(ONE=(1,0))
```

Files: `cgehd2`, `cgelq2`, `cgeql2`, `cgeqr2`, `cgeqr2p`,
`cgeqrt2`, `cgeqrt3`, `cgerq2`, `cpftri`, `cunm2l`, `cunm2r`,
`cunml2`, `cunmr2`, `cupmtr` (complex side);
`sgehd2`, `sgelq2`, `sgeql2`, `sgeqr2`, `sgeqr2p`, `sgeqrt2`,
`sgeqrt3`, `sgerq2`, `sorm2l`, `sorm2r`, `sorml2` (real side).
`clahef_rk/zlahef_rk` splits a combined
`PARAMETER(ONE=…,ZERO=…)` into two single-parameter declarations.

#### 4c. `INTRINSIC` list mismatch (5 files)

One half lists more intrinsics than the other, even though only the
common subset is actually called:

```fortran
! cggbak.f      INTRINSIC MAX
! zggbak.f      INTRINSIC INT, MAX

! slansf.f      INTRINSIC SQRT, ABS
! dlansf.f      INTRINSIC SQRT, ABS, MAX
```

Files: `cggbak/zggbak`, `clanht/zlanht`, `cpoequb/zpoequb`,
`sggbak/dggbak`, `slansf/dlansf`.

#### 4d. `CMPLX(…, KIND=16)` cast on one side only (6 files)

C/Z halves use `CMPLX(expr, KIND=16)` to construct complex values
where the S/D halves assign directly:

```fortran
! chbevd_2stage:   WORK(1) = CMPLX(LWMIN, KIND=16)
! zhbevd_2stage:   WORK(1) = LWMIN

! clalsd:          B(J,I) = CMPLX(RWORK(RE), RWORK(IM), KIND=16)
! zlalsd:          (different code path — no RWORK splitting)
```

Files: `cgedmdq/zgedmdq`, `chbevd_2stage/zhbevd_2stage`,
`cheevd_2stage/zheevd_2stage`, `chgeqz/zhgeqz`,
`clalsd/zlalsd`, `csytri2/zsytri2`.

#### 4e. Type-declaration style differences (9 files)

The `[sc]laqz*` routines use bare `REAL,INTENT(…)` and
`COMPLEX,INTENT(…)` type-specs while the `[dz]laqz*` counterparts
use `REAL(KIND=16),INTENT(…)`. Also, `cgsvj1` has `IMPLICIT NONE`
while `zgsvj1` does not, and `sisnan`/`disnan` have different
function-result variable names.

Files: `claqz0/zlaqz0`, `claqz1/zlaqz1`, `claqz2/zlaqz2`,
`claqz3/zlaqz3`, `slaqz1/dlaqz1`, `slaqz2/dlaqz2`,
`slaqz3/dlaqz3`, `slaqz4/dlaqz4`, `cgsvj1/zgsvj1`,
`sisnan/disnan`, `slaisnan/dlaisnan`.

#### 4f. Formatting / expression-level differences

Semantically identical expressions written slightly differently:

```fortran
! cheequb:  C0 = (-(T*DI)*DI+2*WORK(I)*DI-N*AVG)   ! outer parens
! zheequb:  C0 = -(T*DI)*DI+2*WORK(I)*DI-N*AVG

! ctrsyl3:  DWORK(2,1) = (2*NBB+NBA)
! ztrsyl3:  DWORK(2,1) = 2*NBB+NBA

! clarf1f:  C(2,1)
! zlarf1f:  C(1+1,1)                               ! unfolded constant

! chgeqz:   WORK(1) = 1
! zhgeqz:   WORK(1) = DCMPLX(1)                    ! explicit cast
```

#### 4g. Algorithmic drift — different code paths

A handful of pairs use genuinely different implementations:

```fortran
! slasq2:  IEEE = .FALSE.                          ! hardcoded
! dlasq2:  IEEE = (ILAENV(10,'DLASQ2','N',1,2,3,4).EQ.1)

! clarf1f calls CGER   (single-complex GER, no conjugate)
! zlarf1f calls ZGERC  (conjugate GER) — different BLAS kernel

! clahef_rk uses CGEMM
! zlahef_rk uses ZGEMMTR                           ! different routine
```

The migrator cannot (and should not) close these — they look
symmetric under our character-based canonicalizer only by accident,
but they do call different BLAS kernels.

#### 4h. Workspace-query dummy-argument asymmetry

`sorcsd/dorcsd` (+13): one half allocates a tiny
`REAL DUMMY(1)` to pass to lazy workspace-query calls; the other
passes the real output arrays (which the query-mode call never
touches). Both are correct.

#### Case study: `clahef.f` vs `zlahef.f` — loop-bound off-by-one

In the Bunch-Kaufman pivot-unwinding loops (both in the upper-
triangular and lower-triangular branches), the loop-continuation
check differs between the two halves:

| Source          | Upper branch (label 60) | Lower branch (label 120) |
|-----------------|-------------------------|--------------------------|
| `clahef.f` (C)  | `IF( J.LE.N ) GO TO 60` | `IF( J.GE.1 ) GO TO 120` |
| `zlahef.f` (Z)  | `IF( J.LT.N ) GO TO 60` | `IF( J.GT.1 ) GO TO 120` |

The inner `[CZ]SWAP` call is guarded by `J.LE.N` in **both** halves,
and its body starts with an unconditional `J = J + 1`. When the
C-side does its extra `J == N` iteration, the inner increment makes
`J == N+1`, the swap guard becomes false, and the iteration is a
no-op. The same holds symmetrically for the `J == 1` case in the
lower branch. So the behavioural difference is a benign off-by-one:
clahef executes one additional empty loop turn that zlahef skips.

This is genuine upstream drift — the single-complex and
double-complex halves of `LAHEF` diverged at some point and were
never resynced. The migrator cannot (and should not) paper over it;
the Z-sourced canonical is kept on disk.


## Post-migration `converge` results

The `converge` subcommand is the authoritative post-migration
verifier: it reads each D/Z canonical off disk, re-migrates its
S/C sibling in memory, and compares the two under a light
normalizer (whitespace/case, `END KEYWORD` merging, sorted
`INTRINSIC`/`EXTERNAL` and simple-type-spec bare-identifier
lists, precision-neutral declaration sorting). Precision-prefix
local-variable drift (S↔D, C↔Z swaps in unclassified locals)
is auto-filtered by `_filter_precision_drift`. Every entry that
survives is either migrator slop or genuine upstream drift.

| library    | diverged | missing | character                                  |
|------------|----------|---------|--------------------------------------------|
| BLAS       | 3        | 0       | upstream S≠D asymmetry (see above)         |
| LAPACK     | 431      | 0       | SROUNDUP_LWORK, ILAPREC, upstream drift    |
| BLACS      | 0        | 0       | clean                                      |
| PBLAS      | 61       | 0       | K&R vs ANSI function style + editorial drift |
| ScaLAPACK  | 24       | 0       | upstream S≠D / C≠Z drift                   |

### BLAS (3 diverged)

After auto-filtering precision-prefix drift, only 3 genuine
upstream asymmetries remain — the `sdot`/`ddot`, `scasum`/`dzasum`,
and `srotmg`/`drotmg` cases documented in detail above.

### LAPACK (431 diverged)

All 431 are upstream S/C vs D/Z asymmetries, documented in the
sections above. Categories (a file may belong to multiple):

| category                                   | files |
|--------------------------------------------|-------|
| `QROUNDUP_LWORK` asymmetry (§1)            | 235   |
| Explicit `REAL(int, KIND=16)` casts (§2)   | 83    |
| `ILAPREC` + rfsx declaration drift (§3)    | 10    |
| `PARAMETER(ONE=…)` on one side only (§4b)  | 25    |
| DO-loop / GO TO label renumbering (§4a)    | 10    |
| `CMPLX(…, KIND=16)` cast (§4d)            | 6     |
| `INTRINSIC` list mismatch (§4c)            | 5     |
| Type-declaration style (§4e)               | 9     |
| Formatting / expression drift (§4f)        | ~20   |
| Algorithmic drift (§4g)                    | ~10   |
| Workspace-query dummy-arg (§4h)            | 1     |
| Other misc (see remaining files below)     | ~50   |

Many files overlap categories (e.g., 33 files have both
`QROUNDUP_LWORK` and `REAL(int, KIND=16)` casts). The 431 total
counts each divergent pair once.

### BLACS (0 diverged)

All C-language co-family pairs converge perfectly after migration.

### PBLAS (61 diverged)

All 61 divergences are due to **K&R-style vs ANSI-style function
declarations** in the C sources.  The S/C halves retain K&R-style
parameter declarations while the D/Z halves use ANSI prototypes:

```c
// S half (K&R style):
void pqagemv_ ( TRANS, M, N, ALPHA, ... )
F_CHAR_T TRANS ;
Int * IA , * INCX , ... ;

// D half (ANSI style):
void pqagemv_ ( F_CHAR_T TRANS, Int * M, Int * N, ... )
```

Additionally, `psamax_.c` vs `pdamax_.c` (+12) has editorial
differences: the S half uses direct array indexing while the D half
uses the `Mptr` macro.

This is upstream source asymmetry in ScaLAPACK's PBLAS C sources.

### ScaLAPACK (24 diverged)

24 Fortran pairs diverge due to upstream S≠D or C≠Z asymmetries.
Full file-by-file listing:

#### Precision-dependent algorithmic constants (6 files)

```fortran
! slarre2/dlarre2 and slarre2a/dlarre2a:
!   PERT=4.0 vs PERT=8.0; RTL=HNDRD*EPS vs RTL=SQRT(EPS)
! slarrf2/dlarrf2:
!   LSIGMA-=|LSIGMA|*TWO*EPS vs FOUR*EPS
! sstegr2/dstegr2, sstegr2a/dstegr2a, sstegr2b/dstegr2b:
!   MINRGP=3.0E-3 vs MINRGP=1.0E-3
```

These are deliberate precision-dependent tuning constants.

#### `PB_TOPGET` vs `PB_TOPSET` (2 files)

`pzungql/pcungql` and `pzunml2/pcunml2`: the C half calls
`PB_TOPGET` (query) while the Z half calls `PB_TOPSET` (modify).
Handled via `prefer_source` recipe directive.

#### Variable naming: `TTOPH`/`TTOPV` vs `CONJTOPH`/`CONJTOPV` (2 files)

`pssyttrd/pdsyttrd` (real): uses `TTOPH`/`TTOPV`.
`pzhettrd/pchettrd` (complex): uses `CONJTOPH`/`CONJTOPV` with
`CONJG()`. The names differ but the precision-drift filter doesn't
catch this because it's a semantic rename, not a prefix swap.

#### `PSLAIECT` vs `PDLAIECTB`/`PDLAIECTL` — different algorithm (1 file)

`psstebz/pdstebz` (+46): the S half calls a single `PSLAIECT`
helper for eigenvalue counting; the D half uses two separate
helpers `PDLAIECTB`/`PDLAIECTL` with an additional `IEFLAG`
dispatch. This is the largest single ScaLAPACK divergence and
reflects genuinely different algorithmic implementations for
single vs double precision.

#### Declaration grouping and array-size differences (7 files)

```fortran
! pslacon/pdlacon:  separate ESTWORK(1), TEMP(1), WORK(2) vs combined
! pslawil/pdlawil:  separate BUF(4), H11(1), H12(1)... vs combined
! pslaqr3/pdlaqr3:  TZROWS/TZCOLS workspace calc present in S, absent in D
! pzlawil/pclawil:  declaration order + TAULOC vs TAULOC(1) in args
! pzpotf2/pcpotf2:  XXDOTC helper call vs inline XDOTC function
! pzunmbr/pcunmbr:  PCHK1MAT vs PCHK2MAT external declaration
! pzlarz/pclarz:    TAULOC vs TAULOC(1) in subroutine arguments
```

#### Miscellaneous (8 files)

- `bslaexc/bdlaexc` (+2): `PARAMETER(TEN=10.0E0_16)` vs `TEN=1.0E1_16`
  — literal representation of the same value.
- `pslaed3/pdlaed3` (+4): `IINFO=0` vs `INFO=0`; extra `IF(I+J.LE.N)`
  guard in D half.
- `pslaqr3/pdlaqr3` (+10): D half has `MPI_WTIME` external, S half
  computes `LWK8` via `TZROWS*TZCOLS`, literal `0.0E0_16` vs `0.0E00_16`.
- `pssyevd/pdsyevd` (+2): workspace query `LWORK.EQ.-1 .OR. LIWORK.EQ.-1`
  vs just `LWORK.EQ.-1`.
- `pstrsen/pdtrsen` (+4): `IF/CALL/ENDIF` block vs single-line `IF…CALL`.
- `pzheevd/pcheevd` (+2): `DESCZ,11` vs `DESCZ,12` — argument-position
  constant in PCHK2MAT call.
- `pzheevx/pcheevx` (+2) and `pzhegvx/pchegvx` (+2): duplicate `REAL`
  in INTRINSIC list on Z side.
- `pzlarzc/pclarzc` (+2): `TAULOC` vs `TAULOC(1)` in XGEBR2D call.

### Migrator fixes driven by converge

Running `converge` across BLAS and LAPACK surfaced three bugs in
the migrator that had been masked by the heavy canonicalizer used
for `migrate`/`diverge`:

1. **Rename-pattern cache collision** — `_RENAME_PATTERN_CACHE`
   was keyed by `id(rename_map)`. When `_migrate_with_flang`
   built a fresh per-file filtered rename map, Python reused the
   freed dict's address, so subsequent files hit stale cached
   patterns that didn't include their routines' names. Symptom:
   `qtrsv.f` on disk still carried `SUBROUTINE DTRSV` etc.
   BLAS migrate divergences dropped from 33 to 3 after re-keying
   on `frozenset(upper_map.items())`.

2. **Bare complex literal canonicalization** — C-prefix sources
   wrote `(0.0, 0.0)` (default-kind) while Z-prefix sources wrote
   `(0.0d0, 0.0d0)`; only the latter picked up the `KIND` suffix
   from the existing D/E-exponent rule. `replace_literals` was
   extended to rewrite bare `N.M` to `N.ME0_{kind}`, excluding
   string literals (FORMAT specs) and masking Fortran logical
   operators (`.EQ.`/`.AND.`/…) so `.EQ.0.0` is never re-parsed
   as `.EQ` + `.0` + `.0`. `E+0` exponents also normalize to `E0`.

3. **`CMPLX`/`DCMPLX` intrinsic rewriting** — `DCMPLX(...)` was
   left unrenamed when its inner arguments already contained a
   `KIND=16` from an earlier `DBLE`→`REAL` rewrite (the guard
   against double-annotating was checking for any `KIND=` in
   inner, not just top-level). And bare `CMPLX(a, b)` without a
   `KIND` was passing through unchanged, leaving default-kind
   complex values where kind-16 was intended. Both are fixed:
   `CMPLX` is now in the intrinsic map with `needs_kind=True`
   and the nested-`KIND=` guard uses paren-depth tracking.
