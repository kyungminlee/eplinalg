# Upstream bugs: ScaLAPACK 2.2.3

*Last catalogued: 2026-05-11. See `UPSTREAM_BUGS.md` for the cross-library
index, audit methodology, bug summary table, and how fixes are carried.*

This file collects ScaLAPACK 2.2.3 bugs in the vendored
`external/scalapack-2.2.3/` source (including PBLAS under
`PBLAS/SRC/`). LAPACK and MUMPS bugs are catalogued in
`UPSTREAM_BUGS-LAPACK.md` and `UPSTREAM_BUGS-MUMPS.md` respectively.
Bugs that span both LAPACK and ScaLAPACK (e.g. the XERBLA-string
sweep) live in the index file.

## 2026-05-11 ScaLAPACK residual-divergence audit

Per-pair audit of the 26 ScaLAPACK divergent pairs (see
`doc/scalapack-residual-divergence-categorization.md` for the full
categorization) surfaced **four real upstream bugs**, each with its
own detailed section below:

- `pzungql.f` — PB_TOPGET should be PB_TOPSET (BLACS topology restore)
- `pzunml2.f` — same PB_TOPGET/PB_TOPSET typo
- `pdsyevd.f` — LQUERY missing `.OR. LIWORK.EQ.-1` branch
- `pslaed3.f` — missing `INFO=0` init + missing bounds-guard on INDROW/INDCOL writes

Pre-patch the recipe carried `prefer_source: PCUNGQL, PCUNML2` to
route migrated output around the Z-half topology bugs, plus
`expected_divergences:` entries for all four pairs. The patches
retire both `prefer_source` pins and four whitelist entries. The
`pdsyevd` patch is the only one that changes migrated output (qsyevd
previously inherited the LIWORK-query bug from the D-canonical half);
the other three were already routed-around but are fixed at source
for convergence and future Netlib upstreaming.

Divergent-pair count: **26 → 22** post-patch.

---

### ScaLAPACK 2.2.3: `pzungql.f` / `pzunml2.f` fail to restore BLACS broadcast topology

**Symptom.** After PZUNGQL or PZUNML2 returns, the BLACS process-grid
broadcast topology for the routine's ICTXT is left in the temporary
`'I-ring'` (rowwise) / `' '` (columnwise) state that the routine
installed for its internal communication. The caller's previously
configured topology is silently overwritten. Subsequent BLACS
broadcasts on the same context use the unexpected topology, which
can change reduction order (affecting roundoff) and, on some
broadcast shapes, deadlock when the caller's protocol assumes a
different ring direction.

**Root cause.** Both routines follow the standard save-temporary-
restore idiom for the broadcast topology — save the user's current
topology with `PB_TOPGET`, install `'I-ring'` for the routine's own
broadcasts with `PB_TOPSET`, and restore the user's via `PB_TOPSET`
at the end. The Z-half copy-pastes `PB_TOPGET` into the *restore*
position instead of `PB_TOPSET`. The C-half counterparts have the
correct call.

`pzungql.f:245-251`:

```fortran
*     Save the topologies, install the routine's working topology.
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Rowwise',    ROWBTOP )
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise',    'I-ring' )
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', ' ' )
```

`pzungql.f:292-293` (the bug — should be PB_TOPSET to restore):

```fortran
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Rowwise',    ROWBTOP )
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
```

`pcungql.f:292-293` (correct):

```fortran
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise',    ROWBTOP )
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
```

`pzunml2.f` has the same bug at lines 394-395 with identical
surrounding context. `pcunml2.f:394-395` is correct.

**Affected files.**
- `external/scalapack-2.2.3/SRC/pzungql.f` lines 292-293.
- `external/scalapack-2.2.3/SRC/pzunml2.f` lines 394-395.

**Fix.** Change both `PB_TOPGET` calls at each site to `PB_TOPSET`.
The `ROWBTOP` / `COLBTOP` arguments are already the saved values from
lines 247-248 / 332-333; the routine just needs to push them back
out via `PB_TOPSET` instead of re-reading the current (still-
temporary) values via `PB_TOPGET`.

**Patches.** `recipes/scalapack/patches/pzungql.f.patch`,
`recipes/scalapack/patches/pzunml2.f.patch`.

**Severity.** Caller-visible side effect: callers that rely on a
specific broadcast topology after PZUNGQL / PZUNML2 returns get
silently the wrong topology. May cause perf regression, change
roundoff order, or deadlock depending on the caller's BLACS
protocol. Not a memory-safety bug.

**Why upstream's tests miss it.** ScaLAPACK's test driver constructs
a fresh context per test case and tears it down after, so the
topology leak isn't visible: the broken context is destroyed before
anything else uses it. The bug surfaces only in long-running
applications that reuse a context across many ScaLAPACK calls.

**Upstream report.** File as a PR at
https://github.com/Reference-ScaLAPACK/scalapack — both halves are
trivial three-character fixes, and the C-half is the reference for
correctness.

---

### ScaLAPACK 2.2.3: `pdsyevd.f` LQUERY recognizes only LWORK=-1, not LIWORK=-1

**Symptom.** Calling `PDSYEVD` with `LIWORK = -1` and a positive
LWORK (the LAPACK workspace-query convention asking *just* for the
integer-workspace size) returns INFO = -16 instead of populating
`IWORK(1) = LIWMIN` and returning INFO = 0. The S-half (`pssyevd`)
handles the call correctly. Migrated `qsyevd` previously inherited
the bug because D is the canonical default in the migrator.

**Root cause.** `pdsyevd.f:225` defines the workspace-query predicate
recognizing only LWORK = -1:

```fortran
            LQUERY = ( LWORK.EQ.-1 )
```

The subsequent validation at `pdsyevd.f:253-255` uses `.NOT.LQUERY`
as the guard for both LWORK and LIWORK too-small checks:

```fortran
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -14
            ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -16
```

With LIWORK = -1, the second branch evaluates `LIWORK.LT.LIWMIN`
(true, since -1 < any positive LIWMIN) `.AND. .NOT.LQUERY` (true,
since LQUERY was false). INFO = -16 is set, and the routine returns
via the PXERBLA path without populating `WORK(1)` / `IWORK(1)`.

`pssyevd.f:224` has the correct predicate:

```fortran
            LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
```

**Affected files.**
- `external/scalapack-2.2.3/SRC/pdsyevd.f` line 225.

**Fix.** Replace `LQUERY = ( LWORK.EQ.-1 )` with `LQUERY = (
LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )`. This matches the S-half and the
LAPACK convention: either workspace argument set to -1 indicates a
size-query call.

**Patch.** `recipes/scalapack/patches/pdsyevd.f.patch`.

**Severity.** API contract violation. Standard LAPACK convention is
that *either* LWORK = -1 *or* LIWORK = -1 signals a workspace query.
Callers using the LIWORK-only query form (common when the caller
already knows LWORK is sufficient and only wants to right-size IWORK)
get a confusing INFO = -16 error. This is the same upstream-weakness
family as the documented `pzheevd.f / pcheevd.f miss LIWORK
validation` bug above — both are LIWORK-handling gaps in the
heevd / syevd family, but they fail in opposite directions (pzheevd
admits too-small LIWORK silently; pdsyevd rejects the legitimate
LIWORK = -1 query).

**Why upstream's tests miss it.** Upstream test drivers exercise
workspace queries via LWORK = -1, not LIWORK = -1 alone. The S-half
correctness was likely a happy accident of one developer writing
the inclusive form by habit.

**Upstream report.** File as a PR at
https://github.com/Reference-ScaLAPACK/scalapack — the fix is a
single-line change and the S-half provides the reference for the
correct predicate form.

---

### ScaLAPACK 2.2.3: `pslaed3.f` two bugs — uninitialized INFO + missing bounds guard

**Symptom.** Two independent bugs in the S-half copy of the secular-
equation solver:

1. The output parameter `INFO` is never initialized on the success
   path. Callers that read `INFO == 0` to mean "success" instead
   see whatever stack-garbage value was on entry. The D-half
   (`pdlaed3.f`) correctly initializes `INFO = 0`.

2. The block-cyclic distribution loop writes `INDROW(I+J)` and
   `INDCOL(I+J)` without a bounds check. When N is not a multiple
   of NB, the final outer-loop iteration's inner loop overruns
   the INDROW / INDCOL arrays by `(I + NB - 1) - N` entries. The
   D-half guards with `IF (I+J.LE.N)`.

**Root cause 1 — uninitialized INFO.** `pslaed3.f:156`:

```fortran
*     Test the input parameters.
*
      IINFO = 0
*
*     Quick return if possible
*
      IF( K.EQ.0 )
     $   RETURN
```

The initialization is on the *local* variable `IINFO` (used as the
out-arg of the `SLAED4` calls at lines 209 / 316). The *output
parameter* `INFO` is declared at line 12 but only written when
`IINFO.NE.0` at lines 211 / 318. On the success path INFO is never
written.

`pdlaed3.f:156` correctly initializes the output parameter:

```fortran
      INFO = 0
```

(Note: pdlaed3 doesn't initialize IINFO at all — the redundant
`IINFO = 0` in pslaed3 is dead code from the same edit that
introduced the typo.)

**Root cause 2 — missing bounds guard.** `pslaed3.f:166-173`:

```fortran
      ROW = DROW
      COL = DCOL
      DO 20 I = 1, N, NB
         DO 10 J = 0, NB - 1
            INDROW( I+J ) = ROW
            INDCOL( I+J ) = COL
   10    CONTINUE
         ROW = MOD( ROW+1, NPROW )
         COL = MOD( COL+1, NPCOL )
   20 CONTINUE
```

When `N` is not a multiple of `NB`, the final outer iteration runs
with `I + (NB-1) > N`. The inner loop writes past the end of INDROW
/ INDCOL. The arrays are caller-supplied with documented size N, so
this is a memory-safety bug.

`pdlaed3.f:168-172` has the guard:

```fortran
         DO 10 J = 0, NB - 1
            IF( I+J.LE.N ) THEN
               INDROW( I+J ) = ROW
               INDCOL( I+J ) = COL
            END IF
   10    CONTINUE
```

**Affected files.**
- `external/scalapack-2.2.3/SRC/pslaed3.f` lines 156 and 168-171.

**Fix.** (a) Replace `IINFO = 0` with `INFO = 0` at line 156. (b)
Wrap the INDROW/INDCOL writes in `IF (I+J.LE.N) THEN ... END IF`.
Both changes copy the D-half verbatim.

**Patch.** `recipes/scalapack/patches/pslaed3.f.patch`.

**Severity.** (1) Uninitialized output parameter — undefined
behavior in caller-visible state, although in practice many callers
ignore INFO on the success path. Spec-violation but rarely caught.
(2) Out-of-bounds memory write — undefined behavior, memory
corruption when N % NB != 0. The crashes are non-deterministic
(depends on what follows INDROW / INDCOL in caller memory). May
surface as silent corruption of unrelated data.

**Why upstream's tests miss it.** Bug 1 is masked because most
callers don't check INFO on the success path, and stack memory
sometimes contains 0. Bug 2 is masked because ScaLAPACK's standard
test grid uses N values that are multiples of NB; the failure
condition `N % NB != 0` is outside the regular test matrix.

**Upstream report.** File as a PR at
https://github.com/Reference-ScaLAPACK/scalapack — both fixes are
mechanical (copy from D-half) and the D-half is the reference.


## 2026-05-12 B? deep-dive (ScaLAPACK portion)

The B? items remaining in the `2026-05-11` audit were re-examined; the
ScaLAPACK-side finding is below. See `UPSTREAM_BUGS-LAPACK.md` for the
LAPACK portion (`cgedmdq`, `cgedmd` / `zgedmd` divergences).

### Stays B?: `pslaqr3` vs `pdlaqr3` — LWK8 workspace formula divergence

**Symptom.** `pslaqr3.f` computes `LWK8 = 2*TZROWS*TZCOLS` at line
400, where `TZROWS` and `TZCOLS` are NUMROC computations of the
block-cyclic dimensions of an auxiliary matrix. `pdlaqr3.f` line
398 hardcodes `LWK8 = 0` and never computes TZROWS/TZCOLS — those
variables are declared at line 258 of both halves but only used in
the S-half. LWK8 feeds into the `LWKOPT = MAX(LWK1, ..., LWK8)`
that the routine returns as the optimal-workspace estimate.

**Possible interpretations.**
- pdlaqr3 under-reports LWKOPT — a bug; callers that allocate
  `LWORK = LWKOPT` may hit insufficient-workspace at runtime in
  inner routines.
- pslaqr3 over-reports LWKOPT — harmless conservatism; callers
  allocate extra unused memory.
- LWK8 represents workspace for an optional code path that's
  exercised in S but not D — neither over- nor under-reports.

**Additional finding.** `pdlaqr3.f:273` declares `MPI_WTIME` in its
EXTERNAL list but never calls it. Pure dead-EXTERNAL (D-class), not
related to the LWK8 question but worth noting.

**Affected files.**
- `external/scalapack-2.2.3/SRC/pslaqr3.f` lines 398-400 (LWK8
  computation present); line 258 (TZROWS/TZCOLS declarations,
  used).
- `external/scalapack-2.2.3/SRC/pdlaqr3.f` line 398 (LWK8=0
  hardcoded); line 258 (TZROWS/TZCOLS declared but never used);
  line 273 (dead `MPI_WTIME` in EXTERNAL list).

**Severity.** Depends on interpretation. If pdlaqr3 truly
under-reports LWKOPT, callers querying workspace and allocating
exactly that size hit runtime errors in inner routines — caller-
visible failure on workspace-sized usage. If LWK8 represents an
optional code path or pslaqr3 over-reports, the difference is
harmless. The bug-vs-cosmetic distinction cannot be made by static
reading.

**Why upstream's tests miss them.** ScaLAPACK's pslaqr3/pdlaqr3
test driver likely allocates workspace generously (not via
LWKOPT query) so any under-reporting in pdlaqr3 goes unnoticed.
The dead `TZROWS`/`TZCOLS`/`MPI_WTIME` declarations are pure
compile-time hygiene and never cause a test failure.

**Upstream report.** File as an issue on the Netlib ScaLAPACK
tracker. Resolution requires reading the LWK1..LWK8 algorithm and
matching each LWKi to a specific workspace consumer — which is
ScaLAPACK QR-iteration domain knowledge.

## 2026-05-11 symmetric-fix sweep (ScaLAPACK)

Every D/Z-half ScaLAPACK patch that previously stood alone has been
extended with its S/C-half sibling. Sibling patches generated by
`/tmp/gen_sibling_patch.py` (precision-letter swap derived from
upstream-file comparison, blacklisting `DOUBLE`/`DESCA`/`ZERO`-shape
false positives, then hand-tuned for files whose EXTERNAL list
ordering differs between halves). Patches added:

`bstrexc`, `cpttrsv`, `spttrsv`, `pclascl`, `pslascl`, `pcpttrs`,
`pspttrs`, `pcpbtrsv`, `pspbtrsv`, `psgeequ`, `pcgeequ`, `psposvx`,
`pcposvx`, `pslarz`, `pclarz`, `pslarzb`, `pclarzb`, `pclarzc`,
`pcheevd`.

The "Documented but not patched (S/C non-canonical)" subsections
below — under the INFO=-N, PCHK?MAT, and pzheevd LIWORK entries —
no longer apply: all of those S/C sibling files now carry their
hand-translated patches. `recipes/scalapack.yaml`'s
`asymmetric_patches:` list is correspondingly reduced to a single
entry (`pzunmbr.f.patch` — `pcunmbr.f` already has the fix upstream
so no sibling patch is needed).


## ScaLAPACK 2.2.3: more INFO=-N param-position bugs (extended sweep 2026-05-09)

A second pass on the ScaLAPACK INFO=-N sweep (using a stricter
ScaLAPACK-aware filter that ignores LAPACK's `INFO = -100*PARAM - INDEX`
descriptor convention and multi-arg tests) surfaced three additional
families of real bugs.

### Family A: `?pttrs` / `?pbtrsv` LWORK sentinel-check off-by-one

8 routines × 2 occurrences each = 16 INFO assignments where the LWORK
sentinel-validation `IF(LWORK.LT.-1)` and post-validation
`IF(LWORK.LT.WORK_SIZE_MIN) IF(LWORK.NE.-1)` set INFO to one less than
LWORK's actual signature position.

| Routine | LWORK pos | Old INFO | Fix |
|---------|----------:|---------:|---:|
| `pdpttrs.f`, `pspttrs.f` | 13 | -12 | -13 |
| `pcpttrs.f`, `pzpttrs.f` | 14 | -13 | -14 |
| `pdpbtrsv.f`, `pspbtrsv.f` | 15 | -14 | -15 |
| `pcpbtrsv.f`, `pzpbtrsv.f` | 15 | -14 | -15 |

INFO points to WORK (one before LWORK in the signature) instead of
LWORK itself. **Patched on D/Z** halves; S/C documented (non-canonical).

### Family B: `?lascl` CFROM/CTO INFO codes inherited from LAPACK

`pslascl.f`, `pdlascl.f`, `pclascl.f`, `pzlascl.f` all set:
- `INFO = -4` for bad CFROM (CFROM is at position **2**)
- `INFO = -5` for bad CTO (CTO is at position **3**)

The values -4, -5 are correct in LAPACK's `dlascl.f` because the
LAPACK signature is `(TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO)` —
KL/KU at 2/3 push CFROM to 4 and CTO to 5. ScaLAPACK's `pdlascl`
*drops* KL and KU, making the signature `(TYPE, CFROM, CTO, M, N, A,
IA, JA, DESCA, INFO)` — but kept the LAPACK INFO codes without
re-numbering. **Patched on D/Z** halves with `INFO=-4 → -2` and
`INFO=-5 → -3`. `?lascl` is auxiliary but used internally; bad
CFROM/CTO would report wrong arg.

### Family C: `?ormrz` / `?unmrz` duplicate K-validation (L unchecked)

All 4 halves (`psormrz`, `pdormrz`, `pcunmrz`, `pzunmrz`) have:

```fortran
ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
   INFO = -5
ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN          ! ← duplicate, copy-paste bug
   INFO = -6                                  ! claims to validate L (param 6)
```

The second `IF` was meant to validate L (parameter 6) but the test was
copy-pasted from K (parameter 5). Result: **L is never validated** —
caller passing illegal L (e.g. negative, or larger than M/N depending
on SIDE) slips past argument validation.

Compare LAPACK's `dormrz.f:253-255` which has the correct test:
```fortran
ELSE IF( L.LT.0 .OR. ( LEFT .AND. L.GT.M ) .OR.
$         ( .NOT.LEFT .AND. L.GT.N ) ) THEN
   INFO = -6
```

**Patched on D/Z** halves; restores the LAPACK-style L test. S/C
documented (non-canonical). `pdormrz` and `pzunmrz` already had source
overrides for the post-loop-condition bug from earlier work; the
L-validation fix is added on top.

### Severity

All three families: argument-validation gaps. Family A is purely
diagnostic (LWORK is still validated, just the wrong INFO is reported).
Family B is also diagnostic. **Family C is a real validation gap**:
illegal L values are accepted and the routine proceeds, potentially
performing operations on invalid data.

### Patched on every half (2026-05-11)

- `pspttrs`, `pcpttrs`, `pspbtrsv`, `pcpbtrsv` (LWORK INFO off-by-one)
- `pslascl`, `pclascl` (CFROM/CTO INFO -4/-5 → -2/-3)
- `psormrz`, `pcunmrz` (L never validated) — sibling patches were
  already present in the original 2026-05-09 batch.

### Upstream report

Not yet filed.

---

## ScaLAPACK 2.2.3: `pzheevd.f` / `pcheevd.f` miss LIWORK validation (sweep 2026-05-09)

**Symptom.** Caller passing too-small `LIWORK` (positive value, not
the `-1` query sentinel) bypasses argument validation. The routine
proceeds to the body, where `IWORK(1:LIWMIN)` is written — overrunning
the user's IWORK buffer.

**Root cause.** Argument-validation block in `pzheevd.f:255-289` and
`pcheevd.f:251-285` validates `LWORK` and `LRWORK` but lacks a
matching `LIWORK.LT.LIWMIN` check:

```fortran
ELSE IF( LWORK.LT.LWMIN .AND. LWORK.NE.-1 ) THEN
   INFO = -14
ELSE IF( LRWORK.LT.LRWMIN .AND. LRWORK.NE.-1 ) THEN
   INFO = -16
* No LIWORK check — should be:
* ELSE IF( LIWORK.LT.LIWMIN .AND. LIWORK.NE.-1 ) THEN
*    INFO = -18
ELSE IF( IROFFA.NE.0 ) THEN
   ...
```

**Discovered by** the 2026-05-09 LWORK-validation-gap mechanical sweep,
which audited every routine in LAPACK + ScaLAPACK SRC for `LWORK` /
`LIWORK` / `LRWORK` arguments without matching `.LT.` validation. Of
the 16 routines in the heev/syev/heevd/syevd/heevr/syevr/heevx/syevx
family, **only `pcheevd` and `pzheevd` are missing the LIWORK check**.
All real-precision halves (`pdsyevd`, `pssyevd`) and all eigensolver
variants in the family validate LIWORK correctly. This is a complex-
half-only asymmetry.

**Affected files.**
- `external/scalapack-2.2.3/SRC/pzheevd.f` (line 263 area).
- `external/scalapack-2.2.3/SRC/pcheevd.f` (line 263 area, same bug).

**Fix.** Insert the missing `ELSE IF (LIWORK.LT.LIWMIN .AND.
LIWORK.NE.-1) INFO = -18` branch between the LWORK and LRWORK
validations and `IROFFA.NE.0`. (LIWORK is at parameter position 18 in
both routines' signatures.)

**Patched on both halves** (2026-05-11): `pzheevd.f.patch` adds the
LIWORK check on the Z half; `pcheevd.f.patch` adds it on the C half.
The C-half patch consists of only the LIWORK validation hunk (the
PCHK2MAT 11→12 fix in `pzheevd.f.patch` is already present in
upstream `pcheevd.f`). The migrator's `prefer_source: PZHEEVD` pin
remains for the PCHK2MAT case; convergence is symmetric.

**Severity.** Argument-validation gap. Calling pzheevd with a too-
small LIWORK silently corrupts memory rather than returning
INFO=-18. Same severity class as the documented `pssyevd` LIWORK=-1
LQUERY gap (S03), which is a different aspect of the same upstream
weakness in the heevd / syevd family's workspace handling.

**Upstream report.** Not yet filed.

---

## ScaLAPACK 2.2.3: PCHK?MAT parameter-position drift (mechanical sweep 2026-05-08)

A mechanical sweep parsed every `CALL PCHK1MAT(...)` /
`CALL PCHK2MAT(...)` site (316 total in `scalapack-2.2.3/SRC/`),
extracted the `MAPOS0` / `NAPOS0` / `DESCAPOS0` / `MBPOS0` / `NBPOS0` /
`DESCBPOS0` literal arguments, and verified each one matches the
actual position of the corresponding M / N / DESC arg in the calling
routine's signature.

**53 mismatches in 11 routines, in 4 distinct upstream bug families:**

### Family 1: `?heevr` / `?syevr` — DESCB POS0 21 should be 19

Files: `pcheevr.f`, `pdsyevr.f`, `pssyevr.f`, `pzheevr.f` (all four
halves). Signature ends `..., Z, IZ, JZ, DESCZ, WORK, ...`, putting
DESCZ at position 19. The `PCHK2MAT(..., DESCZ, 21, ...)` call uses
21 — likely lifted from a sibling routine that has additional args
before DESCZ. On illegal DESCZ, PXERBLA reports parameter 21
(WORK or LWORK depending on signature) instead of 19 (DESCZ).

### Family 2: `?hegvx` / `?sygvx` — M/N POS0 4 should be 5

Files: `pchegvx.f`, `pdsygvx.f`, `pssygvx.f`, `pzhegvx.f` (all four).
Signature is `(IBTYPE, JOBZ, RANGE, UPLO, N, ...)` — N at position 5.
Caller forgot to count `IBTYPE` (it's a generalized eigenvalue problem
and IBTYPE selects A·x=λBx vs A·Bx=λx), and used 4 (which is UPLO).
Both PCHK2MAT and PCHK1MAT calls in each routine have the same off-by-1
on every M/N parameter index. On illegal N, PXERBLA misreports as
parameter 4 (UPLO).

### Family 3: `pdtrord` / `pstrord` — every POS0 off by +1

Signature is `(COMPQ, SELECT, PARA, N, T, IT, JT, DESCT, Q, IQ, JQ, DESCQ, ...)`
putting N at 4, DESCT at 8, DESCQ at 12. The PCHK?MAT calls use
5, 9, 13. Looks like a renamed-arg shift that wasn't propagated —
maybe an older signature with an extra leading parameter that was
later removed. Three call sites per routine.

### Family 4: `pzheevd` — DESCB POS0 11 should be 12

Already documented separately as the first PCHK2MAT bug found via
convergence audit (caught earlier this session, patched at
`recipes/scalapack/source_overrides/pzheevd.f`).

### Patches in migrated archive

D/Z-canonical halves patched at `recipes/scalapack/source_overrides/`:
- `pdsyevr.f` (DESCZ POS0 21→19)
- `pzheevr.f` (DESCZ POS0 21→19)
- `pdsygvx.f` (N POS0 4→5 in both calls; B/Z param indices left
  unchanged since they were already correct at 9 and 13)
- `pzhegvx.f` (same)
- `pdtrord.f` (N POS0 5→4, DESCT POS0 9→8, DESCQ POS0 13→12 in all
  three calls)

Wired in `recipes/scalapack.yaml` with matching `prefer_source` pins
(`PDSYEVR`, `PZHEEVR`, `PDSYGVX`, `PZHEGVX`, `PDTRORD`) so the
patched halves win convergence over the still-buggy C/S siblings.

### S/C halves

`pssyevr.f`, `pcheevr.f`, `pssygvx.f`, `pchegvx.f`, `pstrord.f` were
patched alongside their D/Z siblings in the original 2026-05-08
sweep (see `recipes/scalapack/patches/` — sibling files share the
same dataset of POS0 corrections). No further work in this entry.

### Severity

All diagnostic-only. PXERBLA still aborts and rejects bad input;
users just see the wrong parameter number in the error message
("parameter -19 had an illegal value" instead of "parameter -21",
or the inverse). No memory or numerical impact. Same severity class
as the XERBLA-string typos and the original `pzheevd` bug.

### Why upstream's tests miss it

PCHK?MAT parameter-position correctness is invisible unless a caller
deliberately passes an invalid DESC/M/N argument and inspects PXERBLA
output. Test drivers always pass valid args.

### Upstream report

Not yet filed. The `?heevr/?syevr` and `?hegvx/?sygvx` families are
self-contained: 8 file edits cover both bug clusters in single typo
form, making this an attractive low-effort upstream PR.

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

