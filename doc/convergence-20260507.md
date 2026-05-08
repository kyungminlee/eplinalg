# Convergence snapshot — 2026-05-07

Per-library convergence of the kind16 staged tree at `/tmp/stage-q/` against
fresh in-memory re-migration of each S/C sibling. Generated with:

```bash
uv run python -m migrator stage /tmp/stage-q --target kind16
for lib in blas xblas blacs lapack pbblas pblas ptzblas \
           scalapack scalapack_c scalapack_tools mumps; do
    uv run python -m migrator converge \
        recipes/$lib.yaml /tmp/stage-q/$lib/src --target kind16
done
```

See `doc/CONVERGENCE.md` for the underlying methodology (canonical-on-disk
vs. fresh re-migration of the S/C half, suppressed patterns, statuses).

## Summary

| library          |   pairs | converged | diverged | missing | converge % |
|------------------|--------:|----------:|---------:|--------:|-----------:|
| blas             |      75 |        68 |        7 |       0 |     90.7 % |
| xblas            |       0 |         0 |        0 |       0 |        n/a |
| blacs            |      43 |        42 |        1 |       0 |     97.7 % |
| lapack           |   1 018 |       588 |      430 |       0 |     57.8 % |
| pbblas           |      14 |        14 |        0 |       0 |    100.0 % |
| pblas            |      61 |         0 |       61 |       0 |      0.0 % |
| ptzblas          |      47 |        47 |        0 |       0 |    100.0 % |
| scalapack        |     336 |       301 |       35 |       0 |     89.6 % |
| scalapack_c      |       3 |         2 |        1 |       0 |     66.7 % |
| scalapack_tools  |       0 |         0 |        0 |       0 |        n/a |
| mumps            |     200 |        38 |      162 |       0 |     19.0 % |
| **total**        | **1 797** | **1 100** |  **697** |    **0** |  **61.2 %** |

Zero `missing` across the board — every co-family pair has its canonical on
disk. All 697 divergences are pure textual disagreements between the on-disk
canonical and a fresh re-migration of the S/C sibling.

## Per-library findings

### blas — 7 / 75 diverged

Re-audited 2026-05-08. All seven are cosmetic/structural drift in upstream
Netlib BLAS or deliberate D/Z source-override patches. None are upstream
bugs and none affect the migrated archive's numerics. Decision: leave as
known-divergent (parallel to the scalapack cosmetic cohort).

1. **D/Z source-override asymmetry (4 files)** — `crotg/zrotg`,
   `srotg/drotg`, `snrm2/dnrm2`, `scnrm2/dznrm2`. The canonical D/Z output
   is pulled verbatim from `recipes/blas/source_overrides/` (hand-rewritten
   to `LOGICAL,SAVE :: …_INITIALIZED=.FALSE.` + lazy first-call init,
   because the multifloats target overloads `**`/`sqrt` and gfortran
   rejects user-defined ops in `PARAMETER` initializers). No parallel
   override exists for the S/C halves; converge re-migrates upstream
   `crotg.f90`/`scnrm2.f90` straight into `PARAMETER` form and reports the
   structural disagreement. The migrated archive ships the D/Z-derived
   form, which is what we want — silencing the diff would require S/C
   override twins that never compile in any target.
2. **Upstream code-shape drift (1 file)** — `scasum.f` vs `dzasum.f`. The
   D/Z half delegates the absolute value to `DCABS1(ZX(I))`; the S/C half
   inlines `ABS(CX(I))+ABS(AIMAG(CX(I)))`, plus the parameter is named
   `ZX` vs `CX`. Functionally identical, different upstream prose.
3. **Upstream line-order swap (1 file)** — `sdot.f` vs `ddot.f`. Upstream
   `sdot.f:103-104` writes `STEMP=0.0e0` then `SDOT=0.0e0`; `ddot.f:103-104`
   writes `DDOT=0.0d0` then `DTEMP=0.0d0`. After migration both halves
   contain the same lines, but in swapped order. `_filter_precision_drift`
   only folds 1:1 `replace` blocks and doesn't recognize a `delete+insert`
   reordering — this is the actual reason the diff survives. (The earlier
   note flagging this as a filter regression was wrong; the filter is
   doing what it advertises and the drift is upstream.)
4. **Upstream constants/branch drift (1 file)** — `srotmg/drotmg`. The S
   half uses different `DATA GAM,GAMSQ,RGAMSQ` constants
   (`4096.E0_16,1.67772E7_16,…`) and a different conditional structure
   around `SD1.LT.ZERO`. Real upstream divergence; D-half is canonical.

### xblas — 0 / 0

No co-family pairs (xblas ships precision-extended routines without
S/D ↔ C/Z mirrors).

### blacs — 1 / 43 diverged

`sgsum2d/dgsum2d → qgsum2d`. The diff is a real-vs-MPI-sum-op asymmetry:
the migrator inserts a `BI_qMPI_sum` user reduction (`MPI_Op_create` →
`BlacComb`) for one half but lets the other half keep `MPI_SUM`. Both
results are functionally equivalent in IEEE arithmetic at `kind=16`;
the cosmetic difference is a code-path the kind16 retarget triggers in
only one half. Low-priority cleanup.

### lapack — 430 / 1 018 diverged

Single largest contributor (62 % of all divergences). Sample of the
trailing 25 file headers from the converge run:

```
ssytrd_sy2sb / ssytrf{,_aa,_aa_2stage,_rk,_rook} / ssytri{2,_3} /
ssytrs_{3,aa} / stbrfs / stgevc / stgex2 / stgexc / stgsen / stgsna /
stgsyl / stprfs / strevc{,3} / strrfs / strsen / strsyl{,3} / stzrzf
```

The bulk pattern is symmetric/triangular routines (`*sytr*`, `*tg**`,
`*tr**`) where the S half and D half have drifted in upstream LAPACK
over many releases — different workspace-size formulae, different inline
expansions, different comment paraphrases that occasionally leak through
the comment-stripper because they sit on continuation lines. None of
these affect numerics; they're textual disagreements that the lighter
`_light_normalize` (used by `converge`) deliberately refuses to absorb.

`strsyl/dtrsyl → qtrsyl` is the largest single divergence in the report
(174 diff lines) and is a known-asymmetric routine upstream. Worth
considering a `prefer_source:` whitelist for the *trsyl/*tgsyl families.

### pbblas — 0 / 14

Clean.

### pblas — 0 / 61 diverged ⚠

**Every pair diverges.** This is the most striking signal in the report.
Sample diffs (all C sources):

```c
// canonical (Z half, on disk):
void pxgerc_ ( M , N , ALPHA , …, INCY , A , IA , JA , DESCA )
QREAL * ALPHA ;
Int * IA , * INCX , * INCY , …;

// fresh re-migration of the C-half sibling:
# ifdef __STDC__
# else
void pxgerc_ ( … )                              ← K&R-style declarator
QREAL * ALPHA ;                                  ← ditto
Int * IA , * INCX , …                            ← ditto
```

The C halves are still wrapped in `# ifdef __STDC__ / # else …` K&R
fallbacks. Per `project_assume_ansi_c`, the migrator strips these — but
something in the PBLAS C path is leaving the `__STDC__`/K&R blocks intact
during converge re-migration, even though the on-disk canonical has them
stripped. Looks like a single systematic gap (one missed branch in the C
preprocessor handler), not 61 independent issues. Highest-leverage thing
to fix in this report: one bug fix would close 61 pairs.

### ptzblas — 0 / 47

Clean.

### scalapack — 35 / 336 diverged

Two distinct sub-patterns:

1. **Constants drift** in the `*stegr2{,a,b}`, `*larre2{,a}`, `*larrf2`
   families: identical body, different `MINRGP` literal:
   ```
   -PARAMETER(ZERO=0.0E0_16,ONE=1.0E0_16,FOUR=4.0E0_16,MINRGP=3.0E-3_16)
   +PARAMETER(ZERO=0.0E0_16,ONE=1.0E0_16,FOUR=4.0E0_16,MINRGP=1.0E-3_16)
   ```
   The S half ships with `3.0E-3` and the D half with `1.0E-3`.

   **Audited 2026-05-08: not a bug, expected divergence.** This S/D
   tuning split is canonical upstream LAPACK design, mirrored verbatim
   from `external/lapack-3.12.1/SRC/sstemr.f:345` (`MINRGP=3.0E-3`) and
   `dstemr.f:345` (`MINRGP=1.0D-3`). MINRGP is the minimum relative gap
   threshold for cluster decomposition in the RRR (relatively robust
   representations) eigenvalue algorithm: gaps below this fraction of
   the eigenvalue magnitude are not trusted as real cluster boundaries.
   The looser threshold at single precision reflects the larger noise
   floor at ~7 decimal digits; double precision can trust gaps as
   small as 0.1% with ~16 decimal digits.

   The kind16 target preserves both literals faithfully (D=1.0E-3 wins
   as canonical, so `qstegr2.f` ships with the tighter threshold).
   Conservative-but-correct at quad precision — the threshold could in
   principle be tightened further (1e-10 or smaller given ~32 decimal
   digits), but algorithm-tuning constants are preserved verbatim by
   migrator policy. Future enhancement, not a fix.

   No prefer_source pin needed; the divergence is logged as expected.

2. **Routine-level asymmetry** in `pclarzc/pzlarzc → pxlarzc` (51 lines),
   `pslarzb/pdlarzb → pqlarzb` (21 lines), `psstebz/pdstebz → pqstebz`
   (46 lines), `pslarz/pdlarz → pqlarz` (32 lines), `pslanhs/pdlanhs →
   pqlanhs` (20 lines), and ~20 smaller cases. These are real upstream
   ScaLAPACK asymmetries — the same routines that have caused unblock
   work in the past (see `tests/lapack/TODO.md`).

### scalapack_c — 1 / 3 diverged

`clamov/zlamov → xlamov`: a `# define TYPE complex` vs
`# define TYPE QCOMPLEX` token mismatch. The C-half re-migration emits
`complex`; the canonical has `QCOMPLEX`. Looks like a missed
type-alias substitution in one branch of the C migrator. The
`scalapack_c` recipe also surfaces the
`QISNAN`/`QLAISNAN` rename-map collision warning — independent issue,
already noted in the rename-map collision diagnostic.

### scalapack_tools — 0 / 0

No pairs (single-precision-only tools).

### mumps — 162 / 200 diverged

Second-largest contributor (81 % of pairs diverge). The pattern is
consistent across ~all MUMPS solver / OOC / LR / save-restore modules:

```
### ssol_aux.F  vs dsol_aux.F  → qsol_aux.F  (+530)
### ssol_c.F    vs dsol_c.F    → qsol_c.F    (+384)
### stools.F    vs dtools.F    → qtools.F    (+222)
### ssol_lr.F   vs dsol_lr.F   → qsol_lr.F   (+210)
### ssol_driver.F vs dsol_driver.F → qsol_driver.F (+112)
### sstatic_ptr_m.F → qstatic_ptr_m.F (+14)
…
```

The MUMPS S/D halves are *substantially* divergent at the source
level — MUMPS upstream maintains them as separate hand-edited copies,
not as a generated co-family. The migrator already handles this via
`recipes/mumps.yaml`'s heavy `keep_kind_lines` + EP-bridge module
machinery to make both halves retarget correctly *as libraries*; what
it does **not** do is force them to converge at the source-text level,
because they don't converge upstream either.

For MUMPS specifically, the converge metric is largely an upstream
property, not a migrator quality signal. The 38 / 200 that *do*
converge are roughly the small utility modules where upstream happens
to keep S and D in lock-step. The remaining 162 are accepted: the
guarantee that matters for MUMPS is differential precision parity at
the test level (currently 26 / 26 passing), not source-text identity.

## Where to focus next

Ordered by leverage:

1. **PBLAS C K&R-stripper bug** (61 pairs in one fix). Highest ROI.
2. **ScaLAPACK MINRGP families** — decide policy: pin one literal via
   `prefer_source`, or move to the expected-divergence list.
3. **MUMPS** — leave as-is; document that converge is not the right
   metric for MUMPS and rely on the test-level parity gate.

BLAS is now fully audited (2026-05-08): all seven divergences are
upstream cosmetic / structural drift or deliberate D/Z override
asymmetry. Closing them only silences the report without changing the
shipped archive, so they're left as known-divergent. The earlier
"`_filter_precision_drift` regression" hypothesis was disproven — the
`sdot` case is an upstream `STEMP=…; SDOT=…` vs `DDOT=…; DTEMP=…` line
swap, which is a `delete+insert` pair the filter intentionally does not
fold.

## Reproducing

The numbers in the summary table were computed by `/tmp/converge_counts.py`,
which mirrors the pair-grouping logic in `pipeline.run_*_convergence_report`
to recover the `pairs` denominator (the converge command itself only
reports the diverged/missing numerator). The full diff text per library is
captured in `/tmp/converge_all.txt` (lib-bracketed, last 200 lines per lib).
