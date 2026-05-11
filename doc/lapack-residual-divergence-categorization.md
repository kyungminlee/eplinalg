# LAPACK residual divergence categorization (post-audit, 110 pairs)

Per-pair classification of every divergent S↔D or C↔Z pair remaining
after the symmetric-fix sweep and bug-fix audit.

**Categories**

- **C** — Cosmetic: line-position swap, decl-order swap, equivalent expression
- **W** — Comparer/whitespace artifact (no real source difference, fixable in normalizer)
- **D** — Dead declaration (unused variable/EXTERNAL/PARAMETER on one side)
- **N** — Algorithmic-tuning split (intentional precision-specific constant)
- **R** — Refactor (structural difference, functionally equivalent)
- **B** — Real bug worth patching
- **B?** — Ambiguous / minor inconsistency

| Pair | Cat | Note |
|---|---|---|
| cgedmd.f90 vs zgedmd.f90 | B? | Z passes `'R'` to XGEJSV's NTYPE arg (transposed conjugate variant) where C passes `'N'`. Different algorithmic choice for complex JSV; needs domain review. Also `OFL = LAMCH('O')*LAMCH('P')` (C, more conservative) vs `LAMCH('O')` (Z) |
| cgedmdq.f90 vs zgedmdq.f90 | R | Workspace-query structure (LZWORK in LQUERY, separate LIWORK call) |
| cgejsv.f vs zgejsv.f | C | `MAX(OPTWRK,MINWRK)` vs `MAX(MINWRK,OPTWRK)` — arg order |
| cgemqr.f vs zgemqr.f | C | `LWORK.LT.MAX(1,LW)` vs `LWORK.LT.LWMIN` — equivalent forms |
| cgerfsx.f vs zgerfsx.f | C+N | Decl reorder + intentional `ILAPREC('D')` vs `ILAPREC('E')` (S-precision pairs with double for refinement; D/Z pairs with extra) |
| cgesvd.f vs zgesvd.f | C | Decl reorder + `MAX(MINWRK,MAXWRK)` arg-order swap |
| cgesvdx.f vs zgesvdx.f | W | `WORK(1)=MAXWRK` vs `WORK(1)=CMPLX(MAXWRK,ZERO)`. Equivalent. Plus `MINMN.EQ.0` vs `M.EQ.0.OR.N.EQ.0`. Comparer could canonicalize. |
| cgetf2.f vs zgetf2.f | D | C declares `SFMIN` and computes `SFMIN=SLAMCH('S')`. Z doesn't declare/use. (cgetf2 actually uses SFMIN for a scaling check that zgetf2 omits — may actually be a missing safety check on Z) |
| cggbal.f vs zggbal.f | C | `ONE` (named) vs `1` (literal) — Z lacks the ONE param |
| cggev3.f vs zggev3.f | C | C hardcodes `'N','N'` in CGGHD3 call (workspace query in ILV=false branch where JOBVL=JOBVR='N'). Equivalent. |
| cgsvj1.f vs zgsvj1.f | C | Z has `IMPLICIT NONE`, C doesn't. Cosmetic |
| cheevx.f vs zheevx.f | R | Initial-LWKOPT structural difference (already partially patched) |
| cherfsx.f vs zherfsx.f | C+N | Same as cgerfsx |
| chesv_aa_2stage.f vs zhesv_aa_2stage.f | D | Z declares ZERO/ONE PARAMETERs and ILAENV_EP, C doesn't (dead on Z) |
| chetrd_hb2st.F vs zhetrd_hb2st.F | C | Upstream local-var naming asymmetry: SIZEV (Z) vs SICEV (C) — encoded in name with 2nd precision letter |
| chetrf_aa_2stage.f vs zhetrf_aa_2stage.f | D | Z declares unused `LWKOPT`, C doesn't |
| chetri_3.f vs zhetri_3.f | C | Line-position swap of `WORK(1)=LWKOPT` (C wraps in N==0 guard, Z doesn't) |
| chetrs_aa.f vs zhetrs_aa.f | R | `DO WHILE` (C) vs `DO K=...` (Z) loop style |
| clahef_rk.f vs zlahef_rk.f | C | Dead PARAMETER reorder for ZERO/ONE |
| clamswlq.f vs zlamswlq.f | C | Line-position swap of XGEMLQT call |
| claqp2.f vs zlaqp2.f | D | Z declares unused `CONE` PARAMETER |
| claqp2rk.f vs zlaqp2rk.f | D | Z declares unused `CONE` PARAMETER |
| claqz2.f vs zlaqz2.f | W | `REAL :: ` vs `REAL ::` whitespace around `::` |
| clarf1f.f vs zlarf1f.f | D+C | Z declares unused J; `C(2,1)` vs `C(1+1,1)` (arithmetic-equivalent) |
| clartg.f90 vs zlartg.f90 | C | Z spells `XZERO=>XZERO` (no-op rename) explicitly |
| claswlq.f vs zlaswlq.f | D | C declares unused `ILAENV_EP` (dead) |
| clatrs3.f vs zlatrs3.f | C | Line-position swap |
| claunhr_col_getrfnp2.f vs zlaunhr_col_getrfnp2.f | W | `DOUBLE PRECISION CABS1` (Z) vs `REAL CABS1` (C) — comparer should normalize |
| cporfsx.f vs zporfsx.f | C+N | Same as cgerfsx |
| crscl.f vs zrscl.f | D | Z patches added ZSCAL to EXTERNAL; C side hasn't been patched here. **Could be a B** but `CSCAL` already in `EXTERNAL` upstream — Z's added QSCAL is `DSCAL→QSCAL`, real-precision scale not present in C |
| cstemr.f vs zstemr.f | N | MINRGP=3e-3 (C) vs 1e-3 (Z); RTOL1 formula split. Documented algorithmic tuning. |
| csyrfsx.f vs zsyrfsx.f | C+N | Same as cgerfsx |
| csysv_aa_2stage.f vs zsysv_aa_2stage.f | C | Line-position swap of RETURN |
| csyswapr.f vs zsyswapr.f | C | `A(LDA,N)` vs `A(LDA,*)` — equivalent dummy-arg declarations |
| csytrf_aa_2stage.f vs zsytrf_aa_2stage.f | D+C | Z declares unused ZLASWP (commented out call). Paren: `I.EQ.J-1` vs `I.EQ.(J-1)` |
| ctgsen.f vs ztgsen.f | C | EXTERNAL line layout (QLAMCH on same line vs separate) |
| cunbdb.f vs zunbdb.f | D | Z declares unused `ONE` PARAMETER |
| cunbdb1.f vs zunbdb1.f | D | Z declares unused ONE PARAMETER + parameter-order |
| cunbdb2.f vs zunbdb2.f | D | Z declares unused ONE PARAMETER |
| cunbdb3.f vs zunbdb3.f | D | Z declares unused ONE PARAMETER |
| cunbdb4.f vs zunbdb4.f | D | Z declares unused ONE PARAMETER |
| cuncsd.f vs zuncsd.f | C | `LWORKOPT=MAX(...); WORK(1)=LWORKOPT` vs `WORK(1)=MAX(...)` — same result |
| cuncsd2by1.f vs zuncsd2by1.f | C | `DUM(1)` vs `DUM` in BBCSD call — equivalent for assumed-size |
| cungr2.f vs zungr2.f | C | Line-position swap of `A(...)=ONE` |
| cunmlq.f vs zunmlq.f | R | Structural: C uses IF(N=0...) LWKOPT=1 ELSE ENDIF + WORK(1)=LWKOPT; Z uses inline `WORK(1)=1` for the N=0 case |
| cunmrz.f vs zunmrz.f | C | `LWORK.LT.NW` vs `LWORK.LT.MAX(1,NW)` — equivalent since NW≥1 by construction |
| sbdsdc.f vs dbdsdc.f | N | `EPS` (S) vs `0.9*EPS` (D) — tighter convergence tolerance on D |
| sgbrfsx.f vs dgbrfsx.f | C+N | Same as cgerfsx |
| sgbsvxx.f vs dgbsvxx.f | C | `INTEGER` decl order swap (KL/KU position) |
| sgebak.f vs dgebak.f | C | Decl-order swap |
| sgebd2.f vs dgebd2.f | D | S declares unused `ONE`, D doesn't |
| sgedmd.f90 vs dgedmd.f90 | B? | S passes `-1` (query) for LIWORK to GESVDQ; D passes `LIWORK` (actual). `LWRSVQ=INT(RDUMMY(1))` (S) vs `MAX(MWRSVQ, INT(...))` (D). Likely intentional design split, but the MAX wrap suggests D has a safety bump |
| sgedmdq.f90 vs dgedmdq.f90 | B? | S misses `.OR.WNTVCQ` in IF condition; uses `MIN(M,N)` instead of `MINMN`. Could be missing branch on S |
| sgejsv.f vs dgejsv.f | B? | `'L'` (S) vs `'G'` (D) GESVJ first arg at line 1711. Real algorithmic choice. Needs domain expert |
| sgelqt.f vs dgelqt.f | D | S declares unused SGEQRT2, SGEQRT3 in EXTERNAL |
| sgelsd.f vs dgelsd.f | R | Workspace-query refactor (S uses IF block + RETURN; D uses linear + GOTO 1) |
| sgelss.f vs dgelss.f | **B** | Patched on `fix-sgelss-workspace` branch (issue L10). S uses ILAENV estimate for SGELQF where D queries SGELQF directly. PR open. |
| sgeqrfp.f vs dgeqrfp.f | C | Line-position swap of K=MIN(M,N) |
| sgerfsx.f vs dgerfsx.f | C+N | Same as cgerfsx |
| sgesvd.f vs dgesvd.f | C | Line-position swap of MAXWRK |
| sgesvdq.f vs dgesvdq.f | C | Variable rename `LWUNLQ` (S) ↔ `LWORLQ` (D); decl-list reorder |
| sgesvx.f vs dgesvx.f | C | Line-position swap of `WORK(1)=RPVGRW` |
| sgetsls.f vs dgetsls.f | C | Line-position swap of `WSIZEM=1` |
| sgges.f vs dgges.f | **B**+W | **S missing XERBLA from EXTERNAL** (declares calls but not symbol). Plus `IF(...)THEN` vs `IF(...) THEN` whitespace (comparer should normalize) |
| sgges3.f vs dgges3.f | C+W | Line-position swap + `IF(...)THEN` whitespace |
| sggev3.f vs dggev3.f | R | Structural workspace-query reorder (already patched for XERBLA bug) |
| sggevx.f vs dggevx.f | R | BALANC validation rewrite + workspace-min computation rewrite (functionally equivalent) |
| sggqrf.f vs dggqrf.f | C | `LWKOPT=MAX(...); WORK(1)=LWKOPT` vs `WORK(1)=MAX(...)` |
| sggrqf.f vs dggrqf.f | C | Same as sggqrf |
| slaexc.f vs dlaexc.f | C | `4 INFO=1` (S, label on same line) vs `4 CONTINUE / INFO=1` (D) — equivalent |
| slaln2.f vs dlaln2.f | C | Decl-order swap |
| slamswlq.f vs dlamswlq.f | C | Line-position swap of `TR=1` |
| slaqr2.f vs dlaqr2.f | C+W | BETA-save stylistic (already-patched SLARF1L typo). Residual is cosmetic |
| slaqr3.f vs dlaqr3.f | C | BETA-save stylistic |
| slaqz0.f vs dlaqz0.f | C | EXTERNAL list ordering: QLASET position |
| slarf1f.f vs dlarf1f.f | C | `C(2,1)` (S) vs `C(1+1,1)` (D) — arithmetic-equivalent |
| slarf1l.f vs dlarf1l.f | N | D has FIRSTV optimization (algorithmic enhancement on D, not present on S) — kept whitelisted |
| slarre.f vs dlarre.f | N | PERT=4/8, BSRTOL formula, RTL formula — algorithmic tuning splits |
| slarrf.f vs dlarrf.f | N | TWO*EPS vs FOUR*EPS multiplier — tuning split |
| slasq2.f vs dlasq2.f | **B?** | S has `IEEE=.FALSE.` (hardcoded) where D queries `IEEE = ILAENV_EP(10,...).EQ.1`. **Real bug**: S unconditionally disables IEEE-aware path |
| slaswlq.f vs dlaswlq.f | D | S declares unused SGEQRT, STPQRT in EXTERNAL |
| slasyf_rk.f vs dlasyf_rk.f | D | S declares unused JB, JJ in INTEGER list |
| slatrs3.f vs dlatrs3.f | C | Line-position swap of RETURN |
| sopmtr.f vs dopmtr.f | C | Already-patched (DLARF→DLARF1F); residual is AII-save cosmetic |
| sorbdb.f vs dorbdb.f | D | D declares unused `ONE` PARAMETER |
| sorbdb1.f vs dorbdb1.f | D | D declares unused `ONE` PARAMETER |
| sorbdb2.f vs dorbdb2.f | D | D declares unused `ONE` PARAMETER |
| sorbdb3.f vs dorbdb3.f | D | D declares unused `ONE` PARAMETER |
| sorbdb4.f vs dorbdb4.f | D | D declares unused `ONE` PARAMETER (and the I01 LWORK fix was patched earlier) |
| sorcsd.f vs dorcsd.f | **B?** | S passes `DUMMY(1)` (1-element scratch) to inner ORGQR/ORGLQ/ORBDB/BBCSD calls where D passes the actual arrays (`U1`, `THETA`, etc.). S's pattern is wasteful — workspace-query gets bogus argument hints. May affect optimal LWORK estimation |
| sorg2l.f vs dorg2l.f | C | Line-position swap of A(...)=ONE |
| sorgr2.f vs dorgr2.f | C | Line-position swap of A(...)=ONE |
| sormr2.f vs dormr2.f | C | Already-patched cosmetic residual |
| sormrz.f vs dormrz.f | C | Line-position swap of WORK(1)=1 |
| sporfsx.f vs dporfsx.f | C+N | Same as cgerfsx |
| sstemr.f vs dstemr.f | N | MINRGP + RTOL1 tuning split (matches cstemr) |
| ssterf.f vs dsterf.f | D | D declares unused RMAX variable and assigns it once (RMAX=LAMCH('O')) but never reads it |
| sstevr.f vs dstevr.f | C | Decl-order + `ITMP1` swap routine inlined in S, omitted in D (functionally moot if not called) |
| ssyevr.f vs dsyevr.f | R | LWORK validation structure + TEST flag pattern (S uses TEST var, D inlines) + `WORK(1)=26` vs `WORK(1)=1` for N==1 (workspace-return convention — 3 of 4 family use 1) |
| ssyevr_2stage.f vs dsyevr_2stage.f | C | Same TEST pattern as ssyevr |
| ssyrfsx.f vs dsyrfsx.f | C+N | Same as cgerfsx |
| ssysv_aa.f vs dsysv_aa.f | D+W | S declares unused ILAENV_EP. Plus `'SYSV_AA'` (no trailing space, S) vs `'SYSV_AA '` (trailing space, D) — comparer should canonicalize XERBLA banner whitespace |
| ssytrd_sb2st.F vs dsytrd_sb2st.F | C | Local-var name asymmetry SISEV (S) vs DIDEV (D) — known upstream stylistic |
| ssytrd_sy2sb.f vs dsytrd_sy2sb.f | W | `''` (empty string, S) vs `' '` (single space, D) for ILAENV2STAGE OPTS arg — both valid; comparer-fixable |
| ssytrf_aa_2stage.f vs dsytrf_aa_2stage.f | C | `I.EQ.(J-1)` vs `I.EQ.J-1` paren style |
| ssytri2.f vs dsytri2.f | **B?** | S asks ILAENV for SSYTRF blocksize; D asks for DSYTRI2. The routine is SYTRI2; either is plausible but D is more self-consistent |
| ssytrs_3.f vs dsytrs_3.f | C | `DO K=1,N,1` vs `DO K=1,N` — explicit-stride redundancy |
| ssytrs_aa.f vs dsytrs_aa.f | R | Same DO-WHILE vs DO-loop refactor as chetrs_aa |
| stgsyl.f vs dtgsyl.f | C | S has redundant `.AND.NOTRAN` inside an `IF(NOTRAN)` block |
| strsyl3.f vs dtrsyl3.f | C | Line-position swap of inner assignment |

## Summary

| Category | Count | Action |
|---|---:|---|
| **C** (cosmetic) | 60 | No action |
| **D** (dead declaration) | 22 | No action; could cleanup if desired |
| **N** (algorithmic tuning) | 9 | No action — intentional design |
| **R** (refactor, equivalent) | 9 | No action |
| **W** (comparer-fixable) | 7 | Could extend comparer |
| **B?** (ambiguous bug) | 6 | Domain review needed |
| **B** (real bug) | 2 | sgges (missing XERBLA), sgelsd L10 (already patched on fix branch) |

## Real bugs identified for patching

1. **sgges**: missing `XERBLA` in EXTERNAL list (calls XERBLA at error path but doesn't declare it). Mirror of sggev3 fix.

2. **slasq2**: hardcoded `IEEE = .FALSE.` where dlasq2 dynamically queries `ILAENV(10, 'DLASQ2', 'N', ...) .EQ. 1`. This forces single-precision QR to disable IEEE-aware shifts even on IEEE-754 hardware. Actual performance/correctness impact.

## Comparer-fixable normalizations (W cases)

1. `XERBLA('NAME')` vs `XERBLA('NAME ')` — strip trailing-space from XERBLA string literal
2. `''` vs `' '` for ILAENV2STAGE OPTS — collapse empty/single-space string literals
3. `DOUBLE PRECISION X` vs `REAL X` after migration — comparer should equate
4. `IF(...)THEN` vs `IF(...) THEN` — whitespace around `THEN`
5. `CMPLX(MAXWRK, ZERO)` vs `MAXWRK` — explicit cast-to-complex of real value
