# LAPACK 3.12.1 upstream nits (dead declarations)

*Cosmetic-only sibling of [`UPSTREAM_BUGS-LAPACK.md`](UPSTREAM_BUGS-LAPACK.md):
same 2026-05-11 audit (see
`doc/archive/lapack-residual-divergence-categorization.md`), same vendored
tree, but **no numerical, correctness, or interface impact** —
none of these are patched in-tree.*

Every entry below is a declaration in an upstream LAPACK source file
that is never referenced in the routine's body. The migration's
diverge comparator surfaces them because the co-family sibling
correctly omits the dead declaration — making the two halves diverge
on declarations alone. Listed as documentation of the
convergence-report noise floor; patching would be ~23 one-line
patches with zero functional effect.

## D-class entries

| File | Dead symbol | Declared at | Sibling state | Notes |
|---|---|---|---|---|
| `zlaqp2.f` | `CONE` (COMPLEX*16) | line 167-169 | claqp2 omits | Only `ONE` (real) is used; `CONE` never referenced |
| `zlaqp2rk.f` | `CONE` (COMPLEX*16) | line 366-368 | claqp2rk omits | Same pattern |
| `chetrf_aa_2stage.f` | `LWKOPT` (INTEGER) | line 183 | zhetrf_aa_2stage omits | Never assigned or read |
| `chesv_aa_2stage.f` | `ZERO`, `ONE` (COMPLEX) + `ILAENV_EP` (INTEGER) | lines 206-215 | zhesv_aa_2stage omits | Dead constants and unused EXTERNAL |
| `claswlq.f` | `ILAENV` (INTEGER func) | line 190 + EXTERNAL line 192 | zlaswlq omits | Declared but never called |
| `clarf1f.f` | `J` (INTEGER) | line 145 area | (sibling has separate `J` use) | zlarf1f has dead J; clarf1f correctly omits |
| `csytrf_aa_2stage.f` | `ZLASWP` in EXTERNAL | line 192-area | zsytrf_aa_2stage has it; csytrf_aa_2stage doesn't | The actual `CALL ZLASWP` is **commented out** in zsytrf body — both files have the dead-code, but the declaration leaked only in one |
| `zunbdb.f` | `ONE` (COMPLEX*16) | line 310-311 | cunbdb omits | `REALONE` (real) is used; `ONE` (complex) never referenced |
| `zunbdb1.f` | `ONE` (COMPLEX*16) | declaration block | cunbdb1 omits | Same pattern |
| `zunbdb2.f` | `ONE` (COMPLEX*16) | declaration block | cunbdb2 omits | Same |
| `zunbdb3.f` | `ONE` (COMPLEX*16) | declaration block | cunbdb3 omits | Same |
| `zunbdb4.f` | `ONE` (COMPLEX*16) | declaration block | cunbdb4 omits | Same |
| `sgebd2.f` | `ONE` (REAL) | line 205-206 | dgebd2 omits | Declared in PARAMETER block but never used in body |
| `sgelqt.f` | `SGEQRT2`, `SGEQRT3` in EXTERNAL | line 143 | dgelqt omits | sgelqt's body calls only `SGELQT3` and `SLARFB` — the other two LQT/QRT helpers are dead |
| `slaswlq.f` | `SGEQRT`, `STPQRT` in EXTERNAL | EXTERNAL list | dlaswlq omits | Same pattern — dead helpers in declaration only |
| `slasyf_rk.f` | `JB`, `JJ` (INTEGER) | line 286 | dlasyf_rk omits | INTEGER decl-list includes JB / JJ but the body never uses them |
| `dorbdb.f` | `ONE` (DOUBLE PRECISION) | declaration block | sorbdb omits | (Sibling-asymmetric: D declares, S omits) |
| `dorbdb1.f` | `ONE` (DOUBLE PRECISION) | declaration block | sorbdb1 omits | Same |
| `dorbdb2.f` | `ONE` (DOUBLE PRECISION) | declaration block | sorbdb2 omits | Same |
| `dorbdb3.f` | `ONE` (DOUBLE PRECISION) | declaration block | sorbdb3 omits | Same |
| `dorbdb4.f` | `ONE` (DOUBLE PRECISION) | declaration block | sorbdb4 omits | Same |
| `dsterf.f` | `RMAX` (DOUBLE PRECISION) | line 112 | ssterf omits | Declared, **assigned** at line 148 (`RMAX = DLAMCH('O')`), but never read afterward |
| `ssysv_aa.f` | `ILAENV_EP` (INTEGER func) | declaration block | dsysv_aa omits | Declared but never called |
| `zgetf2.f` | `SFMIN` (DOUBLE PRECISION) + `DLAMCH` external | line 129, 133, 167 | cgetf2 omits both | See discussion below; safety lives in ZRSCL |

## zgetf2 SFMIN — a D-class case that initially looked like a bug

**zgetf2.f line 129 + 167**: declares `DOUBLE PRECISION SFMIN`,
assigns `SFMIN = DLAMCH('S')`, then never reads it. cgetf2 (the
sibling) correctly omits both the declaration and the assignment.

This was initially flagged as a potential missing-safety-check bug
because sgetf2 / dgetf2 use SFMIN to gate the pivot scaling at
line 184:

```fortran
IF( ABS(A(J,J)) .GE. SFMIN ) THEN
   CALL SSCAL( M-J, ONE / A(J,J), A(J+1,J), 1 )    ! fast reciprocal
ELSE
   DO I = 1, M-J                                    ! safe per-element
      A(J+I,J) = A(J+I,J) / A(J,J)
   END DO
END IF
```

— making it look like zgetf2 has the SFMIN declaration but lost
the IF branch around it. But on closer reading, both **cgetf2 and
zgetf2** delegate the pivot scaling to a dedicated routine
(`CRSCL` / `ZRSCL`, "reciprocal scaling") which internally handles
the underflow-safe split. The S/D halves inline what C/Z do via a
helper call. So:

- The safety check **is** present on every half — just structured
  differently between the real and complex variants.
- zgetf2's `SFMIN = DLAMCH('S')` is genuinely dead code, left over
  from a refactor that moved the safe-scaling logic into ZRSCL but
  didn't remove the now-unused threshold value.
- cgetf2 was cleaned up correctly.

Categorized as a D-class dead declaration alongside the others —
no UPSTREAM_BUGS entry needed.

## Aggregate

- **Total dead-declaration cases**: 24
- **Aggregate cosmetic impact**: each contributes ~1-3 lines to the
  per-pair diff in the no-whitelist divergence report
- **Patch cost if cleanup is wanted**: ~23 one-line patches
  (`patches/<file>.f.patch`) with hunk header `@@ -<DECL>,1 +<DECL>,1 @@`
  removing the dead identifier from the declaration list, or
  removing entire dead PARAMETER lines

## Why not patch them

LAPACK upstream-evolution typically corrects these on a case-by-case
basis when surrounding code changes. Patching them in-tree would:

1. Add patch maintenance for zero functional gain
2. Create churn each time upstream bumps the source files
3. Not reduce divergence-report value — the *pair* still appears in
   the report only if it had non-cosmetic divergence too, and a
   one-line cosmetic diff in an otherwise-clean pair is easy to
   recognise

The comparer-side alternative (teach `_canonicalize_for_compare` to
strip orphan PARAMETER/EXTERNAL declarations) is also unattractive:
it would require parsing the routine body to know what's actually
referenced, which is much more work than the alleviation justifies.

## See also

- `doc/UPSTREAM_BUGS.md` — real bugs (numerical, correctness, or
  interface mismatch) found in the same audit
- `doc/archive/lapack-residual-divergence-categorization.md` — full per-pair
  classification of all 110 remaining divergent pairs
