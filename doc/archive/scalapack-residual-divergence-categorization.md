# ScaLAPACK residual-divergence categorization (2026-05-11)

Per-pair audit of the 26 residual divergent pairs in
`recipes/scalapack.yaml` after the symmetric-fix sweep (see
`doc/UPSTREAM_BUGS.md` §2026-05-11 symmetric-fix sweep).

Reproducer:

```
uv run python -m migrator diverge recipes/scalapack.yaml --target kind16 --no-whitelist --context 0
```

Audit categories: **B** real bug / **B?** ambiguous / **W**
comparer-fixable / **C** cosmetic / **N** algorithmic tuning.

## Summary

| Pair | Class | Fixed by |
|---|---|---|
| `pzungql/pcungql` | B | `pzungql.f.patch` |
| `pzunml2/pcunml2` | B | `pzunml2.f.patch` |
| `pdsyevd/pssyevd` | B | `pdsyevd.f.patch` |
| `pslaed3/pdlaed3` | B (two bugs in one file) | `pslaed3.f.patch` |
| `pslaqr3/pdlaqr3` | B? | (deferred — needs domain review) |
| `bslaexc/bdlaexc` | W | (deferred — comparer) |
| `psgesvd/pdgesvd` | W | (deferred — comparer) |
| `pzgesvd/pcgesvd` | W | (deferred — comparer) |
| `pzheev/pcheev` | W | (deferred — comparer) |
| `pclarz/pzlarz` | C | — |
| `pclarzc/pzlarzc` | C | — |
| `pzlawil/pclawil` | C | — |
| `pslacon/pdlacon` | C | — |
| `pslawil/pdlawil` | C | — |
| `pssyttrd/pdsyttrd` | C | — |
| `pzhettrd/pchettrd` | C | — |
| `pzlattrs/pclattrs` | C | — |
| `pstrsen/pdtrsen` | C | — |
| `pzpotf2/pcpotf2` | C | — |
| `psstebz/pdstebz` | C | — |
| `slarre2/dlarre2` | N | — |
| `slarre2a/dlarre2a` | N | — |
| `slarrf2/dlarrf2` | N | — |
| `sstegr2/dstegr2` | N | — |
| `sstegr2a/dstegr2a` | N | — |
| `sstegr2b/dstegr2b` | N | — |

**26 total: 4 B (fixed), 1 B? (deferred), 4 W, 11 C, 6 N.**

Post-patch divergent-pair count: **22**.

## Real bugs (B) — fixed in this audit

See `doc/UPSTREAM_BUGS.md` §2026-05-11 ScaLAPACK residual-divergence
audit for the per-bug write-ups. Headline:

1. `pzungql.f:292-293` and `pzunml2.f:394-395` — `PB_TOPGET` instead
   of `PB_TOPSET` to restore BLACS broadcast topology. Z half buggy;
   C half correct. Pre-patch the recipe routed around these via
   `prefer_source: PCUNGQL, PCUNML2`; the patches let those entries
   retire.
2. `pdsyevd.f:225` — `LQUERY = (LWORK.EQ.-1)` missing the
   `.OR. LIWORK.EQ.-1` branch. Workspace queries via LIWORK=-1 alone
   incorrectly error. D half buggy; S half correct. Migrated qsyevd
   inherited the bug because D is the canonical default — this patch
   is the only one in the audit that changes migrated output.
3. `pslaed3.f:156` — `IINFO=0` (local) where it should be `INFO=0`
   (output parameter); leaves INFO undefined on successful return.
   And lines 168-171 lack the `IF(I+J.LE.N)` bounds guard on the
   INDROW/INDCOL writes — out-of-bounds write on the final
   outer-loop iteration when N isn't a multiple of NB. S half buggy;
   D half correct.

## B? deferred (1)

- **`pslaqr3` / `pdlaqr3`** — `LWK8 = 2*TZROWS*TZCOLS` on S vs
  `LWK8 = 0` on D. LWK8 feeds the LWKOPT max across LWK1..LWK8.
  Either D under-reports workspace (real bug) or S over-reports
  (cosmetic). Plus an extraneous `MPI_WTIME` in D's EXTERNAL list
  (likely dead). Needs QR-iteration domain inspection to decide
  which is canonical.

## W comparer-fixable (4)

Local-variable names embedding the precision letter that the
existing `_canonicalize_for_compare` regex doesn't reach:

- `bslaexc/bdlaexc` — `PARAMETER(TEN=10)` vs `PARAMETER(TEN=1.0E1)`.
  Integer literal vs typed-real literal of the same value.
- `psgesvd/pdgesvd` — locals `WPSLANGE`, `WPSGEBRD`, `WSBDSQR`,
  `WPSORMBR*` (S) vs `WPDLANGE`, etc. (D). The leading `W` blocks
  the existing `[SDCZ]` → `@` regex.
- `pzgesvd/pcgesvd` — same `WPZ*` / `WPC*` pattern.
- `pzheev/pcheev` — same `SIZEPZ*` / `SIZEPC*` pattern.

Deferred — comparer normalization at this depth (folding the
embedded letter inside intra-routine local names) is a bigger
regex-engineering project; not pursued in this audit.

## C cosmetic (11)

Pure upstream stylistic drift, no semantic difference:

- `pclarz/pzlarz`, `pclarzc/pzlarzc`, `pzlawil/pclawil` —
  `TAULOC(1)` (array element by reference) vs `TAULOC` (array by
  reference). Identical under Fortran pass-by-reference.
- `pslacon/pdlacon` — three separate `REAL` declarations vs one
  combined `REAL ESTWORK(1),TEMP(1),WORK(2)` line.
- `pslawil/pdlawil` — declaration block split / migrator emit
  artifact (`REAL BUF(4),1,1,1,1`).
- `pssyttrd/pdsyttrd`, `pzhettrd/pchettrd` — local rename
  `TTOPH`/`TTOPV` vs `CONJTOPH`/`CONJTOPV`. For the real
  symmetric / complex Hermitian split, conjugate-vs-identity is
  semantically equivalent.
- `pzlattrs/pclattrs` — single `REAL XMAX(1)` declaration moved.
- `pstrsen/pdtrsen` — multi-line `IF/CALL/ENDIF` vs single-line
  `IF(...) CALL`.
- `pzpotf2/pcpotf2` — `XXDOTC` (subroutine form, output arg) vs
  `XDOTC` (function form). Different ScaLAPACK BLAS API call shape,
  equivalent result.
- `psstebz/pdstebz` — D-half refactors `PSLAIECT` into split
  `PQLAIECTB`/`PQLAIECTL` with extra `IEFLAG` branch. Different
  algorithmic dispatch; numerically equivalent.

## N algorithmic tuning (6 — the MRRR family, S/D only)

Precision-tuned numerical constants in the MRRR (Multiple Relatively
Robust Representations) eigensolver:

| File | Tuning |
|---|---|
| `slarre2/dlarre2` | `PERT=4`/`8`; `RTL=HNDRD*EPS` vs `SQRT(EPS)` |
| `slarre2a/dlarre2a` | same as above |
| `slarrf2/dlarrf2` | `*TWO*EPS` vs `*FOUR*EPS` (sigma perturbation) |
| `sstegr2/dstegr2` | `MINRGP=3.0E-3` vs `1.0E-3` (relative gap) |
| `sstegr2a/dstegr2a` | same MINRGP |
| `sstegr2b/dstegr2b` | same MINRGP |

These mirror the same upstream tuning split documented in
`doc/lapack-residual-divergence-categorization.md` for
LAPACK's `slarre/dlarre/slarrf/dlarrf` family. Accepted as upstream
design.

## See also

- `doc/UPSTREAM_BUGS.md` §2026-05-11 ScaLAPACK residual-divergence audit
- `doc/lapack-residual-divergence-categorization.md` — same audit
  shape applied to LAPACK
- `doc/mumps-residual-divergence-categorization.md` — same for MUMPS
