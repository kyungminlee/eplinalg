# Upstream bugs in vendored Netlib sources

*Last catalogued: 2026-05-11*

This document catalogues bugs found in the vendored upstream sources
(`external/lapack-3.12.1/`, `external/scalapack-2.2.3/`,
`external/MUMPS_5.8.2/`) that the migrator works around without
editing `external/`. Each entry records the symptom, root cause, and
the in-tree workaround.

Per-library catalogues:

- [`lapack.md`](lapack.md) — LAPACK 3.12.1
- [`scalapack.md`](scalapack.md) — ScaLAPACK 2.2.3 (and PBLAS)
- [`mumps.md`](mumps.md) — MUMPS 5.8.2
- [`lapack-nits.md`](lapack-nits.md) — LAPACK 3.12.1 cosmetic dead declarations (catalogued, not patched)

The shared sections below — audit methodology, the consolidated bug
summary table, the "how fixes are carried" notes, and the
cross-library XERBLA-string sweep that touches both LAPACK and
ScaLAPACK files — remain in this index.

## Audit methodology

Bugs surfaced via five distinct audit passes:

1. **Convergence-diff audit** — for every co-family pair (S/D and
   C/Z), light-normalize both halves and report textual divergences.
   Catches single-precision-letter typos and one-half-only fixes.
   Yield: 11 bugs across BLAS/LAPACK/ScaLAPACK/MUMPS/PBLAS.
2. **Differential-precision tests** — runtime kind16-vs-Netlib-quad
   reference comparison. Catches numerical-correctness bugs that
   convergence misses (both halves wrong same way). Yield: pre-existing
   ScaLAPACK bugs (orbdb3 stride, lanhs underestimate, geequ axis,
   posvx LWMIN, syevx LQUERY, larz family, ungql/unml2 topology, etc.)
   plus the MUMPS bad-input SIGSEGV.
3. **XERBLA-string mechanical sweep** — verify each `CALL XERBLA('NAME', ...)`
   matches the enclosing SUBROUTINE name. 1 888 sites, 7 distinct typos
   (`'xORBDB'` placeholder unsubstituted in 4 halves; `bdtrexc`/`bstrexc`
   `'DTREXC'` instead of `'B?TREXC'`; `?pttrsv` reports `'?PTTRS'` not
   `'?PTTRSV'`).
4. **EXTERNAL-vs-CALL mechanical sweep** — verify every CALL target
   appears in EXTERNAL and vice versa. 462 routines flagged; 18 real
   bugs (8 patched, 10 S/C-only documented).
5. **PCHK?MAT param-position mechanical sweep** (ScaLAPACK only) —
   verify literal `MAPOS0`/`NAPOS0`/`DESCAPOS0` arguments to
   PCHK1MAT/PCHK2MAT match actual position of M/N/DESC in caller's
   signature. 316 sites, 53 mismatches in 11 routines (4 distinct
   bug families). **17 % bug rate — highest yield audit run.**
6. **INFO=-N param-position mechanical sweep** (LAPACK SRC) — for
   each `IF (test) THEN; INFO = -N` arg-validation block, identify
   the single signature arg referenced in the test and verify N
   matches its position. 11 352 sites, 9 distinct bugs found
   (single-arg LWORK / LRWORK validations whose INFO points to a
   different argument). Uncovered the striking `?ggsvd3` case where
   `INFO = -24` references INFO itself (an OUT-only argument the
   user can never pass illegally).
7. **Doc-header `\brief \b NAME` consistency** (LAPACK) — verify
   each `\brief \b NAME` matches the file's first SUBROUTINE/FUNCTION
   declaration. 1 847 files, 1 real bug (`zrscl.f` documents itself
   as `ZDRSCL`, the helper it calls).
8. **XERBLA-string trailing-space asymmetry** — verify halves agree
   on whether the routine name in `XERBLA('NAME', ...)` carries the
   trailing-space pad. 1 confirming hit (already-known `SSYSV_AA`).
9. **PXERBLA pattern audit** — categorize all 522 PXERBLA call sites
   in ScaLAPACK by third-argument shape (`-INFO`, `-IINFO`, descriptive
   suffixes, literal integer). Re-run SRNAME-vs-enclosing check with a
   PXERBLA-aware regex (the original sweep missed PXERBLA's
   `(ICTXT, 'NAME', ...)` form). 1 new bug: `psgebal.f:386` reports
   `'PDGEBAL'` (S-half typo). Verified `?laconsb` family's literal-
   integer calls (`PXERBLA(..., 10)`) correctly index LWORK at
   parameter position 10.
10. **LWORK-validation-gap mechanical sweep** — find routines with
    `LWORK`/`LIWORK`/`LRWORK` in signature that lack a matching
    `.LT.` validation. Most hits are LAPACK auxiliary routines (whose
    callers validate). 1 real bug in public-facing code:
    `pzheevd`/`pcheevd` miss LIWORK validation while all 14 sibling
    eigensolvers in the heev*/syev* family validate it correctly
    (complex-half asymmetry).
11. **Doc-header `\brief \b NAME` consistency** (already covered in #7
    above). 1 bug.
12. **Build-system registration drift** — verify CMakeLists.txt
    listings match the actual file set in SRC. LAPACK 2 055 SRC files
    + 56 in DEPRECATED, all registered correctly. ScaLAPACK 678/678
    perfect match. **No drift, no bugs.**

## Bug summary

| ID | Bug | Severity | Patched? | Audit |
|----|-----|---------:|:---:|-------|
| L01 | `?orbdb3.f` `?ROT` stride bug (LDX11 vs LDX21) | Memory | ✓ all 4 halves | runtime |
| L02 | `dgejsv.f` LWORK overflow in inner DGESVJ | Memory | ✓ D | conv |
| L03 | `?trsyl3.f` D/C/Z miss LIWORK / LDSWORK validation | Diag | ✓ D, Z | conv |
| L04 | Stale/missing EXTERNAL (5+ routines, 1st batch) | Adv | ✓ 4 D/Z | conv |
| L05 | `dlaswlq.f` accepts NB=0 | Validation | ✓ D | conv |
| L06 | `sggev3.f` query gates on ILVL not ILV | Workspace | — (S) | conv |
| L07 | `slarf1l.f` missing LASTV.GT.0 guards | Memory (OOB) | — (S) | conv |
| L08 | `slarrf.f` 2×eps shift safety vs dlarrf 4×eps | Numerical | — (S) | conv |
| L09 | `slasq2.f` IEEE hard-coded false | Performance | — (S) | conv |
| L10 | `sgelss.f` ILAENV-estimated workspace | Workspace | — (S) | conv |
| L11 | `sgejsv.f` JOBA='L' should be 'G' | Numerical | — (S) | conv |
| L12 | `zla_syrfsx_extended.f` reports 'ZLA_HERFSX_EXTENDED' | Diag | ✓ Z | conv |
| **X01** | `?orbdb.f`/`?unbdb.f` `'xORBDB'` placeholder (all 4 halves) | Diag | ✓ D, Z | xerbla |
| **X02** | `bdtrexc.f`/`bstrexc.f` line 173 `'DTREXC'` typo | Diag | ✓ D | xerbla |
| **X03** | `?pttrsv.f` reports `'?PTTRS'` not `'?PTTRSV'` (all 4) | Diag | ✓ D, Z | xerbla |
| **E01** | 8 D/Z EXTERNAL drifts (dlarf1l, zlarf1l, dlatrs3, zlatrs3, zgelss, dorbdb4, dorgr2, zrscl) | Adv | ✓ 8 | extern |
| **E02** | 10 S/C-half EXTERNAL drifts (clarf1l, slarf1l, slatrs3, clatrs3, cgelss, sorbdb4, sorgr2, crscl, dopmtr, clarf1f) | Adv | — | extern |
| **P01** | `?heevr/?syevr` DESCZ POS0 21 vs 19 (all 4 halves) | Diag | ✓ D, Z | pchk |
| **P02** | `?hegvx/?sygvx` N POS0 4 vs 5 (all 4 halves) | Diag | ✓ D, Z | pchk |
| **P03** | `pdtrord/pstrord` all POS0 off by +1 | Diag | ✓ D | pchk |
| **P04** | `pzheevd.f` DESCZ POS0 11 vs 12 | Diag | ✓ Z | pchk |
| **I01** | `?orbdb4`/`?unbdb4` INFO=-14 (WORK) for LWORK test, should be -15 | Diag | ✓ D, Z | info-n |
| **I02** | `?ggsvd3` INFO=-24 (INFO itself!) for LWORK test, should be -22 | Diag | ✓ D, Z | info-n |
| **I03** | `?orcsd`/`?uncsd` INFO=-22 (LDU2) for LWORK test, should be -28 | Diag | ✓ D, Z | info-n |
| **I04** | `zuncsd` INFO=-24 (LDV1T) for LRWORK test, should be -30 | Diag | ✓ Z | info-n |
| **I05** | `?laqz0` INFO=-19 (RWORK) for LWORK test, should be -18 | Diag | ✓ Z | info-n |
| **I06** | `?laqz2` INFO=-26 (RWORK) for LWORK test, should be -25 | Diag | ✓ Z | info-n |
| **D01** | `zrscl.f` doc header `\brief \b ZDRSCL` should be ZRSCL | Doc | ✓ Z | brief |
| **X04** | `psgebal.f:386` PXERBLA reports `'PDGEBAL'` (line 225 correct) | Diag | — (S) | pxerbla |
| **W01** | `pzheevd`/`pcheevd` miss LIWORK validation (memory bug; complex-half-only asymmetry vs heevd/syevd siblings) | Memory | ✓ Z, — C | lwork-gap |
| **I07** | `?pttrs` / `?pbtrsv` LWORK sentinel-check INFO off by 1 (8 routines, all 4 halves) | Diag | ✓ D, Z (4 routines) | info-n |
| **I08** | `?lascl` CFROM/CTO INFO codes -4/-5 inherited from LAPACK signature; should be -2/-3 (4 halves) | Diag | ✓ D, Z | info-n |
| **I09** | `?ormrz` / `?unmrz` duplicate K-validation (L never checked; 4 halves) | Validation | ✓ D, Z | info-n |
| S01 | `pzunmbr.f` EXTERNAL declares PCHK1MAT, body calls PCHK2MAT | Adv | ✓ Z | conv |
| S02 | `pssyevd.f` LQUERY misses LIWORK=-1 | Validation | — (S) | conv |
| S03 | `pslaed3.f` clobbers user INFO | Validation | — (S) | conv |
| S04 | `p{d,z}atrmv_.c` ALPHA hardcoded one (UPLO=L, TRANS=T/C) | Numerical | ✓ D, Z | runtime |
| S05 | `p?lanhs.f` NPROW=1 underestimate | Numerical | ✓ D, Z | runtime |
| S06 | `p?lanhs.f` IAROW double-advance | Numerical | ✓ D, Z | runtime |
| S07 | `p?geequ.f` column-scale wrong axis | Numerical | ✓ D, Z | runtime |
| S08 | `p?posvx.f` LWMIN too small | Validation | ✓ D, Z | runtime |
| S09 | `pdsyevx.f`/`pzheevx.f` LQUERY-path early-write | Memory | ✓ D, Z | runtime |
| S10 | `pdtrsen.f` IWORK early-write during LQUERY | Memory | ✓ D | runtime |
| S11 | `pdlarzb.f` PBDTRAN N-arg vs LV-buffer mismatch | Numerical | ✓ D | runtime |
| S12 | `p?ormrz.f` post-loop condition copy-paste | Numerical | ✓ all | runtime |
| S13 | `p?larz.f` MPV/NQV undersizing (heap overrun) | Memory | ✓ all | runtime |
| S14 | `p?larz.f` ZAXPY stride bug | Numerical | ✓ all | runtime |
| S15 | `pzlarz.f`/`pzlarzc.f` SIDE='L' missing ZLACGV; ZGERC vs ZGERU | Numerical | ✓ Z | runtime |
| S16 | `p?larz.f`/`p?larzc.f` PBxTRNV reads M where L stored | Memory | ✓ all | runtime |
| S17 | `pzungql.f`/`pzunml2.f` PB_TOPGET should be PB_TOPSET | Restore | ✓ Z | runtime |
| M01 | MUMPS 5.8.2 SIGSEGV on bad-input cases | Robustness | — | runtime |

**Severity legend:**
- *Memory*: out-of-bounds read/write, possible corruption.
- *Numerical*: wrong result returned for valid input.
- *Validation*: invalid input not rejected (silent bad behavior).
- *Workspace*: under-allocation in size queries; overflow possible.
- *Diag*: only the diagnostic message is wrong; rejection still works.
- *Adv*: advisory only (Fortran linker resolves regardless).
- *Performance*: correct result but slower than necessary.
- *Restore*: incorrect post-call cleanup of communicator state.
- *Robustness*: graceful error reporting required by spec, but routine SIGSEGVs.

## How fixes are carried

Recipes carry a ``patches:`` list (see ``recipes/README.md``) that
applies declarative diffs to a staged copy of the upstream sources
before migration. Patches are written in upstream shape
(``DOUBLE PRECISION`` types, ``pd*``/``pz*`` symbol names,
``dgemm`` call sites, …) so a single patch produces correctly
renamed/promoted output for every target. The standard-precision
archive built from the unmodified ``external/`` tree is unaffected
— only the migrated extended-precision archive carries the fix.

When the convergence picker would otherwise pick the un-fixed C/S
half over the patched D/Z half, the recipe's ``prefer_source`` field
pins the correct canonical (the rank picker doesn't recognize
ScaLAPACK's ``pd*``/``pz*`` two-letter prefix — first character is
always ``P`` — so it sorts alphabetically by file name).

---


## LAPACK 3.12.1 + ScaLAPACK 2.2.3: XERBLA routine-name string typos (sweep 2026-05-08)

A mechanical XERBLA-string sweep across all four upstream Fortran
trees (LAPACK SRC + INSTALL, ScaLAPACK SRC + TOOLS + REDIST + PBLAS,
BLAS, MUMPS) found 21 `CALL XERBLA('NAME', -INFO)` / `PXERBLA(...)`
sites where the reported routine name doesn't match the enclosing
SUBROUTINE. After triaging out intentional forwarding patterns
(`?sytrd_2stage` / `?hetrd_2stage` reporting their inner SY2SB / SB2ST
sub-routines on inner-call errors; `?la_heamv` reporting the related
`?HEMV` for documentation purposes), seven sites are real typos:

**Patched in migrated archive (D/Z-canonical halves):**

| File | Line | Bug | Fix |
|------|-----:|-----|-----|
| `dorbdb.f` | 385 | `'xORBDB'` (lowercase x — script-generation placeholder, never substituted) | `'DORBDB'` |
| `zunbdb.f` | 388 | `'xORBDB'` | `'ZUNBDB'` |
| `bdtrexc.f` (ScaLAPACK) | 173 | `'DTREXC'` — but eight *other* XERBLA calls in the same routine (lines 242–342) correctly report `'BDTREXC'`, so this one is internally inconsistent | `'BDTREXC'` |
| `dpttrsv.f` (ScaLAPACK) | 100 | `'DPTTRS'` — routine name is `DPTTRSV`, doesn't call DPTTRS internally | `'DPTTRSV'` |
| `zpttrsv.f` (ScaLAPACK) | 115 | `'ZPTTRS'` | `'ZPTTRSV'` |

Carried in `recipes/lapack/patches/{dorbdb,zunbdb}.f` and
`recipes/scalapack/patches/{bdtrexc,dpttrsv,zpttrsv}.f`,
wired in the matching recipe files.

**Documented but not patched (S/C-half non-canonical, never reaches
migrated archive):**

| File | Line | Bug |
|------|-----:|-----|
| `sorbdb.f` | 383 | `'xORBDB'` — should be `'SORBDB'` |
| `cunbdb.f` | 386 | `'xORBDB'` — should be `'CUNBDB'` |
| `bstrexc.f` (ScaLAPACK) | 173 | `'DTREXC'` — should be `'BSTREXC'`. Worse than the bdtrexc case because the precision letter is also wrong (D where the file is S). |
| `spttrsv.f` (ScaLAPACK) | 100 | `'SPTTRS'` — should be `'SPTTRSV'` |
| `cpttrsv.f` (ScaLAPACK) | 115 | `'CPTTRS'` — should be `'CPTTRSV'` |

**Severity.** All diagnostic-only — XERBLA prints the supplied string
and aborts; the routine still rejects the bad input. Users who
encounter `XERBLA: 'xORBDB' Parameter -7 had an illegal value` will
be confused (no LAPACK routine is called `xORBDB`); same with
`'BSTREXC' → 'DTREXC'` (looking up DTREXC docs gets you the wrong
arg-list reference). No memory-safety or numerical impact.

**Why upstream's tests miss them.** Reference test drivers always
pass valid arguments, so the XERBLA-error path never fires in CI.
The strings are only seen by users who pass invalid arguments and
read the diagnostic output.

**Already documented separately:** `zla_syrfsx_extended.f:496` reports
`'ZLA_HERFSX_EXTENDED'` (Hermitian) inside the symmetric routine —
caught earlier and patched in `recipes/lapack/patches/`.

**Re-sweep with PXERBLA-aware regex (2026-05-09):** The earlier sweep
matched `XERBLA('NAME', ...)` form only, missing ScaLAPACK's
`PXERBLA(ICTXT, 'NAME', INFO)` form (where the name is the *second*
argument). Re-running with `\bXERBLA\s*\(\s*'…'` and
`\bPXERBLA\s*\(\s*\w+\s*,\s*'…'` regexes surfaced one additional
mismatch:

- `psgebal.f:386` reports `'PDGEBAL'` (double-precision) inside the
  single-precision routine `PSGEBAL`. Same flavor as `bstrexc.f`
  reporting `'DTREXC'`. The other PXERBLA call in psgebal at line
  225 correctly reports `'PSGEBAL'`, so this is a copy-paste typo
  from `pdgebal.f` (which has `'PDGEBAL'` correctly at both lines).
  S-half-only; not patched (non-canonical). Documented for upstream.

**Upstream report.** Not yet filed.

---
