# Convergence snapshot — 2026-05-09

Per-library convergence of the kind16 staged tree at `/tmp/stage-q/`
against fresh in-memory re-migration of each S/C sibling. Generated
with:

```bash
uv run python -m migrator stage /tmp/stage-q --target kind16
for lib in blas xblas blacs lapack pbblas pblas ptzblas \
           scalapack scalapack_c scalapack_tools mumps; do
    uv run python -m migrator converge \
        recipes/$lib.yaml /tmp/stage-q/$lib/src --target kind16
done
```

Numbers below are RAW (with each recipe's `expected_divergences:` /
`defer_all_divergences:` whitelist temporarily disabled). After the
whitelist is reapplied, every library returns rc=0 from the
convergence gate.

See `doc/CONVERGENCE.md` for methodology and
`doc/refactor-20260509.md` for the pipeline reshape that produced
this snapshot.

## Summary

| library          |   pairs | converged | diverged | missing | converge % |
|------------------|--------:|----------:|---------:|--------:|-----------:|
| blas             |      75 |        68 |        7 |       0 |     90.7 % |
| xblas            |       0 |         0 |        0 |       0 |        n/a |
| blacs            |      43 |        42 |        1 |       0 |     97.7 % |
| lapack           |   1 018 |       857 |      161 |       0 |     84.2 % |
| pbblas           |      14 |        14 |        0 |       0 |    100.0 % |
| pblas            |      61 |         0 |       61 |       0 |      0.0 % |
| ptzblas          |      47 |        47 |        0 |       0 |    100.0 % |
| scalapack        |     336 |       297 |       39 |       0 |     88.4 % |
| scalapack_c      |       3 |         2 |        1 |       0 |     66.7 % |
| scalapack_tools  |       0 |         0 |        0 |       0 |        n/a |
| mumps            |     200 |        38 |      162 |       0 |     19.0 % |
| **total**        | **1 797** | **1 365** |  **432** |    **0** |  **76.0 %** |

Zero `missing` across the board.

## Delta vs the pre-refactor snapshot (2026-05-07)

| library     | 2026-05-07 | 2026-05-09 | delta | driver                                                                                                       |
|-------------|-----------:|-----------:|------:|--------------------------------------------------------------------------------------------------------------|
| blas        |          7 |          7 |     — | unchanged — all upstream-design drift (sdot line-swap, srotmg constants, scasum prose, B01 D/Z override)     |
| blacs       |          1 |          1 |     — | unchanged                                                                                                    |
| lapack      |        430 |    **161** |  **−269** | patches now apply at the migrator-input layer (Phases B+C) so cosmetic D/Z↔S/C asymmetry collapsed across many headers, INTRINSIC lists, dead PARAMETER blocks |
| pblas       |         61 |         61 |     — | migrator-side K&R re-emergence — needs C-pass fix                                                            |
| scalapack   |         35 |     **39** |  **+4** | net of: 13 new D/Z-only patches that introduced text drift vs unpatched S/C minus 9 PS/PC sibling patches added in commit `ad97572f` |
| scalapack_c |          1 |          1 |     — | migrator-side `#define TYPE complex` rename gap                                                              |
| mumps       |        162 |        162 |     — | upstream maintains S/D as separately hand-edited copies                                                      |
| **total**   |    **697** |    **432** | **−265** |                                                                                                              |

The two-thirds shrink in LAPACK is the headline result of the
refactor: by routing migration through `build/staged-sources/<lib>/`
(post-patch) instead of raw `external/`, every patch's effect is
seen by both halves at re-migration time. Many divergences in the
2026-05-07 cohort were actually "we patched D, S still produces a
slightly different text." After Phase B+C those collapse without
needing a sibling patch.

The +4 in ScaLAPACK is the symptom that motivated Phase D: D/Z-only
patches that fix real bugs DO introduce text drift between the
patched D/Z half and the unpatched S/C half. As PS/PC sibling
patches are added (the auto-translator landed 9 in commit
`ad97572f`), the count will shrink. 19 PD/PZ patches in
`recipes/scalapack.yaml`'s `asymmetric_patches:` still need
hand-translated PS/PC siblings.

## Per-library findings

### blas — 7 / 75 diverged

Audit unchanged from 2026-05-07. All 7 are upstream cosmetic /
structural drift; none affect numerics. Logged in
`recipes/blas.yaml`'s `expected_divergences:` (DDOT / DNRM2 / DROTG
/ DROTMG / DZASUM / DZNRM2 / ZROTG).

The B01 patches (`recipes/blas/patches/B01__nrm2-rotg-lazy-init.patch`)
fix a multifloats build issue on D/Z; no S/C sibling exists because
the workaround is target-specific and S/C halves don't pick the
canonical for kind16/multifloats targets. Logged in `blas.yaml`'s
`asymmetric_patches:`.

### xblas — 0 / 0

No co-family pairs.

### blacs — 1 / 43 diverged

`sgsum2d_/dgsum2d_ → qgsum2d_`. The B01 patch installs a
user-defined MPI op (`MPI_Op_create(BI_dMPI_sum, ...)`) on the
D-half to bypass Intel MPI 2021.17's broken `MPI_SUM` on
`MPI_REAL16` / `MPI_LONG_DOUBLE`. The S half doesn't need it (no
MPI op corruption at fp32). Whitelisted as expected (`DGSUM2D_`)
and listed in `asymmetric_patches:`.

### lapack — 161 / 1 018 diverged

The big winner of the refactor: 430 → 161, a 63 % reduction.

The remaining 161 break down approximately:

- **Dead PARAMETER ONE block** (~38 pairs): `dgehd2`, `dgelq2`,
  `dgeql2`, `dgeqr2`, `dgeqr2p`, `dgerq2`, `dorm2l`, `dorm2r`,
  `dorml2` (D-side) and the matching Z-side cohort. The patches
  strip the dead block on D/Z so the migrated text matches the
  S/C half's already-stripped form. Recorded under
  `one_sided_cleanup:`. The remaining diff is text the patch
  doesn't quite normalize to S/C's exact form (line-spacing /
  comment differences).
- **Header / INTRINSIC / literal-style drift** (~120 pairs):
  routine prologues, `INTRINSIC` lists, comment paraphrasing
  that drifts between upstream halves. Not patched (would touch
  hundreds of files for cosmetic-only gain).
- **Real-bug patches whose S/C sibling still needs hand-translation**
  (~3 pairs): `dtrsyl3` / `ztrsyl3` are missing `strsyl3`
  sibling; `dorcsd` / `zuncsd` / `cuncsd` are missing `sorcsd`
  sibling.

All logged in `recipes/lapack.yaml`'s `expected_divergences:`
(161 stems enumerated).

### pbblas — 0 / 14

Clean.

### pblas — 61 / 61 diverged

Unchanged from 2026-05-07. The PBLAS C re-migration path leaves
K&R fallback blocks (`#ifdef __STDC__ / #else …`) intact while the
on-disk canonical strips them. Migrator-side bug; not addressed by
this refactor (out of scope per `doc/refactor-20260509.md`).
`recipes/pblas.yaml` has `defer_all_divergences: true`; the field
should be flipped back to enumerated `expected_divergences:` once
the C-pass K&R stripper is fixed.

### ptzblas — 0 / 47

Clean.

### scalapack — 39 / 336 diverged

Up 4 from 2026-05-07. Composition:

- **MINRGP tuning split** (≈12 pairs): `pdstegr2`, `pdstegr2a`,
  `pdstegr2b`, `pdlarre2`, `pdlarre2a`, `pdlarrf2`. Confirmed
  upstream-design (S uses `MINRGP = 3.0E-3`, D uses `1.0E-3`).
  Listed in `expected_divergences:`.
- **PD/PZ patches missing PS/PC siblings** (≈19 pairs): real
  upstream bugs in geequ / posvx / larzb / larz / lascl / pbtrsv
  / pttrs / heevd / unmbr / etc. The auto-translator (in commit
  `ad97572f`) couldn't generate the PS/PC versions cleanly;
  hand-translation is needed. Listed in `asymmetric_patches:`.
- **Other upstream asymmetries** (≈8 pairs): `pdstebz`, `pdlarz`,
  `pdlanhs`, `pdlaqr3`, `pdtrsen`, `pclawil`, etc. — algorithmic
  S/D drift the migrator cannot reconcile.

All 39 stems logged in `recipes/scalapack.yaml`'s
`expected_divergences:`.

### scalapack_c — 1 / 3 diverged

Unchanged. `clamov_/zlamov_ → xlamov_`. Migrator C-pass leaves
`# define TYPE complex` in the C-half re-migration where the
canonical has `QCOMPLEX`. `defer_all_divergences: true` — flip
back to enumerated when the migrator gap is fixed.

### scalapack_tools — 0 / 0

No pairs.

### mumps — 162 / 200 diverged

Unchanged. MUMPS upstream maintains S and D halves as separately
hand-edited copies; convergence is a property of upstream, not
the migrator. The 38 / 200 that converge are the small utility
modules where upstream happens to keep S and D in lock-step. The
guarantee that matters for MUMPS is differential-precision parity
at the test level (currently 26/26). `defer_all_divergences: true`.

## Audit status (2026-05-09)

| library         | pre-refactor | post-refactor | gate status | follow-on work                                              |
|-----------------|-------------:|--------------:|:-----------:|-------------------------------------------------------------|
| blas            |            7 |             7 |     ✓       | none — upstream-design                                      |
| blacs           |            1 |             1 |     ✓       | none — kind16-target-specific                               |
| lapack          |          430 |           161 |     ✓       | hand-translate `strsyl3`, `sorcsd` siblings (3 stems)      |
| pblas           |           61 |            61 |     ✓ defer | migrator C-pass K&R stripper fix                           |
| scalapack       |           35 |            39 |     ✓       | hand-translate 19 PS/PC sibling patches (`asymmetric_patches:`) |
| scalapack_c     |            1 |             1 |     ✓ defer | migrator C-pass `#define TYPE complex` rename                |
| mumps           |          162 |           162 |     ✓ defer | none (upstream property)                                    |

CI gate (`migrator converge`) returns rc=0 for all 9 libraries:
the kind16 staged tree at `/tmp/stage-q` passes the convergence
check.

End-to-end test status (`ctest --test-dir build -j8`, kind16,
linux-impi preset): **1151/1151 pass** in 209s — no regression
from the refactor (verified at commit `ad97572f`).

## Reproducing

```bash
uv run python -m migrator stage /tmp/stage-q --target kind16
for lib in blas blacs lapack pblas pbblas ptzblas \
           scalapack scalapack_c mumps; do
    uv run python -m migrator converge \
        recipes/$lib.yaml /tmp/stage-q/$lib/src --target kind16 \
        > /tmp/conv_$lib.out 2>&1
    echo "$lib: rc=$? $(tail -1 /tmp/conv_$lib.out)"
done
```

Build + test (Intel MPI canonical):

```bash
cd /tmp/stage-q
cmake --preset linux-impi
cmake --build build -j8
ctest --test-dir build -j8 --output-on-failure
```

## Cross-references

- `doc/CONVERGENCE.md` — methodology (`converge` vs `diverge`,
  status semantics, suppressed patterns).
- `doc/convergence-20260507.md` — pre-refactor snapshot.
- `doc/refactor-20260509.md` — pipeline reshape: prepare →
  patch → migrate → diverge → pick → build.
- `doc/UPSTREAM_BUGS.md` — bug catalogue with patch filename
  cross-references.
- `recipes/<lib>.yaml` — `expected_divergences:`,
  `defer_all_divergences:`, `asymmetric_patches:`,
  `one_sided_cleanup:` fields.
