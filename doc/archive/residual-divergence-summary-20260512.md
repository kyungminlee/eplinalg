# Residual divergence summary across all libraries (2026-05-12)

Cross-library tally after the 2026-05-11 LAPACK / MUMPS / ScaLAPACK
residual-divergence audits and the 2026-05-12 comparer-wiring session
(`_filter_precision_drift` wired into the divergence-report path;
`_sort_decl` switched to a precision-neutral key;
numeric-literal scientific → decimal collapse for integer-valued
exponents).

Reproducer (per-library counts):

```
for lib in blas blacs lapack pblas scalapack mumps; do
    echo -n "$lib: "
    uv run python -m migrator diverge recipes/$lib.yaml \
        --target kind16 --no-whitelist --context 0 \
        2>&1 | grep -c '^###'
done
```

## Per-category breakdown

| Library | C | D | N | R | W | B? | B | **Total** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| blas | 2 | 0 | 1 | 0 | 0 | 0 | 0 | **3** |
| blacs | 0 | 0 | 0 | 0 | 0 | 0 | 0 | **0** |
| lapack | 60 | 22 | 9 | 9 | 5 | 1 | 0 | **106** |
| pblas | 0 | 0 | 0 | 0 | 0 | 0 | 0 | **0** |
| scalapack | 11 | 0 | 6 | 0 | 0 | 1 | 0 | **18** |
| mumps | 4 | 0 | 6 | 0 | 0 | 0 | 0 | **10** |
| **Total** | **77** | **22** | **22** | **9** | **5** | **2** | **0** | **137** |

## Categories

- **C — cosmetic upstream drift.** Label numbering, comment prose,
  declaration ordering, string trailing whitespace. Both halves are
  numerically and behaviorally identical; the diff is purely textual.
  **No action.** Examples: `scasum/dzasum` prose comment block,
  `sdot/ddot` DO/CONTINUE label numbers, MUMPS `?ana_driver` banner
  text.
- **D — dead declarations.** Variables / PARAMETER constants /
  EXTERNAL names declared on one half but never referenced. Each
  contributes ~1-3 lines to the per-pair diff. Catalogued in
  `doc/UPSTREAM_NITS-LAPACK.md`. **No action** — would require 22 one-line
  patches with zero functional effect.
- **N — algorithmic tuning constants.** Numeric constants
  deliberately tuned per precision. The MRRR family (`larre`,
  `larrf`, `lasq2`, `stegr2*`) uses different `MINRGP`, `PERT`,
  `RTL`, `FOUR*EPS` vs `TWO*EPS`. BLAS `srotmg` has
  single-precision-specific underflow constants. MUMPS `KEEP(122)`
  / `KEEP(421)` have per-precision defaults. **No action — upstream
  design.**
- **R — refactor / structurally equivalent.** Same computation
  expressed differently across halves: MUMPS MPI_REDUCE temp-var
  staging (S/C buffer through `TMPTIME` + master-guard, D/Z reduce
  directly), ScaLAPACK `pzpotf2` function-vs-subroutine BLAS form,
  declaration block splits. Both produce identical results. **No
  action.**
- **W — comparer-fixable, deferred.** Drift the comparer
  *could* fold but doesn't yet. The 5 remaining lapack W's need
  wider canonicalization than the 2026-05-12 wiring covered. MUMPS
  W's (two pairs from the original audit) and ScaLAPACK W's (four
  pairs) were folded by the comparer changes this session.
- **B? — ambiguous bug, needs domain review.** Real source-level
  divergence whose correctness cannot be determined by static
  reading alone. The 3 remaining (lapack `cgedmd`/`zgedmd`,
  `sgejsv`/`dgejsv`; scalapack `pslaqr3` LWK8) live in newer or
  domain-specific routines (GEDMD/GEJSV, QR-iteration workspace
  tuning). Resolution path: file as **issues** at Netlib LAPACK /
  Netlib ScaLAPACK trackers (not PRs — we don't know which half is
  right). See `doc/UPSTREAM_BUGS.md` §B? reporting note.
- **B — real bug.** Confirmed upstream bug with a known correct
  fix. **Zero remaining**: every B identified by the audits is
  patched (`recipes/<lib>/patches/`) and documented at the
  ##-level in `doc/UPSTREAM_BUGS.md`.

## Audit-actionable surface

| Class | Count | Status |
|---|---:|---|
| B | 0 | All patched |
| B? | 3 | Deferred — file as Netlib issues |
| W | 5 | Deferred — wider comparer work |
| **Audit-actionable** | **8 / 138** | **5.8%** |
| C + D + N + R | 130 | Non-actionable upstream design |

The audit is at long-tail. Future per-pair diffing of the
remaining 130 non-actionable pairs is unlikely to yield bugs; new
discoveries from here will require either deeper algorithmic
context (the 3 deferred B?) or different audit modalities
(Valgrind/ASan on the test suite) rather than static source
comparison.

## What the 2026-05-12 comparer work changed

| Change | Effect |
|---|---|
| Wire `_filter_precision_drift` into `run_divergence_report` (it was defined but never called) | Folds line-level precision-letter token drift uniformly |
| Switch `_sort_decl` to precision-neutral sort key (S→D, C→Z before comparing) | Lets co-family local-variable lists like `WPSLANGE,WPSGEBRD,WSBDSQR` align with `WPDLANGE,WPDGEBRD,WDBDSQR` after sorting so the line-level filter recognizes them |
| Collapse `<int>.<zeros>?E<+?int>` scientific notation to expanded integer when exp ≤ 18 | Folds `1.0E1` ↔ `10` form drift |

Pre-session counts (post 2026-05-11 audit, pre 2026-05-12 comparer work):

```
blas: 3, blacs: 0, lapack: 109, pblas: 0,
scalapack: 22, mumps: 10 — total 144
```

Post-session counts:

```
blas: 3, blacs: 0, lapack: 107, pblas: 0,
scalapack: 18, mumps: 10 — total 138
```

**Net: 6 pairs folded.** No source patches in this session — purely
comparer normalization.

## See also

- `doc/lapack-residual-divergence-categorization.md` — per-pair
  lapack classification (130 → 122 → ... → 107)
- `doc/scalapack-residual-divergence-categorization.md` — per-pair
  scalapack classification (26 → 22 → 18)
- `doc/mumps-residual-divergence-categorization.md` — per-pair
  mumps classification (162 → 14 → 10)
- `doc/UPSTREAM_BUGS.md` — bug catalogue with patches and Netlib
  upstream-report guidance
- `doc/UPSTREAM_NITS-LAPACK.md` — D-class dead declarations across all
  libraries
- `doc/CONVERGENCE.md` — overall convergence framework
