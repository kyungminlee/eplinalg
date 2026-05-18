# MUMPS residual-divergence categorization (2026-05-11)

After comparer canonicalization (`_canonicalize_for_compare` plus
precision-letter / kind-suffix masking), MUMPS's 162 raw S↔D / C↔Z
divergent pairs collapse to **10 residual pairs** with non-trivial
diff content.

The reduction from 14 → 10 reported pairs happened in this audit
when two W-class items were collapsed into the comparer:

1. The numeric-literal kind regex was extended from `(\d)[dD](...)` to
   `([\d.])[dD](...)` so `1.D0` (a typed-real literal with no leading
   digit) normalizes alongside `1D0`. This cleared the
   `sana_aux_par` ↔ `dana_aux_par` divergence.
2. A library-root prefix-letter fold was added: any `[SDCZ]`
   immediately preceding `MUMPS_` / `LAPACK` / `BLAS` / `BLACS` /
   `SCALAPACK` collapses to `@`, regardless of preceding-character
   context. This catches precision-letter routine names embedded
   inside error-string literals where the existing word-boundary
   regex deliberately did not reach. This cleared both
   `?mumps_comm_buffer` pairs.

See `src/migrator/pipeline.py` `_canonicalize_for_compare` steps 1
and 5_str.

Reproducer:

```
uv run python -m migrator diverge recipes/mumps.yaml --target kind16 --no-whitelist
```

Audit categories follow the same scheme used for LAPACK
(`doc/lapack-residual-divergence-categorization.md`):
**C** cosmetic / **N** algorithmic tuning / **W** comparer-fixable /
**B?** ambiguous / **B** real bug.

## Summary

| Pair (both halves) | Symptom | Class |
|---|---|---|
| `?ana_driver` | Storage-banner string differs (`'% SINGLE PRECISION STORAGE'` vs `'% REAL STORAGE'`; `'% SINGLE COMPLEX STORAGE'` vs `'% COMPLEX STORAGE'`) | C |
| `?fac_distrib_distentry` | Trailing-space difference in error-message string literal | C |
| `?fac_driver` | MPI_REDUCE intermediate-variable pattern: S/C use `TMPTIME` / `TMPFLOP` scratch + `IF(MYID.EQ.MASTER)` guard + `DBLE()` cast on `RINFOG`; D/Z reduce directly into the final variable without guard | N |
| `?ini_defaults` | Algorithmic tuning constants: `KEEP(122)=15`/`150` and `KEEP(421)=1000`/`500` (C/Z), `=3000`/`650` (S/D) | N |
| `?tools` | C-API split: S/C route through `MUMPS_SETRVAL_ADDR_C` with a `REAL` temporary; D/Z call `MUMPS_SETDVAL_ADDR_C` directly with the `DOUBLE PRECISION` value | N |

10 pairs total: 4 C, 6 N. **Zero B / B? findings — no upstream bug,
no migrator bug.** The two W pairs from the initial audit are now
folded into the comparer.

## Per-pair detail

### C: cosmetic upstream string drift

- **`?ana_driver.F`** — banner-text difference in the analysis
  driver's `%FILE` dump. Each half writes its own marketing label
  (`'SINGLE PRECISION STORAGE'` vs `'REAL STORAGE'`, etc.). String
  content only, no numerical or interface impact.
- **`?fac_distrib_distentry.F`** — `'**ERROR ALLOCATING REAL BUFFER
  FOR MATRIX DISTRIBUTION '` has a trailing space on one half but
  not the other. Pure typo-level upstream drift.

No action; expected to remain on the divergence report indefinitely.

### W (resolved in comparer): precision-letter / kind-suffix masking

Two W-class items were folded into `_canonicalize_for_compare`
(see "Summary" above for the regex changes):

- **`?mumps_comm_buffer.F`** — error string `'... IN SMUMPS_BUF_…'` vs
  `'... IN DMUMPS_BUF_…'`. Resolved by the `[SDCZ](?=MUMPS_|LAPACK|…)`
  fold at step 5_str.
- **`?ana_aux_par.F`** — `SYMMETRY = 1.D0` (S side) vs `= 1` (D side);
  `* 100.D0` vs `* 100`. Resolved by extending the `(\d)[dD]`
  exponent regex to `([\d.])[dD]` at step 1, which then chains into
  the existing `1.0` → `1` fold.

Both are now invisible on the divergence report.

### N: algorithmic tuning / design asymmetry

- **`?fac_driver.F`** — S/C halves stage MPI_REDUCE results through
  `TMPTIME` / `TMPFLOP` scratch variables and gate the assignment on
  `MYID.EQ.MASTER`; D/Z reduce directly into `TIMEMAX` / `FLOPMAX`.
  Plus S/C use `DBLE(ID%RINFOG(3))` to widen the input to
  `MUMPS_COMPUTE_GLOBAL_GAINS`; D/Z pass `ID%RINFOG(3)` directly.
  Both forms produce identical output on the master rank (the only
  rank that reads the values); the IF-guard variant just leaves the
  prior value untouched on non-master ranks, the direct-write
  variant overwrites with `MPI_REDUCE`'s undefined-on-non-root
  sentinel. Upstream is hand-maintained S↔D asymmetric, not a bug.

- **`?ini_defaults.F`** — default values for `KEEP(122)` (block
  granularity) and `KEEP(421)` (analysis-blocking auto-tune
  threshold) are deliberately precision-tuned: complex precisions
  get smaller defaults than real, single smaller than double. Pure
  algorithmic-tuning constants; expected to differ.

- **`?tools.F`** — `MUMPS_UPDATE_PROGRESS` routes its `OPELI`
  (DOUBLE PRECISION) flop counter into the C-side progress display
  via a precision-tagged C entry point: S/C cast down to `REAL` and
  call `MUMPS_SETRVAL_ADDR_C`; D/Z call `MUMPS_SETDVAL_ADDR_C`
  directly. The C runtime exposes both entries; the upstream split
  is deliberate ABI design (single and double progress storage are
  registered separately on the C side). Not unifiable without
  collapsing the C entry points.

## Aggregate conclusion

- The 162 raw MUMPS divergences reduce to 10 substantive residual
  pairs after comparer canonicalization (including the two W-class
  folds added by this audit).
- None are migrator bugs; none are upstream bugs.
- 4 are pure cosmetic upstream string drift, accepted as-is.
- 6 are intentional upstream algorithmic / ABI asymmetry between
  S/C and D/Z halves: MPI_REDUCE staging, tuning constants, and
  the C-progress-display ABI split.

Net migrator impact: **zero source patches needed**; two comparer
normalization rules added. The existing `defer_all_divergences: true`
in `recipes/mumps.yaml` continues to be the correct stance for MUMPS.

## See also

- `doc/archive/convergence-20260507.md` §mumps — original 162/200 figure and
  why upstream maintains S/D as hand-edited copies.
- `doc/lapack-residual-divergence-categorization.md` — same audit
  applied to LAPACK's 122-pair residual.
- `recipes/mumps.yaml` — `defer_all_divergences: true` rationale.
