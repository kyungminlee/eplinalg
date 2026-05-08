# Convergence

## What it is

Upstream BLAS / LAPACK / ScaLAPACK / MUMPS ship every routine in two
co-family halves: `S`/`D` for real, `C`/`Z` for complex. After type
migration, both halves of a co-family are supposed to map to the same
target file — `sgemm.f` and `dgemm.f` both become `qgemm.f` for
`kind16`. **Convergence** is the property that the two migrated bodies
are textually identical (up to whitespace, comments, and a handful of
documented cosmetic equivalences). When they aren't, that is a
**divergence** — a sign that the migrator missed something, or that
upstream's two halves genuinely differ in a way the migrator cannot
reconcile.

The migrator emits the D/Z half as canonical and discards the S/C half
silently when they converge. Divergences are surfaced through two
verifiers.

## Two reports

| command       | input                           | what it checks                                                                 |
|---------------|---------------------------------|--------------------------------------------------------------------------------|
| `diverge`     | recipe (no on-disk output)      | re-migrates every co-family pair in memory, compares with `_canonicalize_for_compare` (semantic — uppercase, KIND= stripping, REAL/CMPLX cast stripping, prefix-token folding) |
| `converge`    | recipe **+ migrated output dir** | reads on-disk canonical, re-migrates the S/C sibling in memory, compares with `_light_normalize` (whitespace + END-keyword merging only — no semantic equivalences) |

`diverge` answers "do the two upstream halves agree after migration in
principle?" `converge` answers "does the file the migrator actually
wrote to disk still agree with a fresh re-migration of its S/C
sibling?" The latter is the stricter check: anything that would have
required the canonicalizer's semantic rewrites to compare equal counts
as a real bug, because the migrator should have produced identical
output on its own.

## Generated from pipeline outputs

**Yes.** Both reports are first-class pipeline commands and consume
exactly the artifacts the migration pipeline produces:

```bash
# Stage all libraries (produces the canonical migrated tree)
uv run python -m migrator stage /tmp/stage-q --target kind16

# Per-library convergence against the staged outputs
uv run python -m migrator converge recipes/lapack.yaml \
    /tmp/stage-q/lapack --target kind16

# Or, without needing an output dir, the in-memory pair check:
uv run python -m migrator diverge recipes/lapack.yaml --target kind16
```

`migrator run` chains migrate → converge → verify → build, so
convergence is checked on every full-pipeline invocation and exits
non-zero when any pair diverges or any expected output is missing on
disk. CI gates on the exit code.

The reports themselves are returned as structured lists by
`pipeline.run_convergence_report` / `pipeline.run_divergence_report` —
see `src/migrator/pipeline.py` — so any custom rollup (per-target
matrix, per-library counts, machine-readable export) can be assembled
from those return values without re-running the migration.

## What gets suppressed

`converge` deliberately tolerates a small set of cosmetic
upstream-drift patterns so the report focuses on real bugs:

- **Comments** — stripped before comparison; D/Z and S/C halves
  routinely paraphrase the same prose differently.
- **Recipe-declared local renames** (`local_renames:` in the YAML) —
  applied to *both* halves at compare time only. Used for non-precision
  identifier drift like ScaLAPACK's `CR`/`CI` ↔ `ZR`/`ZI` in `pzlattrs`
  or `CONJTOPH` ↔ `TTOPH`.
- **Pure precision-prefix local-variable drift** — if a hunk differs
  only at S↔D / C↔Z positions inside identifier tokens, it's silently
  dropped (these are unclassified locals, not symbol-table entries —
  any properly-classified symbol has the same migrated name in both
  halves).

`diverge` is more permissive (it includes the full
`_canonicalize_for_compare` suite — KIND= suffix stripping, REAL/CMPLX
cast unwrapping, declaration-list sorting, etc.) because its purpose
is to find pairs that disagree *semantically*, not pairs that the
migrator failed to collapse.

## Statuses

`converge` entries carry a `status`:

- `diverged` — both files exist; their normalized texts differ.
- `missing` — the on-disk canonical isn't present in `output_dir`
  (usually a misconfigured recipe or a `skip_files` mistake).

`diverge` has no status field — every entry is a real semantic
disagreement.

## When divergences are expected

Some pairs are known to diverge for upstream reasons (e.g. a routine
where the C-half intentionally calls a different helper than the
Z-half). Those are tracked manually — there's no in-recipe whitelist
yet; they are reviewed each time the report is regenerated. If a
pair becomes a permanent expected divergence, the standard fix is to
make the recipe pin one half via `prefer_source:` and add a
`local_renames:` entry to fold the cosmetic drift, so the report
returns to clean.

## Related

- `doc/ARCHITECTURE.md` — overall pipeline.
- `pipeline._canonicalize_for_compare` / `pipeline._light_normalize` —
  the two normalization functions, with extensive in-source comments
  on every rule and why it exists.
- `pipeline._filter_precision_drift` — the S↔D / C↔Z prefix-folding
  step applied only to `converge`.
