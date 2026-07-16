# tests/pbblas — TODO

Resolved items have moved to `CHANGELOG.md`.

## By-design exclusions

### `pbdtrget` / `pbdtrsrt` / `pbdtrst1` (and complex variants)

Internal helpers used by `pbdtran` / `pbdtrnv` to gather, sort, and
stitch the diagonal-block representation of a column-block matrix.
Their inputs are intermediate buffers laid out by the caller and have
no standalone semantic meaning — testing them in isolation would
require fabricating a contrived data layout that may not match what the
caller actually produces. These are exercised transitively through the
`pb[dz]tran` / `pb[dz]trnv` test families (baseline + replicated +
non-square-grid + NZ-offset). Reopen if a regression motivates digging
in.

## Cross-tree follow-ups

(None. The NUMROC / ICEIL / ILCM cross-tree dependency was resolved
on the migrator side — see CHANGELOG entry "scalapack_tools recipe".)
