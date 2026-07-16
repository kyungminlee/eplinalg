# tests/blacs — TODO

Resolved items have moved to `CHANGELOG.md`. No live items.

## Cross-tree follow-ups (NOT to be edited from here)

- The migrator's BLACS recipe stages `blacs_setup_.c` but no public
  Fortran callers in the test tree need it. If a future test wants to
  exercise the alternate setup path, add an interface block in
  `target_blacs_body.fypp`.
- `target_blacs` does not declare wrappers for `cgesd2d` / `zgesd2d`
  etc. (the original-precision entry points). All wrappers go through
  the migrated `q` / `x` / `i` entries. Add original-precision wrappers
  if a regression test is wanted to confirm parity.
