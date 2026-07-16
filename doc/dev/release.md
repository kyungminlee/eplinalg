# Release

1. Merge `develop` → `main`; create an annotated tag `vX.Y.Z` on the
   merge commit; push the tag.
2. The `v*` tag triggers `.github/workflows/release.yml`: stage each
   target, emit convergence reports, build the matrix, run the test
   slot, and publish 25 assets (one combined archive + 24 per-combo)
   with generated release notes.
3. List versions with `git tag --sort=-v:refname` — plain `git tag`
   sorts lexically (`v0.13.0` before `v0.8.0`).

Validation coverage of a release (CI matrix + off-CI MUMPS sweep):
[test.md](test.md).

Deprecations ride one release cycle: symbols renamed in `vX.Y` keep
deprecated aliases until `vX.(Y+1)`; note the window in the release
notes and remove the aliases after.
