This directory is a workspace for type migration of BLAS, LAPACK, BLACS, and ScaLAPACK in respective directories.
It contains  scripts that extracts all the necessary information, does the type migration, and checks for correctness, and produce reports,
compile, and then build the libraries in extended precision.

## Per-library sidecars

Each library may carry an optional `<library>/` sidecar directory alongside its
top-level `<library>.yaml` recipe. The sidecar holds library-specific
artifacts — override tables, constant maps, or line-level manifests — that
would clutter the YAML.

### `mumps/keep-kind.manifest`

A line-level **keep-kind manifest**: a list of `path:lineno` entries naming
`DOUBLE PRECISION` declarations in MUMPS's `d*`/`z*` files that must *not*
be promoted to the target precision. This is MUMPS-specific because MUMPS
overloads `DOUBLE PRECISION` to mean two different things inside its
per-arithmetic files:

1. **Working precision** (to be promoted).
2. **Arithmetic-agnostic DP** — timing, flop counters, memory/load
   estimates, and external-library ABI buffers (Scotch, `MPI_WTIME`, …) —
   which must remain `DOUBLE PRECISION` in every build. Promoting it would
   break the MPI_WTIME return-kind contract and similar ABIs, and would
   create divergence vs. the s/c-migrated output.

A purely syntactic migration rule cannot tell these two apart: both are
spelled `DOUBLE PRECISION`. The manifest is the exception list.

**Generation rule (one rule, applied only to paired files).** For each
`d*`/`z*` file, a DP declaration is keep-kind iff the same declaration
(whitespace-normalized) also appears in the `s*`/`c*` partner file. If
both the s-author and d-author agreed on DP, that agreement is meaningful:
it signals intent independent of working precision.

Regenerate after any MUMPS upstream bump:

```
python scripts/mumps_sweep_keep_kind.py
```

**Shared files are not in the manifest.** Un-prefixed files
(e.g. `tools_common.F`, `mumps_load.F`, `lr_stats.F`) are linked into all
four arithmetic builds and therefore cannot depend on working precision;
every floating-point declaration in them is arithmetic-agnostic by
construction. Rather than protecting them line-by-line, the recipe lists
each shared file in `copy_files` so the whole file is copied through
unchanged. The sweep script reports the list of shared files it detects
(run with `--verbose`) — paste that list into the recipe.

**Not** the right mechanism for:

- *Whole files* to skip or copy through unchanged — use `skip_files` /
  `copy_files` in the YAML (e.g. `double_linked_list.F` and every shared
  un-prefixed file belong in `copy_files`).
- *Routines* with no arithmetic variant — same, via `skip_files`.

See `doc/guide/developer.md` §3 for the definition of divergence that motivates
this mechanism.
