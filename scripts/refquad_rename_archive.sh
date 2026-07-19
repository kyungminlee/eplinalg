#!/bin/bash
# Rename every Fortran-mangled symbol (lowercase + trailing underscore)
# in a static archive from `<name>_` to `<name>_quad_`. Used as a
# POST_BUILD step on refblas_quad and reflapack_quad so they coexist in
# the same link with the standard `blas` / `lapack` archives.
#
# Why rename everything (not just public BLAS / LAPACK entries):
# upstream Fortran routines have many internal helpers (dladiv1_,
# dladiv2_, dlartg_, …) that aren't in any ref_quad_*.f90 interface
# block but are called from the public routines. If those internals
# stayed un-renamed, kind4 / kind8 baseline links would pull
# dladiv.f.o from BOTH archives and ld would error on multiple
# definition. Renaming uniformly is also kind to migrated test cycles
# that already pull std blas / lapack on their PUBLIC chain — the
# native helpers from std blas / lapack satisfy any leftover lsame_ /
# ilaenv_ / dlamch_ references.
#
# objcopy --redefine-syms renames both definitions and references, so
# cross-archive calls stay consistent: reflapack_quad's
# `U dnrm2_` reference is renamed to `U dnrm2_quad_`, then resolved
# from refblas_quad which now exports `T dnrm2_quad_`.
#
# Usage: refquad_rename_archive.sh <archive.a> [<sibling.a> ...] \
#                                   [-x <exclude.a> ...]
#
# Renames every Fortran-mangled symbol *defined* in <archive.a> or any
# sibling archive, EXCEPT names defined in any -x archive. References
# to renamed names inside <archive.a> are renamed in lockstep so
# cross-archive calls (reflapack_quad → refblas_quad's dnrm2_ →
# dnrm2_quad_) stay consistent.
#
# Exclude archives carry "common" routines that the migrated build
# expects to find under their ORIGINAL names — e.g. lapack_common
# provides iparam2stage_ (patched), la_constants_ep, mumps_abort_, …
# and these are referenced by reflapack_quad. Renaming those refs to
# `*_quad_` would leave them dangling in the migrated link cycle (no
# archive defines `iparam2stage_quad_`). The exclude list pulls the
# names out of the rename set so the references stay un-renamed and
# resolve normally.
#
# Snapshot mechanism: each archive's defined Fortran symbols are
# written to ``<archive>.quad-symbols`` *before* the rename. Sibling
# archives may already be renamed by the time this script runs on
# them, in which case `nm` returns `*_quad_` names; we read the
# pre-rename snapshot file instead.
#
# A stamp file alongside the archive marks completion so re-runs are
# no-ops (objcopy is non-idempotent without the guard — running twice
# yields `<name>_quad_quad_`).

set -e
ar=$1
shift
if [ -z "$ar" ] || [ ! -f "$ar" ]; then
    echo "usage: $0 <archive.a> [<sibling.a> ...] [-x <exclude.a> ...]" >&2
    exit 2
fi

# Split remaining args into siblings and excludes. Avoid bash arrays
# so the script runs under dash / sh too (the CMake call sites invoke
# `sh ${SCRIPT}` directly).
siblings_file=$(mktemp)
excludes_file=$(mktemp)
while [ $# -gt 0 ]; do
    case "$1" in
        -x)
            shift
            if [ -z "$1" ]; then
                echo "$0: -x needs an archive argument" >&2; exit 2
            fi
            printf '%s\n' "$1" >> "$excludes_file"
            ;;
        *)
            printf '%s\n' "$1" >> "$siblings_file"
            ;;
    esac
    shift
done

stamp="${ar}.quad-renamed"
snapshot="${ar}.quad-symbols"
if [ -f "$stamp" ] && [ "$stamp" -nt "$ar" ]; then
    exit 0
fi

# nm_pairs_of <archive> — emit `<orig> <new>` rename pairs for the
# archive's defined Fortran symbols (T = text, W = weak). The single
# owner of the rename policy: Fortran-mangled symbols (<name>_ ->
# append _quad_) and gfortran module-procedure symbols
# (__<mod>_MOD_<proc> -> append _quad). Module procs do not carry the
# trailing-underscore convention, but still need renaming so quad-
# promoted module bodies (la_xisnan, la_constants) do not collide with
# the std archive native versions of the same modules at link time.
nm_pairs_of() {
    nm --defined-only "$1" 2>/dev/null \
        | awk 'NF >= 3 {
                sec = $2; name = $3;
                if (sec != "T" && sec != "W") next;
                if (name ~ /^[a-z_][a-z0-9_]*_$/) {
                    bare = name; sub(/_$/, "", bare);
                    if (bare ~ /_quad$/) next;
                    print name " " bare "_quad_";
                } else if (name ~ /^__[a-z_][a-z0-9_]*_MOD_[a-z_][a-z0-9_]*$/) {
                    if (name ~ /_quad$/) next;
                    print name " " name "_quad";
                }
            }'
}

# Step 1 — scan: snapshot this archive's defined Fortran symbols
# BEFORE renaming. A sibling rename invoked later (when this archive
# has already been renamed to *_quad_) reads the pre-rename snapshot
# instead of nm-ing the archive directly.
nm_pairs_of "$ar" | sort -u > "$snapshot"

# Step 2 — build the rename set: this archive's defs + every sibling's
# defs (read from sibling snapshot if present, otherwise nm). objcopy
# --redefine-syms then renames T/W definitions AND U references in
# this archive that match any name in the rename set. References to
# symbols OUTSIDE this set (e.g. iparam2stage_, supplied by
# lapack_common from the migrated build, or libgfortran intrinsics
# whose names don't fit the pattern) stay un-renamed.
syms=$(mktemp)
exclset=$(mktemp)
trap 'rm -f "$syms" "$exclset" "$siblings_file" "$excludes_file"' EXIT

# Snapshot lines are now `<orig_symbol> <new_symbol>` pairs covering
# both Fortran-mangled (`<name>_` -> `<name>_quad_`) and gfortran
# module-procedure (`__<mod>_MOD_<proc>` -> `__<mod>_MOD_<proc>_quad`)
# patterns. Collect exclude *original* names and filter them out of
# the rename pairs.
{
    while IFS= read -r ex; do
        if [ -f "$ex" ]; then
            nm_pairs_of "$ex" | awk '{print $1}'
        fi
    done < "$excludes_file"
} | sort -u > "$exclset"

# Build the rename set: this archive's snapshot + sibling snapshots
# (read from sibling .quad-symbols if present, else nm-extracted),
# minus pairs whose original name is in the exclude set.
{
    cat "$snapshot"
    while IFS= read -r sib; do
        sib_snapshot="${sib}.quad-symbols"
        if [ -f "$sib_snapshot" ]; then
            cat "$sib_snapshot"
        elif [ -f "$sib" ]; then
            nm_pairs_of "$sib"
        fi
    done < "$siblings_file"
} | sort -u \
    | awk -v exclfile="$exclset" '
        BEGIN {
            while ((getline line < exclfile) > 0) excl[line] = 1;
        }
        $1 in excl { next }
        { print $1, $2 }
    ' > "$syms"

if [ ! -s "$syms" ]; then
    echo "$0: no Fortran-mangled symbols found in $ar" >&2
    touch "$stamp"
    exit 0
fi

n=$(wc -l < "$syms")
objcopy --redefine-syms="$syms" "$ar"
touch "$stamp"
echo "$0: renamed $n symbols in $(basename "$ar") (-> *_quad_)"
