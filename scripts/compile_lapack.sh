#!/usr/bin/env bash
# Parallel compile sweep over a migrated multifloats source tree.
# Usage: scripts/compile_lapack.sh [source-dir] [mod-dir]
#
# Reports total ok/fail counts and writes failing filenames to
# /tmp/lapack_failures.txt for follow-up analysis.

set -uo pipefail

SRC_DIR="${1:-/tmp/mf_lapack}"
MOD_DIR="${2:-/tmp}"
WORKERS="${WORKERS:-$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4)}"

REPO="$(cd "$(dirname "$0")"/.. && pwd)"
MULTIFLOATS_BUILD="$REPO/external/multifloats/build"

# Build the helper modules first (single-threaded; fast).
gfortran -c -ffree-line-length-132 \
    -I "$MULTIFLOATS_BUILD" -J "$MOD_DIR" \
    "$REPO/external/lapack-3.12.1/SRC/la_constants_mf.f90" \
    -o "$MOD_DIR/la_constants_mf.o" || exit 1
gfortran -c -ffree-line-length-132 \
    -I "$MULTIFLOATS_BUILD" -J "$MOD_DIR" \
    "$REPO/external/lapack-3.12.1/SRC/la_xisnan_mf.f90" \
    -o "$MOD_DIR/la_xisnan_mf.o" || exit 1

mkdir -p "$MOD_DIR/objs"
: > /tmp/lapack_failures.txt

compile_one() {
    f="$1"
    if [[ "$f" == *.f90 ]]; then
        gfortran -c -ffree-line-length-132 \
            -I "$2" -I "$3" \
            "$f" -o "$3/objs/$(basename "$f").o" 2>/dev/null
    else
        gfortran -c -ffixed-line-length-72 \
            -I "$2" -I "$3" \
            "$f" -o "$3/objs/$(basename "$f").o" 2>/dev/null
    fi
    if [ $? -eq 0 ]; then
        echo OK
    else
        echo "FAIL $(basename "$f")"
    fi
}
export -f compile_one

find "$SRC_DIR" -maxdepth 1 \( -name '*.f' -o -name '*.f90' \) -print0 \
    | xargs -0 -n1 -P "$WORKERS" -I {} \
        bash -c 'compile_one "$@"' _ {} "$MULTIFLOATS_BUILD" "$MOD_DIR" \
    > /tmp/compile_results.txt

ok=$(grep -c '^OK$' /tmp/compile_results.txt || echo 0)
fail=$(grep -c '^FAIL ' /tmp/compile_results.txt || echo 0)
grep '^FAIL ' /tmp/compile_results.txt | sed 's/^FAIL //' > /tmp/lapack_failures.txt
echo "Total: $ok ok / $fail failed"
