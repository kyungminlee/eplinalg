#!/usr/bin/env bash
# Compile sweep over migrated multifloats ScaLAPACK source tree.
#
# Usage:  scripts/compile_scalapack.sh [fortran-dir] [c-dir]
#
# Pipeline:
#   1. Run migrator for ScaLAPACK Fortran and C.
#   2. Compile Fortran files with gfortran.
#   3. Compile C files as C++ with mpicxx.
#   4. Report ok/fail counts.

set -uo pipefail

REPO="$(cd "$(dirname "$0")"/.. && pwd)"
F_DIR="${1:-/tmp/mf_scalapack}"
C_DIR="${2:-/tmp/mf_scalapack_c}"
MOD_DIR="/tmp/mf_scalapack_mod"
WORKERS="${WORKERS:-$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4)}"

MF_BRIDGE_INC="$REPO/external/multifloats/include"
MF_HH_DIR="$REPO/external"
MF_MOD="$REPO/external/multifloats/build"
MPICXX="${MPICXX:-mpicxx}"

export OMPI_CC="${OMPI_CC:-gcc-15}"
export OMPI_CXX="${OMPI_CXX:-g++-15}"
export OMPI_FC="${OMPI_FC:-gfortran-15}"

# ---- Step 1: Migrate ----

echo "[1/4] Migrating ScaLAPACK Fortran..."
(cd "$REPO/src" && uv run --project .. python -m migrator migrate \
    "$REPO/recipes/scalapack.yaml" "$F_DIR" --target multifloats) 2>&1 | tail -5

echo "[2/4] Migrating ScaLAPACK C..."
(cd "$REPO/src" && uv run --project .. python -m migrator migrate \
    "$REPO/recipes/scalapack_c.yaml" "$C_DIR" --target multifloats) 2>&1 | tail -5

# ---- Step 2: Build helper modules ----

mkdir -p "$MOD_DIR/objs"
gfortran -c -ffree-line-length-132 \
    -I "$MF_MOD" -J "$MOD_DIR" \
    "$REPO/external/lapack-3.12.1/SRC/la_constants_mf.f90" \
    -o "$MOD_DIR/la_constants_mf.o" || exit 1
gfortran -c -ffree-line-length-132 \
    -I "$MF_MOD" -J "$MOD_DIR" \
    "$REPO/external/lapack-3.12.1/SRC/la_xisnan_mf.f90" \
    -o "$MOD_DIR/la_xisnan_mf.o" || exit 1

# ---- Step 3: Compile Fortran ----

echo "[3/4] Compiling Fortran files..."
f_ok=0; f_fail=0
: > /tmp/scalapack_f_failures.txt
for f in "$F_DIR"/*.f "$F_DIR"/*.f90; do
    [ -f "$f" ] || continue
    base="$(basename "$f")"
    if [[ "$f" == *.f90 ]]; then
        gfortran -c -ffree-line-length-132 \
            -I "$MF_MOD" -I "$MOD_DIR" \
            "$f" -o "$MOD_DIR/objs/$base.o" 2>/dev/null
    else
        gfortran -c -ffixed-line-length-72 \
            -I "$MF_MOD" -I "$MOD_DIR" \
            "$f" -o "$MOD_DIR/objs/$base.o" 2>/dev/null
    fi
    if [ $? -eq 0 ]; then
        f_ok=$((f_ok+1))
    else
        f_fail=$((f_fail+1))
        echo "$base" >> /tmp/scalapack_f_failures.txt
    fi
done
echo "  Fortran: $f_ok ok / $f_fail failed"

# ---- Step 4: Compile C as C++ ----

echo "[4/4] Compiling C files as C++ with mpicxx..."
c_ok=0; c_fail=0
: > /tmp/scalapack_c_failures.txt
cp "$MF_BRIDGE_INC/multifloats_bridge.h" "$C_DIR/"
for f in "$C_DIR"/*.c; do
    [ -f "$f" ] || continue
    base="$(basename "$f")"
    out=$($MPICXX -O2 -std=c++17 -DAdd_ \
        -fpermissive -Wno-write-strings \
        -I "$MF_HH_DIR" -I "$MF_BRIDGE_INC" -I "$C_DIR" \
        -include multifloats_bridge.h \
        -c "$f" -o "$MOD_DIR/objs/$base.o" 2>&1)
    if [ $? -eq 0 ]; then
        c_ok=$((c_ok+1))
    else
        c_fail=$((c_fail+1))
        echo "$base" >> /tmp/scalapack_c_failures.txt
    fi
done
echo "  C/C++:   $c_ok ok / $c_fail failed"

echo ""
echo "ScaLAPACK total: $((f_ok + c_ok)) ok / $((f_fail + c_fail)) failed"
