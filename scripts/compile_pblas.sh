#!/usr/bin/env bash
# Compile sweep over migrated multifloats PBLAS source tree.
#
# Usage:  scripts/compile_pblas.sh [output-dir]
#
# Pipeline:
#   1. Run migrator for PTZBLAS, PBBLAS (Fortran) and PBLAS (C).
#   2. Copy PTZBLAS + PBBLAS Fortran objects into PBLAS output.
#   3. Compile Fortran kernels with gfortran.
#   4. Compile C entry points + PTOOLS as C++ with mpicxx.
#   5. Report ok/fail counts.

set -uo pipefail

REPO="$(cd "$(dirname "$0")"/.. && pwd)"
OUT_DIR="${1:-/tmp/mf_pblas_full}"
WORKERS="${WORKERS:-$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4)}"

MF_BRIDGE_INC="$REPO/external/multifloats/include"
MF_HH_DIR="$REPO/external"
MF_MOD="$REPO/external/multifloats/build"
MF_MPI_SRC="$REPO/external/multifloats/src/multifloats_mpi.cpp"
MPICC="${MPICC:-mpicc}"
MPICXX="${MPICXX:-mpicxx}"

export OMPI_CC="${OMPI_CC:-gcc-15}"
export OMPI_CXX="${OMPI_CXX:-g++-15}"
export OMPI_FC="${OMPI_FC:-gfortran-15}"

PTZBLAS_DIR="/tmp/mf_ptzblas"
PBBLAS_DIR="/tmp/mf_pbblas"
PBLAS_C_DIR="/tmp/mf_pblas_c"

echo "[1/5] Migrating PTZBLAS (Fortran kernels)..."
(cd "$REPO/src" && uv run --project .. python -m migrator migrate \
    "$REPO/recipes/ptzblas.yaml" "$PTZBLAS_DIR" --target multifloats) 2>&1 | tail -3

echo "[2/5] Migrating PBBLAS (Fortran kernels)..."
(cd "$REPO/src" && uv run --project .. python -m migrator migrate \
    "$REPO/recipes/pbblas.yaml" "$PBBLAS_DIR" --target multifloats) 2>&1 | tail -3

echo "[3/5] Migrating PBLAS C entry points + PTOOLS..."
(cd "$REPO/src" && uv run --project .. python -m migrator migrate \
    "$REPO/recipes/pblas.yaml" "$PBLAS_C_DIR" --target multifloats) 2>&1 | tail -3

echo "[4/5] Compiling Fortran kernels with gfortran..."
mkdir -p "$OUT_DIR/objs"
f_ok=0; f_fail=0
for d in "$PTZBLAS_DIR" "$PBBLAS_DIR"; do
    for f in "$d"/*.f; do
        [ -f "$f" ] || continue
        gfortran -c -ffixed-line-length-72 \
            -I "$MF_MOD" \
            "$f" -o "$OUT_DIR/objs/$(basename "$f").o" 2>/dev/null
        if [ $? -eq 0 ]; then
            f_ok=$((f_ok+1))
        else
            f_fail=$((f_fail+1))
        fi
    done
done
echo "  Fortran: $f_ok ok / $f_fail failed"

echo "[5/5] Compiling C/C++ entry points + PTOOLS with mpicxx..."
BLACS_DEFS="-DCSAMEF77=1 -DInt=int -DUseMpi2=1 -DAdd_"
c_ok=0; c_fail=0
# Copy bridge header into PBLAS output for override includes
cp "$MF_BRIDGE_INC/multifloats_bridge.h" "$PBLAS_C_DIR/"
for f in "$PBLAS_C_DIR"/*.c; do
    [ -f "$f" ] || continue
    out=$($MPICXX -O2 -std=c++17 -DAdd_ \
        -fpermissive -Wno-write-strings \
        -I "$MF_HH_DIR" -I "$MF_BRIDGE_INC" -I "$PBLAS_C_DIR" \
        -include multifloats_bridge.h \
        -c "$f" -o "$OUT_DIR/objs/$(basename "$f").o" 2>&1)
    if [ $? -eq 0 ]; then
        c_ok=$((c_ok+1))
    else
        c_fail=$((c_fail+1))
    fi
done
echo "  C/C++:   $c_ok ok / $c_fail failed"

total_ok=$((f_ok + c_ok))
total_fail=$((f_fail + c_fail))
echo ""
echo "PBLAS total: $total_ok ok / $total_fail failed"
