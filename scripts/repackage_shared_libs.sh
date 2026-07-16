#!/bin/bash
# Repackage one extended-precision pair's static archives as per-archive
# shared libraries, for MKL (LP64) + Intel MPI consumers.
#
#   usage: repackage_shared_libs.sh <pair> <install-lib-dir> <output-dir>
#
#     pair             qx | ey | mw
#     install-lib-dir  lib/ of an eplinalg install tree holding the pair's
#                      static archives (PIC objects — the shipped archives
#                      are compiled -fPIC)
#     output-dir       where lib*.so land (created; scratch in <out>/pruned)
#
# Environment: mpicc on PATH targeting the MPI ABI the archives were
# tagged for; MKLROOT set (MKL stays the provider of every
# standard-precision symbol).  CC/MPICC override the drivers.
#
# One .so per static archive, named after the archive base (libqxmumps.so,
# libqxblas.so, libmumps_common.so, libepblacs_common.so, ...).
# Exceptions and fixups:
#
#   * ordering engines (pord/metis/scotch/ptscotch/esmumps) are merged INTO
#     libmumps_common.so; a GENERATED version script exports exactly
#     def(mumps_common) + the renamed engine entry points that the pair
#     MUMPS archive references directly (scotchf*_mumps_, metis_mumps_*),
#     hiding every bare engine symbol (SCOTCH_context*, gk_*, PORD, ...)
#   * Scotch error-handler prune: drop library_error.c.o + library_errcom.c.o
#     from scotch/ptscotch copies, keep library_error_exit.c.o (the object
#     static links select)
#   * cross-archive duplicate members get ONE home:
#       {c}{c}dotc/{c}{c}dotu  -> keep in ptzblas (it calls them), prune the
#                                 scalapack copies
#       {c}symv/{c}syr         -> keep in lapack, prune the ptzblas copies
#                                 (nothing in ptzblas calls them)
#   * family symbols (blas, lapack, scalapack, mumps APIs) are all
#     exported — no filtering on family libs
#   * every link: -Wl,--no-undefined, -soname, -rpath $ORIGIN; deps on the
#     already-built .so set + MKL/gfortran/quadmath are --as-needed so
#     NEEDED lists stay minimal.  NEVER --allow-multiple-definition.
#
# Two flags are MANDATORY on every .so link here and cannot be supplied by
# a consumer link line afterwards (see doc/user/mkl-coexistence.md,
# "Repackaging archives as shared libraries"):
#   --no-define-common   never allocate Fortran COMMON blocks in a .so
#   -z now               eager PLT binding; lazy resolution corrupts live
#                        FP vector state inside the library's own PLT
set -uo pipefail

if [ $# -ne 3 ]; then
  sed -n '2,12p' "$0" | sed 's/^# \{0,1\}//'
  exit 2
fi
PAIR=$1; L=$2; OUTDIR=$3
case $PAIR in qx|ey|mw) ;; *) echo "unknown pair '$PAIR' (want qx|ey|mw)"; exit 2 ;; esac
[ -d "$L" ] || { echo "install lib dir '$L' not found"; exit 2; }
: "${MKLROOT:?set MKLROOT (source the oneAPI environment)}"

CC=${CC:-gcc}
MPICC=${MPICC:-mpicc}
MKL="-L$MKLROOT/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64"
TAIL_COMMON="-lgfortran -lquadmath -lstdc++ -lpthread -lm -ldl"

D=$OUTDIR; mkdir -p "$D"
P=$OUTDIR/pruned; rm -rf "$P"; mkdir -p "$P"
R=${PAIR:0:1}; C=${PAIR:1:1}
BUILT=""   # .so files built so far, newest first (link deps)

# member_of <archive> <symbol> -> archive member defining symbol
member_of() {
  nm -A --defined-only "$1" 2>/dev/null | \
    awk -F: -v s="$2" '{ n=split($0,f," "); if (f[n]==s) { print $2; exit } }'
}

# link_so <soname-base> <driver> <archives...>
link_so() {
  local base=$1 drv=$2; shift 2
  local so="$D/lib$base.so"
  # --no-define-common: never allocate Fortran COMMON blocks in a .so;
  # they must resolve as imports of their real definer (Intel MPI's
  # mpipriv*/mpifcmb* sentinel commons -> libmpifort.so, quad/multifloats
  # op handles -> the C bridge .so).  A second allocation silently forks
  # MPI_IN_PLACE/MPI_BOTTOM detection across the boundary.
  # -z now: eager PLT binding.  Lazy first-call resolution through glibc's
  # _dl_runtime_resolve corrupts live FP state in calls whose arguments/
  # returns are multi-limb reals in vector registers, inside the library's
  # OWN PLT slots — an eager executable does not protect a lazy library
  # (observed: flaky, ASLR-layout-deterministic wrong double-double
  # correction limbs in the w family at np=4; -z now suppresses it 100%).
  $drv -shared -o "$so" -Wl,-soname,"lib$base.so" -Wl,-rpath,'$ORIGIN' \
    -Wl,--whole-archive "$@" -Wl,--no-whole-archive \
    ${DEFCOMMON--Wl,--no-define-common} -Wl,-z,now ${EXTRA:-} \
    -Wl,--push-state,--as-needed $BUILT $MKL $TAIL_COMMON -Wl,--pop-state \
    -Wl,--no-undefined 2> "$D/err-$base.log"
  local rc=$?
  if [ $rc -ne 0 ]; then
    echo "== $PAIR lib$base.so: LINK FAILED rc=$rc"
    sort -u "$D/err-$base.log" | head -12
    exit 1
  fi
  BUILT="$so $BUILT"
  echo "   $PAIR lib$base.so  $(ls -lh "$so" | awk '{print $5}')  NEEDED: $(readelf -d "$so" | awk '/NEEDED/{gsub(/[\[\]]/,"",$5); printf "%s ",$5}')"
}

echo "== $PAIR: building per-archive .so set in $D"

# ---- 1. runtime bridges -----------------------------------------------------
case $PAIR in
  qx)
    link_so quad_mpi   "$MPICC" "$L"/libquad_mpi-*.a
    link_so quad_mpi_f "$MPICC" "$L"/libquad_mpi_f-*.a
    ;;
  mw)
    link_so multifloats       "$CC"    "$L"/libmultifloats-nolto.a
    link_so multifloatsf      "$CC"    "$L"/libmultifloatsf-*.a
    link_so multifloats_mpi   "$MPICC" "$L"/libmultifloats_mpi-*.a
    link_so multifloats_mpi_f "$MPICC" "$L"/libmultifloats_mpi_f-*.a
    ;;
esac

# ---- 2. ep_ commons ---------------------------------------------------------
link_so epblacs_common     "$MPICC" "$L"/libepblacs_common-*.a
link_so epptzblas_common   "$CC"    "$L"/libepptzblas_common-*.a
link_so eppblas_common     "$CC"    "$L"/libeppblas_common-*.a
# ep_sltimer00_ (SLtimer state) has no B/D definer anywhere and no user
# outside this archive — let this .so allocate it (skip --no-define-common;
# ld's -d does NOT override an earlier --no-define-common, they are
# independent flags and inhibit wins)
DEFCOMMON="" \
link_so epscalapack_common "$CC"    "$L"/libepscalapack_common-*.a

# ---- 3. family compute stack ------------------------------------------------
link_so ${PAIR}blas   "$CC" "$L"/lib${PAIR}blas-*.a
link_so xblas_common  "$CC" "$L"/libxblas_common.a
link_so ${PAIR}xblas  "$CC" "$L"/lib${PAIR}xblas.a
link_so ${PAIR}lapack "$CC" "$L"/lib${PAIR}lapack-*.a

# ptzblas: prune {c}symv/{c}syr (lapack is their home; ptzblas never calls them)
cp "$L"/lib${PAIR}ptzblas-*.a "$P/ptzblas.a"
for s in ${C}symv_ ${C}syr_; do
  m=$(member_of "$P/ptzblas.a" "$s")
  [ -n "$m" ] && ar d "$P/ptzblas.a" "$m"
done
link_so ${PAIR}ptzblas "$CC" "$P/ptzblas.a"

link_so ${PAIR}blacs  "$MPICC" "$L"/lib${PAIR}blacs-*.a
link_so ${PAIR}pblas  "$CC"    "$L"/lib${PAIR}pblas.a
link_so ${PAIR}pbblas "$CC"    "$L"/lib${PAIR}pbblas-*.a

# scalapack: prune {c}{c}dotc/{c}{c}dotu (ptzblas is their home and calls them)
cp "$L"/lib${PAIR}scalapack-*.a "$P/scalapack.a"
for s in ${C}${C}dotc_ ${C}${C}dotu_; do
  m=$(member_of "$P/scalapack.a" "$s")
  [ -n "$m" ] && ar d "$P/scalapack.a" "$m"
done
link_so ${PAIR}scalapack "$CC" "$P/scalapack.a"

# ---- 4. libmumps_common.so = glue + ordering engines ------------------------
cp "$L"/libptscotch_mumps-*.a "$P/ptscotch.a"
ar d "$P/ptscotch.a" library_error.c.o library_errcom.c.o
cp "$L"/libscotch_mumps.a "$P/scotch.a"
ar d "$P/scotch.a" library_error.c.o 2>/dev/null || true
ar d "$P/scotch.a" library_errcom.c.o 2>/dev/null || true

# version script: def(mumps_common) + engine syms the pair mumps archive
# references directly; everything else (bare engine internals) local
MC=$(ls "$L"/libmumps_common-*.a)
MU=$(ls "$L"/lib${PAIR}mumps-*.a)
nm --defined-only "$MC" 2>/dev/null | awk 'NF==3 && $2 ~ /^[TDBRGWV]$/{print $3}' | sort -u > "$P/def-mc"
nm --defined-only "$P/ptscotch.a" "$P/scotch.a" "$L"/libpord_mumps.a "$L"/libesmumps_mumps.a "$L"/libmetis_mumps.a 2>/dev/null | \
  awk 'NF==3 && $2 ~ /^[TDBRGWV]$/{print $3}' | sort -u > "$P/def-eng"
nm --undefined-only "$MU" 2>/dev/null | awk '{print $NF}' | sort -u > "$P/und-mu"
comm -12 "$P/und-mu" "$P/def-eng" > "$P/eng-needed"
{
  echo "{"
  echo "global:"
  cat "$P/def-mc" "$P/eng-needed" | sort -u | sed 's/$/;/'
  echo "local: *;"
  echo "};"
} > "$P/mumps_common.map"
echo "   $PAIR mumps_common version script: $(wc -l < "$P/def-mc") glue syms + $(wc -l < "$P/eng-needed") engine entry points exported"

EXTRA="-Wl,--version-script=$P/mumps_common.map" \
link_so mumps_common "$MPICC" "$MC" "$L"/libpord_mumps.a "$P/ptscotch.a" "$P/scotch.a" "$L"/libesmumps_mumps.a "$L"/libmetis_mumps.a
unset EXTRA

# ---- 5. pair MUMPS ----------------------------------------------------------
link_so ${PAIR}mumps "$MPICC" "$MU"

# ---- validation ---------------------------------------------------------
bare=$(readelf --dyn-syms -W "$D/libmumps_common.so" | \
       awk 'NF>=8 && $7 != "UND" && $6 == "DEFAULT" && ($5=="GLOBAL"||$5=="WEAK") {print $8}' | \
       grep -icE '^(SCOTCH_|METIS_|gk_|esmumps|errexit|buildelement|errorprint)' || true)
fam=$(readelf --dyn-syms -W "$D/lib${PAIR}blas.so" | grep -cE " (${R}gemm_|${C}gemm_)$" || true)
slp=$(readelf --dyn-syms -W "$D/lib${PAIR}scalapack.so" | grep -cE " p${R}getrf_$" || true)
mmp=$(readelf --dyn-syms -W "$D/lib${PAIR}mumps.so" | grep -ci "mumps" || true)
# cross-.so duplicate exports (defined, default visibility) — should be none
for so in $BUILT; do
  readelf --dyn-syms -W "$so" | \
    awk 'NF>=8 && $7 != "UND" && $6 == "DEFAULT" && $5 == "GLOBAL" {print $8}'
done | sort | uniq -d | grep -vE '^(__bss_start|_edata|_end)$' > "$D/dup-exports"
# no .so may allocate the Intel MPI Fortran sentinel commons
sentinels=0
for so in $BUILT; do
  n=$(readelf --dyn-syms -W "$so" | awk 'NF>=8 && $7!="UND" && $5=="GLOBAL"' | grep -cE ' (mpipriv|mpifcmb)' || true)
  sentinels=$((sentinels + n))
done
# every .so must carry DF_BIND_NOW except the designated COMMON owner too
noweager=0
for so in $BUILT; do
  readelf -d "$so" | grep -q BIND_NOW || { echo "   FATAL: $so lacks BIND_NOW"; noweager=1; }
done
echo "== $PAIR SUMMARY: bare-engine-exports=$bare  ${PAIR}blas gemm exports=$fam  p${R}getrf exported=$slp  mumps-syms in lib${PAIR}mumps.so=$mmp  cross-so dup exports=$(wc -l < "$D/dup-exports")  sentinel-common-allocs=$sentinels"
[ -s "$D/dup-exports" ] && { echo "   DUPS:"; head -10 "$D/dup-exports"; }
[ "$sentinels" -ne 0 ] && { echo "   FATAL: MPI sentinel commons allocated in a .so"; exit 1; }
[ "$noweager" -ne 0 ] && exit 1
echo "REPACKAGE $PAIR: DONE"
