# ---------------------------------------------------------------------------
# mmsolve — a standalone consumer of an *installed* eplinalg release.
#
# This is NOT part of the eplinalg build tree. It find_package()es a released
# eplinalg (the genuine MUMPS s/c/d/z solvers + ordering libraries), links Intel
# MKL as the ScaLAPACK/BLACS/PBLAS backend the genuine solvers expect, and MPI.
#
# Configure like:
#   cmake -S examples/mmsolve -B build \
#         -DCMAKE_PREFIX_PATH="<eplinalg-install>;<mkl-root>/lib/cmake" \
#         -DEPLINALG_MUMPS_PACKAGE=emumps \
#         -DMKL_MPI=openmpi
#
# The MKL BLACS flavor (MKL_MPI) MUST match the MPI the eplinalg release was
# built against: an OpenMPI release needs MKL_MPI=openmpi (libmkl_blacs_openmpi),
# an Intel-MPI release needs MKL_MPI=intelmpi. Mismatch => link/runtime breakage.
# ---------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.20)

# C for the program, Fortran because the genuine solver objects are gfortran-
# compiled .F — the final link must pull the gfortran runtime (see below).
project(mmsolve LANGUAGES C Fortran)

# The genuine dzmumps/scmumps solvers ship *inside* the ${LIB_PREFIX}mumps
# package of a release (e.g. kind10 -> emumps, kind16 -> qmumps). Point this at
# whichever package your install provides; the genuine targets are the same.
set(EPLINALG_MUMPS_PACKAGE "emumps" CACHE STRING
    "Name of the installed eplinalg MUMPS package to find (e.g. emumps, qmumps, dmumps)")

# --- eplinalg (genuine MUMPS s/c/d/z + ordering) --------------------------
find_package(${EPLINALG_MUMPS_PACKAGE} REQUIRED)

# --- Intel MKL: the s/c/d/z ScaLAPACK/BLACS/PBLAS backend -----------------
# Genuine dzmumps/scmumps reference only the *public* ScaLAPACK/BLACS API
# (pdgetrf_, pzpotrf_, blacs_gridinit_, …), which MKL exports — so MKL cleanly
# substitutes for the in-tree parallel stack (that stack is BUILD_INTERFACE-only
# in the export and never leaks into the consumer link).
set(MKL_INTERFACE   lp64       CACHE STRING "")   # 32-bit integers (matches MUMPS_INT)
set(MKL_THREADING   sequential CACHE STRING "")   # avoid oversubscribing under MPI
set(MKL_MPI         openmpi    CACHE STRING "MKL BLACS flavor; must match eplinalg's MPI")
set(ENABLE_SCALAPACK ON        CACHE BOOL   "")
find_package(MKL CONFIG REQUIRED)                 # provides MKL::MKL_SCALAPACK

# --- MPI ------------------------------------------------------------------
find_package(MPI REQUIRED COMPONENTS C Fortran)

# --- the executable -------------------------------------------------------
add_executable(mmsolve mmsolve.c mmio_min.c)

target_link_libraries(mmsolve PRIVATE
    eplinalg::dzmumps_full          # double + double-complex solver (+ mumps_common, ordering)
    eplinalg::scmumps_full          # single + single-complex solver
    MKL::MKL_SCALAPACK              # ScaLAPACK + PBLAS + BLACS + core MKL
    MPI::MPI_Fortran MPI::MPI_C)

# The genuine solver archives are gfortran-compiled .F objects. Drive the final
# link with the Fortran driver so libgfortran/libquadmath come in, and disable
# LTO so the C main and the prebuilt Fortran objects don't fight over IR.
set_target_properties(mmsolve PROPERTIES LINKER_LANGUAGE Fortran)
target_compile_options(mmsolve PRIVATE -fno-lto)
target_link_options(mmsolve PRIVATE -fno-lto)
