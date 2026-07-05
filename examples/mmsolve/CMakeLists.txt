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

# --- eplinalg extended kind10 e/y (long double) closure -------------------
# The MUMPS package carries the extended emumps solver, but its imported target
# pulls a per-precision parallel stack whose sibling packages eplinalg
# deliberately does NOT auto-load (see epblas-parallel/CONTEXT.md, "No
# sibling-deps inside eplinalg"): a consumer must find_package every piece of
# the closure itself, and — because find_package(emumps) validates that its
# referenced targets already exist — the whole closure must be found *before*
# the MUMPS package below.
#
# The arithmetic-free reference plumbing (scalapack/blacs/pblas/pbblas/ptzblas/
# blas/lapack) is named by the typed archives as their default backend; we load
# it only so the imported targets resolve. At link time MKL (ordered first)
# satisfies every one of those symbols, so no reference member is ever pulled.
foreach(_p eplinalgStdBlas eplinalgStdLapack eplinalgStdBlacs eplinalgStdPtzblas
           eplinalgStdPbblas eplinalgStdPblas eplinalgStdScalapack)
  find_package(${_p} REQUIRED)
endforeach()
# The kind10 e/y typed archives (long-double ScaLAPACK/BLACS/PBLAS/LAPACK/BLAS),
# in dependency order (each before its dependents):
foreach(_p eblas exblas eblacs eptzblas elapack
           epbblas epblas escalapack)
  find_package(${_p} REQUIRED)
endforeach()

# --- eplinalg (genuine MUMPS s/c/d/z + extended emumps + ordering) ---------
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

# Flatten MKL::MKL_SCALAPACK's transitive interface to absolute .so paths. The
# link line (below) must place MKL *ahead* of the in-tree static archives so MKL
# wins every symbol it provides — both the genuine s/c/d/z typed routines and the
# arithmetic-free plumbing (blacs_gridinit_, descinit_, numroc_, …). But CMake
# relocates library *targets* to the end of the link line, so passing MKL as the
# MKL::MKL_SCALAPACK target lets the in-tree reference plumbing (which the e/y
# typed archives name) capture those type-agnostic symbols first. Resolving MKL
# to raw file paths lets us pin it at the front, where raw items keep their order.
set(_mkl_files "")
set(_mkl_work MKL::MKL_SCALAPACK)
set(_mkl_seen "")
while(_mkl_work)
  list(POP_FRONT _mkl_work _cur)
  if(_cur IN_LIST _mkl_seen)
    continue()
  endif()
  list(APPEND _mkl_seen "${_cur}")
  if(TARGET ${_cur})
    get_target_property(_loc ${_cur} IMPORTED_LOCATION)
    if(_loc)
      list(APPEND _mkl_files "${_loc}")
    endif()
    get_target_property(_deps ${_cur} INTERFACE_LINK_LIBRARIES)
    if(_deps)
      list(APPEND _mkl_work ${_deps})
    endif()
  elseif(_cur MATCHES "\\.(so|a)$")
    list(APPEND _mkl_files "${_cur}")
  endif()
endwhile()
if(NOT _mkl_files)
  message(FATAL_ERROR "mmsolve: could not resolve MKL::MKL_SCALAPACK to .so files")
endif()
# Runtime path for the raw MKL .so (the MKL::MKL_SCALAPACK target would have added
# this itself; we replicate it since we bypass the target on the link line).
list(GET _mkl_files 0 _mkl_first)
get_filename_component(MMSOLVE_MKL_LIBDIR "${_mkl_first}" DIRECTORY)

# --- MPI ------------------------------------------------------------------
find_package(MPI REQUIRED COMPONENTS C Fortran)

# --- the executable -------------------------------------------------------
add_executable(mmsolve mmsolve.c mmio_min.c)

# We link the solver archives by file path (see below), which carries no usage
# requirements — so pull the MUMPS C headers (?mumps_c.h / ?mumps_struc.h) from
# the imported target's include interface explicitly.
target_include_directories(mmsolve PRIVATE
    $<TARGET_PROPERTY:eplinalg::emumps,INTERFACE_INCLUDE_DIRECTORIES>)

# Link order is load-bearing, and it is spelled out explicitly with archive file
# paths (MKL flattened to .so above; the eplinalg archives via $<TARGET_FILE:…>)
# rather than the packages' umbrella targets. Three reasons:
#
#   1. MKL must come FIRST. MKL's shared ScaLAPACK/BLACS provide the genuine
#      s/c/d/z typed routines (pdgetrf_, pzpotrf_, …) AND the arithmetic-free
#      plumbing (blacs_gridinit_, descinit_, numroc_, Cblacs_gridinfo). Placed
#      first, MKL (a .so) satisfies every one of those from the archives that
#      follow, so the in-tree reference members are never pulled and the whole
#      genuine + type-agnostic layer runs on MKL, as mandated. (A library *target*
#      here would be relocated by CMake to the end of the line and lose to the
#      in-tree plumbing — hence the flatten-to-.so above.)
#   2. The extended (e/y) typed archives name the in-tree reference backend
#      (eplinalg::scalapack/blacs/pblas/…) as a usage requirement, and the
#      reference ScaLAPACK archive itself DEFINES pdgetrf_/pzgetrf_. Referencing
#      the archives by file path drags NO interface, so no reference member can
#      sneak ahead of MKL and capture a genuine symbol.
#   3. The three genuine/extended solvers all cycle with the single mumps_common
#      runtime, so each vendor "*_full" umbrella wraps {solver, mumps_common} in
#      its own $<LINK_GROUP:RESCAN,…>; CMake forbids one target (mumps_common)
#      living in several distinct link groups, so we build one shared --start-group
#      by hand instead.
#
# The resulting line, in order:
#   MKL   — shared libs FIRST; win the genuine typed + all type-agnostic plumbing.
#   [group] genuine dz + sc + extended e/y solvers + shared mumps_common
#           — mutually recursive; pdgetrf_/pzgetrf_ (->MKL), pegetrf_/pypotrf_/
#             Cegamx2d (long double, ->e/y group), plumbing (->MKL).
#   [group] e/y typed long-double ScaLAPACK/LAPACK/BLACS/PBLAS/BLAS
#           — provide pegetrf_/pypotrf_/… (MKL has no long double); their
#             type-agnostic refs resolve to MKL above (a .so satisfies later refs).
#   [group] ordering: pord + esmumps + metis + scotch — mumps_common's analysis deps.
#   MPI, libm — libm supplies the long-double C math (cabsl, csqrtl, logl, …).
#
# MKL's OpenMPI BLACS is the netlib reference BLACS recompiled (it exports the BI_*
# internals), so the long-double reference routines share MKL's BLACS context
# representation: no ABI clash, no --allow-multiple-definition.
target_link_libraries(mmsolve PRIVATE
    -Wl,--no-as-needed                         # keep MKL .so even where nothing yet references them
    ${_mkl_files}                              # MKL FIRST: genuine s/c/d/z backend + type-agnostic plumbing
    -Wl,--start-group
      $<TARGET_FILE:eplinalg::dzmumps>       # double + double-complex genuine solver
      $<TARGET_FILE:eplinalg::scmumps>       # single + single-complex genuine solver
      $<TARGET_FILE:eplinalg::emumps>        # extended kind10 e/y (long double) solver
      $<TARGET_FILE:eplinalg::mumps_common>  # shared MUMPS runtime (Fortran common + C runtime)
    -Wl,--end-group
    -Wl,--start-group
      $<TARGET_FILE:eplinalg::escalapack>       # e/y long-double ScaLAPACK (pegetrf_, pypotrf_, …)
      # NUMROC/ICEIL/ILCM (formerly escalapack_tools) are now folded into
      # scalapack_common, which is dropped under MKL — MKL supplies the
      # type-agnostic numroc_/iceil_/ilcm_ here (linked first, above).
      $<TARGET_FILE:eplinalg::elapack>          # e/y long-double LAPACK
      $<TARGET_FILE:eplinalg::eblacs>           # e/y long-double BLACS (Cegamx2d, …)
      $<TARGET_FILE:eplinalg::epblas>           # e/y long-double PBLAS
      $<TARGET_FILE:eplinalg::epbblas>          # e/y long-double PB-BLAS
      $<TARGET_FILE:eplinalg::eptzblas>         # e/y long-double trapezoidal BLAS
      $<TARGET_FILE:eplinalg::eblas>            # e/y long-double BLAS
      $<TARGET_FILE:eplinalg::exblas>           # e/y long-double extra-precision BLAS
    -Wl,--end-group
    -Wl,--start-group
      $<TARGET_FILE:eplinalg::pord>             # MUMPS ordering deps (pulled by mumps_common)
      $<TARGET_FILE:eplinalg::esmumps>
      $<TARGET_FILE:eplinalg::metis>
      $<TARGET_FILE:eplinalg::scotch>
      $<TARGET_FILE:eplinalg::scotcherr>
    -Wl,--end-group
    MPI::MPI_Fortran MPI::MPI_C
    m)                             # long-double libm: cabsl, csqrtl, logl, log10l, lroundl

# The genuine solver archives are gfortran-compiled .F objects. Drive the final
# link with the Fortran driver so libgfortran/libquadmath come in, and disable
# LTO so the C main and the prebuilt Fortran objects don't fight over IR.
set_target_properties(mmsolve PROPERTIES LINKER_LANGUAGE Fortran)
target_compile_options(mmsolve PRIVATE -fno-lto)
target_link_options(mmsolve PRIVATE -fno-lto)

# We linked MKL as raw .so paths (not the MKL::MKL_SCALAPACK target), so add its
# runtime search path ourselves for a self-contained executable.
target_link_options(mmsolve PRIVATE "-Wl,-rpath,${MMSOLVE_MKL_LIBDIR}")
