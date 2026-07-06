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
#         -DMKL_MPI=openmpi
#
# The MKL BLACS flavor (MKL_MPI) MUST match the MPI the eplinalg release was
# built against: an OpenMPI release needs MKL_MPI=openmpi (libmkl_blacs_openmpi),
# an Intel-MPI release needs MKL_MPI=intelmpi. Mismatch => link/runtime breakage.
# ---------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.20)

# C for the program, Fortran because the genuine solver objects are gfortran-
# compiled .F — the final link must pull the gfortran runtime (see below).
# CXX because the double-double C++ core (package `multifloats`) is compiled
# from C++ and its objects reference the C++ runtime (libstdc++); enabling CXX
# makes CMake link that runtime into the final Fortran-driver executable. The
# release ships the core as the portable, untagged libmultifloats.a (built
# non-LTO), so find_package(multifloats) resolves regardless of the consumer's
# exact C++ compiler — no LTO compiler-tag matching required.
project(mmsolve LANGUAGES C CXX Fortran)

# --- sequential (libmpiseq) release ---------------------------------------
# A `seq` release (MUMPS_LIBSEQ_RELEASE=ON) links the libmpiseq MPI stub rather
# than a real MPI. Its archives are still *compiled* against the (Intel) MPI
# mpi.h — libmpiseq.a is ABI-compatible — so we keep find_package(MPI) for the
# headers but replace the MPI *link* libraries with the installed mpiseq archive.
# Set -DMMSOLVE_SEQ=ON and point CMAKE_PREFIX_PATH at a seq install.
option(MMSOLVE_SEQ "Consume a sequential libmpiseq (seq-tag) eplinalg release" OFF)

# A seq release tags its MPI-dependent archives `seq` (MUMPS_LIBSEQ_RELEASE=ON),
# but this consumer still finds a real MPI below purely for the mpi.h headers —
# so the eplinalg Config packages' vendor probe would otherwise resolve the
# detected vendor's tag (e.g. intelmpi-2021.18) and fail to find the seq builds.
# EPLINALG_MPI_TAG forces the MPI flavor tag those Configs match against.
if(MMSOLVE_SEQ)
  set(EPLINALG_MPI_TAG "seq" CACHE STRING "Force the eplinalg MPI flavor tag (seq release)")
endif()

# --- eplinalg extended kind10 e/y (long double) closure -------------------
# The MUMPS package carries the extended eymumps solver, but its imported target
# pulls a per-precision parallel stack whose sibling packages eplinalg
# deliberately does NOT auto-load: a consumer must find_package every piece of
# the closure itself, and — because find_package(eymumps) validates that its
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
foreach(_p eyblas eyxblas eyblacs eyptzblas eylapack
           eypbblas eypblas eyscalapack)
  find_package(${_p} REQUIRED)
endforeach()

# --- eplinalg extended kind16 q/x (__float128) closure --------------------
# Same shape as e/y: the qxmumps solver names a per-precision parallel stack it
# does not auto-load, so we find every typed archive of the closure before the
# MUMPS package. quad_mpi(_f) are the first-party MPI bridges that register the
# __float128 MPI_Op / MPI_COMPLEX32 datatypes (transparent find_dependency's of
# qxblacs/qxscalapack/qxmumps, but we also link them directly below because
# mmsolve.c calls quad_mpi_init()).
foreach(_p qxblas qxxblas qxblacs qxptzblas qxlapack
           qxpbblas qxpblas qxscalapack)
  find_package(${_p} REQUIRED)
endforeach()

# --- eplinalg extended multifloats m/w (double-double) closure ------------
# The header-instantiated template ops of the double-double C++ core dedupe into
# each m-typed archive (weak/COMDAT), but its compiled extern "C" functions
# (adddd / cadddd / sindd / csqrtdd …) live only in the core archive and are
# shipped as the standalone `multifloats` package (linked below).
# multifloats_mpi(_f) bridge the double-double MPI_Op. The double-double
# numeric-constant / isnan helper modules (la_constants_mw / la_xisnan_mw)
# ride inside the m/w LAPACK archive itself.
foreach(_p mwblas mwxblas mwblacs mwptzblas mwlapack
           mwpbblas mwpblas mwscalapack)
  find_package(${_p} REQUIRED)
endforeach()

# --- first-party MPI bridges + double-double helper modules ---------------
# Linked directly (mmsolve.c calls quad_mpi_init() / multifloats_mpi_init()
# after MPI_Init — the genuine blacs_pinfo_ wins the link and never registers
# the q/x __float128 and m/w double-double custom MPI ops otherwise).
# multifloatsf is the double-double *Fortran* module (module `multifloats`,
# symbols `__multifloats_MOD_*`). The m/w typed archives reference it as
# *undefined* symbols (the eplinalg build links it BUILD_INTERFACE-only, so
# unlike the header-only C++ core it is NOT baked into the archives). The
# consumer must link it.
foreach(_p multifloatsf quad_mpi quad_mpi_f multifloats_mpi multifloats_mpi_f)
  find_package(${_p} REQUIRED)
endforeach()
# The double-double C++ core (package `multifloats`, target
# multifloats::multifloats). find_package(multifloatsf) already pulls it in via
# find_dependency, but request it explicitly so its archive can be placed on the
# link line: multifloatsf's module procedures call straight into its extern "C"
# adddd / cadddd / … which are defined nowhere else.
find_package(multifloats REQUIRED)

# --- eplinalg (genuine MUMPS s/c/d/z + extended e/q/m solvers + ordering) --
# All three typed MUMPS packages share the same genuine dz/sc, mumps_common and
# ordering leaves (factored into standalone packages), so finding all three is
# collision-free.
foreach(_p eymumps qxmumps mwmumps)
  find_package(${_p} REQUIRED)
endforeach()

# --- Intel MKL: the s/c/d/z ScaLAPACK/BLACS/PBLAS backend -----------------
# Genuine dzmumps/scmumps reference only the *public* ScaLAPACK/BLACS API
# (pdgetrf_, pzpotrf_, blacs_gridinit_, …), which MKL exports — so MKL cleanly
# substitutes for the in-tree parallel stack (that stack is BUILD_INTERFACE-only
# in the export and never leaks into the consumer link).
#
# A seq (libmpiseq) release is the exception: MKL's BLACS ships only in real-MPI
# flavors (libmkl_blacs_intelmpi/openmpi/…), each of which itself references the
# full real-MPI C runtime (MPI_Init_thread, MPI_Alltoallv, …) that libmpiseq does
# NOT provide — so pulling MKL ScaLAPACK/BLACS into a seq link is unresolvable.
# For seq we therefore take MKL as BLAS+LAPACK ONLY (MKL::MKL, no ScaLAPACK) —
# still the numeric workhorse (dgetrf_/dpotrf_ etc.), so MKL is not dropped — and
# let libmpiseq supply the ScaLAPACK/BLACS entry points (blacs_gridinit_, pdgetrf_,
# descinit_, numroc_, …) as its self-contained sequential stubs: at np=1 MUMPS
# routes the root node through the sequential LAPACK path and never calls them, so
# the stubs exist purely to satisfy the linker.
set(MKL_INTERFACE   lp64       CACHE STRING "")   # 32-bit integers (matches MUMPS_INT)
set(MKL_THREADING   sequential CACHE STRING "")   # avoid oversubscribing under MPI
set(MKL_MPI         openmpi    CACHE STRING "MKL BLACS flavor; must match eplinalg's MPI")
if(MMSOLVE_SEQ)
  set(ENABLE_SCALAPACK OFF CACHE BOOL "")         # no MKL BLACS/ScaLAPACK for seq
  set(_mkl_root_target MKL::MKL)                   # BLAS + LAPACK only
else()
  set(ENABLE_SCALAPACK ON  CACHE BOOL "")
  set(_mkl_root_target MKL::MKL_SCALAPACK)
endif()
find_package(MKL CONFIG REQUIRED)                 # provides MKL::MKL[_SCALAPACK]

# Flatten the chosen MKL root target's transitive interface to absolute .so paths.
# The link line (below) must place MKL *ahead* of the in-tree static archives so
# MKL wins every symbol it provides — both the genuine s/c/d/z typed routines and
# the arithmetic-free plumbing (blacs_gridinit_, descinit_, numroc_, …). But CMake
# relocates library *targets* to the end of the link line, so passing MKL as the
# target lets the in-tree reference plumbing (which the e/y typed archives name)
# capture those type-agnostic symbols first. Resolving MKL to raw file paths lets
# us pin it at the front, where raw items keep their order.
set(_mkl_files "")
set(_mkl_work ${_mkl_root_target})
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
  message(FATAL_ERROR "mmsolve: could not resolve ${_mkl_root_target} to .so files")
endif()
# Runtime path for the raw MKL .so (the MKL::MKL_SCALAPACK target would have added
# this itself; we replicate it since we bypass the target on the link line).
list(GET _mkl_files 0 _mkl_first)
get_filename_component(MMSOLVE_MKL_LIBDIR "${_mkl_first}" DIRECTORY)

# --- MPI ------------------------------------------------------------------
# Needed in every configuration for the mpi.h / mpif.h headers (mmsolve.c and
# the solver archives were compiled against them). For a seq release we take
# only the *include* dirs from MPI and link libmpiseq instead (below).
find_package(MPI REQUIRED COMPONENTS C Fortran)
if(MMSOLVE_SEQ)
  # libmpiseq ships in the seq install as its own package (export mpiseqTargets
  # -> config `mpiseq`, target eplinalg::mpiseq). It carries the MPI stub
  # symbols MKL BLACS / the genuine solvers reference at link time.
  find_package(mpiseq REQUIRED)
endif()

# --- the executable -------------------------------------------------------
add_executable(mmsolve mmsolve.c mmio_min.c)

# We link the solver archives by file path (see below), which carries no usage
# requirements — so pull the MUMPS C headers (?mumps_c.h / ?mumps_struc.h) from
# the imported targets' include interfaces explicitly. eymumps carries the
# genuine s/c/d/z + kind10 e/y headers; qxmumps the kind16 q/x; mwmumps the
# multifloats m/w.
target_include_directories(mmsolve PRIVATE
    $<TARGET_PROPERTY:eplinalg::eymumps,INTERFACE_INCLUDE_DIRECTORIES>
    $<TARGET_PROPERTY:eplinalg::qxmumps,INTERFACE_INCLUDE_DIRECTORIES>
    $<TARGET_PROPERTY:eplinalg::mwmumps,INTERFACE_INCLUDE_DIRECTORIES>)

# For a seq release we do NOT link the real-MPI targets (they'd drag the actual
# MPI runtime); we link libmpiseq instead. But mmsolve.c still #includes mpi.h,
# so pull the MPI headers as include-only usage requirements here.
if(MMSOLVE_SEQ)
  target_include_directories(mmsolve PRIVATE
      $<TARGET_PROPERTY:MPI::MPI_C,INTERFACE_INCLUDE_DIRECTORIES>
      $<TARGET_PROPERTY:MPI::MPI_Fortran,INTERFACE_INCLUDE_DIRECTORIES>)
  # A seq build does NOT link MKL's ScaLAPACK/BLACS (its BLACS flavors need a
  # real MPI runtime libmpiseq can't provide — see the MKL block above). MKL
  # supplies only the sequential BLAS/LAPACK numerics; the ScaLAPACK/BLACS layer
  # is libmpiseq's own MUMPS-libseq dummies, backed by the in-tree genuine
  # *reference* BLACS/ScaLAPACK (compiled -seq against libmpiseq) for the few
  # internals libmpiseq lacks. All of this goes in ONE --start-group with
  # libmpiseq FIRST — the order is load-bearing:
  #
  #   * libmpiseq's mpi.f.o is monolithic: it defines the MPI stubs AND dummy
  #     s/c/d/z entry points (blacs_gridinit_, descinit_, numroc_, pXgetrf_,
  #     pXpotrf_, pXscal_, infog2l_, …). Pulling it in for the MPI stubs drags
  #     those dummies too. They are the standard MUMPS-libseq stubs and are
  #     *correct* at np=1 (numroc_ returns the global n, descinit_ fills the
  #     descriptor, the root is factored by sequential LAPACK so pXgetrf_ is
  #     never called) — so we WANT them to win.
  #   * Scanning libmpiseq first means blacs_gridinit_/descinit_/numroc_/pXgetrf_
  #     are already defined when the reference archives are reached, so their
  #     colliding members (blacs_init_.c.o, descinit.f.o, numroc.f.o, pXgetrf.f.o)
  #     are NOT extracted — no multiple-definition, no --allow-multiple-definition.
  #   * The reference archives still contribute the members libmpiseq lacks:
  #       blacs      — BI_MringComb/BI_AuxBuff/BI_UpdateBuffs/… the typed
  #                    m/w · e/y · q/x BLACS archives reference (1×1-grid global
  #                    ops); their MPI refs resolve back to libmpiseq in-group.
  #       scalapack  — iceil_/ilcm_ and other type-agnostic helpers not in libseq.
  #       pblas/pbblas/ptzblas — any parallel-BLAS internals the above pull in.
  #     (Members whose symbol libmpiseq already defined are simply skipped.)
  set(_mmsolve_mpi_link
      -Wl,--start-group
        $<TARGET_FILE:eplinalg::mpiseq>       # FIRST: MPI stubs + s/c/d/z libseq dummies (win, correct at np=1)
        $<TARGET_FILE:eplinalg::blacs>        # reference BLACS: BI_* internals the typed BLACS need
        $<TARGET_FILE:eplinalg::scalapack>    # reference ScaLAPACK: iceil_/ilcm_/… libseq lacks
        $<TARGET_FILE:eplinalg::pblas>
        $<TARGET_FILE:eplinalg::pbblas>
        $<TARGET_FILE:eplinalg::ptzblas>
      -Wl,--end-group)
else()
  set(_mmsolve_mpi_link MPI::MPI_Fortran MPI::MPI_C)
endif()

# PT-Scotch is the parallel-ordering leaf: an MPI release (stacks ① ②) ships
# ptscotch so mmsolve can exercise ICNTL(28)=2 parallel analysis; a seq
# release (stack ③, no-PTSCOTCH) omits it. Link it only when the install
# actually provides it (it arrives via eplinalgPtscotch, a transparent
# find_dependency of the real-MPI flavors of the MUMPS packages above).
if(NOT MMSOLVE_SEQ AND TARGET eplinalg::ptscotch)
  set(_mmsolve_ptscotch
      $<TARGET_FILE:eplinalg::ptscotch>)
else()
  set(_mmsolve_ptscotch "")
endif()

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
      $<TARGET_FILE:eplinalg::eymumps>        # extended kind10 e/y (long double) solver
      $<TARGET_FILE:eplinalg::qxmumps>        # extended kind16 q/x (__float128) solver
      $<TARGET_FILE:eplinalg::mwmumps>        # extended multifloats m/w (double-double) solver
      $<TARGET_FILE:eplinalg::mumps_common>  # shared MUMPS runtime (Fortran common + C runtime)
    -Wl,--end-group
    -Wl,--start-group
      $<TARGET_FILE:eplinalg::eyscalapack>       # e/y long-double ScaLAPACK (pegetrf_, pypotrf_, …)
      # NUMROC/ICEIL/ILCM live in scalapack_common, which is dropped
      # under MKL — MKL supplies the type-agnostic numroc_/iceil_/ilcm_
      # here (linked first, above).
      $<TARGET_FILE:eplinalg::eylapack>          # e/y long-double LAPACK
      $<TARGET_FILE:eplinalg::eyblacs>           # e/y long-double BLACS (Cegamx2d, …)
      $<TARGET_FILE:eplinalg::eypblas>           # e/y long-double PBLAS
      $<TARGET_FILE:eplinalg::eypbblas>          # e/y long-double PB-BLAS
      $<TARGET_FILE:eplinalg::eyptzblas>         # e/y long-double trapezoidal BLAS
      $<TARGET_FILE:eplinalg::eyblas>            # e/y long-double BLAS
      $<TARGET_FILE:eplinalg::eyxblas>           # e/y long-double extra-precision BLAS
    -Wl,--end-group
    -Wl,--start-group
      $<TARGET_FILE:eplinalg::qxscalapack>       # q/x __float128 ScaLAPACK (pqgetrf_, pxpotrf_, …)
      $<TARGET_FILE:eplinalg::qxlapack>          # q/x __float128 LAPACK
      $<TARGET_FILE:eplinalg::qxblacs>           # q/x __float128 BLACS
      $<TARGET_FILE:eplinalg::qxpblas>           # q/x __float128 PBLAS
      $<TARGET_FILE:eplinalg::qxpbblas>          # q/x __float128 PB-BLAS
      $<TARGET_FILE:eplinalg::qxptzblas>         # q/x __float128 trapezoidal BLAS
      $<TARGET_FILE:eplinalg::qxblas>            # q/x __float128 BLAS
      $<TARGET_FILE:eplinalg::qxxblas>           # q/x __float128 extra-precision BLAS
    -Wl,--end-group
    -Wl,--start-group
      $<TARGET_FILE:eplinalg::mwscalapack>       # m/w double-double ScaLAPACK (pmgetrf_, pwpotrf_, …)
      $<TARGET_FILE:eplinalg::mwlapack>          # m/w double-double LAPACK
      $<TARGET_FILE:eplinalg::mwblacs>           # m/w double-double BLACS
      $<TARGET_FILE:eplinalg::mwpblas>           # m/w double-double PBLAS
      $<TARGET_FILE:eplinalg::mwpbblas>          # m/w double-double PB-BLAS
      $<TARGET_FILE:eplinalg::mwptzblas>         # m/w double-double trapezoidal BLAS
      $<TARGET_FILE:eplinalg::mwblas>            # m/w double-double BLAS
      $<TARGET_FILE:eplinalg::mwxblas>           # m/w double-double extra-precision BLAS
    -Wl,--end-group
    -Wl,--start-group
      # First-party MPI bridges + double-double helpers. quad_mpi(_f) register
      # the __float128 MPI ops/datatypes; multifloats_mpi(_f) the double-double
      # ones; multifloatsf the double-double Fortran arithmetic module itself
      # (__multifloats_MOD_*), referenced-undefined by all of the m/w archives.
      $<TARGET_FILE:eplinalg::multifloatsf>
      $<TARGET_FILE:eplinalg::quad_mpi>
      $<TARGET_FILE:eplinalg::quad_mpi_f>
      $<TARGET_FILE:eplinalg::multifloats_mpi>
      $<TARGET_FILE:eplinalg::multifloats_mpi_f>
      $<TARGET_FILE:multifloats::multifloats>  # C++ core: extern "C" adddd/cadddd/sindd/…
    -Wl,--end-group
    -Wl,--start-group
      $<TARGET_FILE:eplinalg::pord>             # MUMPS ordering deps (pulled by mumps_common)
      $<TARGET_FILE:eplinalg::esmumps>
      $<TARGET_FILE:eplinalg::metis>
      $<TARGET_FILE:eplinalg::scotch>
      ${_mmsolve_ptscotch}                      # PT-Scotch parallel ordering (MPI builds only)
    -Wl,--end-group
    ${_mmsolve_mpi_link}           # real MPI, or libmpiseq for a seq release
    m                              # long-double libm: cabsl, csqrtl, logl, log10l, lroundl
    quadmath                       # q/x __float128 math (from libquadmath)
    stdc++)                        # m/w double-double C++ core baked into the m-typed archives

# The genuine solver archives are gfortran-compiled .F objects. Drive the final
# link with the Fortran driver so libgfortran/libquadmath come in, and disable
# LTO so the C main and the prebuilt Fortran objects don't fight over IR.
set_target_properties(mmsolve PROPERTIES LINKER_LANGUAGE Fortran)
target_compile_options(mmsolve PRIVATE -fno-lto)
target_link_options(mmsolve PRIVATE -fno-lto)

# We linked MKL as raw .so paths (not the MKL::MKL_SCALAPACK target), so add its
# runtime search path ourselves for a self-contained executable.
target_link_options(mmsolve PRIVATE "-Wl,-rpath,${MMSOLVE_MKL_LIBDIR}")
