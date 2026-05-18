# Hand-written parallel BLAS overlay for extended-precision targets.
#
# The overlay static archive `${LIB_PREFIX}blas_parallel` is ALWAYS
# built (no gate) so the kernels are continually exercised by CI
# alongside the migrated baseline.
#
# The `PARALLEL_BLAS` option (default ON) only controls whether the
# public CMake target `${LIB_PREFIX}blas` is rewired into an INTERFACE
# composite that links `_parallel` before `_serial` — see
# doc/parallel-blas-design.md §1. With `PARALLEL_BLAS=OFF`, the
# overlay archive still builds (and the consistency / fuzz / bench
# tests still run against it directly) but downstream consumers
# (${LIB_PREFIX}lapack etc.) link the serial migrated archive only.
# This is the knob for A/B'ing overall perf.

# PARALLEL_BLAS option is declared at the top of cmake/CMakeLists.txt
# so the BLAS section there can branch on it before this subdir runs.

if(BASELINE_BUILD)
    message(STATUS "parallel_blas: skipping under BASELINE_BUILD (no migrated target)")
    return()
endif()

set(_per_target_dir ${CMAKE_CURRENT_SOURCE_DIR}/${TARGET_NAME})
if(NOT IS_DIRECTORY ${_per_target_dir})
    message(STATUS "parallel_blas: no kernels for target '${TARGET_NAME}' "
                   "at ${_per_target_dir} — skipping")
    return()
endif()

# Detect both C and CXX so per-target subdirs (kind10/kind16 in C,
# multifloats in C++) can pick whichever they need.
find_package(OpenMP COMPONENTS C CXX)

# Shared knob across per-target subdirs.
set(PARALLEL_BLAS_MARCH "native" CACHE STRING
    "value passed to -march for parallel BLAS kernels (empty to disable)")

message(STATUS "parallel_blas: building ${LIB_PREFIX}blas_parallel from ${_per_target_dir}")
add_subdirectory(${_per_target_dir} parallel_blas_${TARGET_NAME})
