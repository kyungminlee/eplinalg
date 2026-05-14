# Hand-written parallel BLAS overlay for extended-precision targets.
#
# The overlay static archive `${LIB_PREFIX}blas_parallel` is ALWAYS
# built (no gate) so the kernels are continually exercised by CI
# alongside the migrated baseline.
#
# The `PARALLEL_BLAS` option (default ON) only controls whether the
# public CMake target `${LIB_PREFIX}blas` is rewired into an INTERFACE
# composite that links `_parallel` before `_serial` — see
# doc/parallel-blas-20260513.md §1. With `PARALLEL_BLAS=OFF`, the
# overlay archive still builds (and the consistency / fuzz / bench
# tests still run against it directly) but downstream consumers
# (${LIB_PREFIX}lapack etc.) link the serial migrated archive only.
# This is the knob for A/B'ing overall perf.

option(PARALLEL_BLAS
    "Wire parallel BLAS overlay into the public \${LIB_PREFIX}blas target" ON)

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

# OpenMP is optional at step 0 (the passthrough stubs don't use it);
# step 4+ depends on it. Detect now so the per-target CMakeLists can
# link OpenMP::OpenMP_C unconditionally if it's needed there.
find_package(OpenMP COMPONENTS C)

message(STATUS "parallel_blas: building ${LIB_PREFIX}blas_parallel from ${_per_target_dir}")
add_subdirectory(${_per_target_dir} parallel_blas_${TARGET_NAME})
