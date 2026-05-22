# OpenBLAS-derived BLAS overlay (epopenblas).
#
# Builds a separate STATIC archive `${LIB_PREFIX}blas_epopenblas` that
# is NOT wired into the public ${LIB_PREFIX}blas composite. Used purely
# for A/B'ing OpenBLAS-D/Z-ports-at-kind10 against the migrated baseline
# via tests/blas_epopenblas/.
#
# Structure mirrors src/parallel_blas/ but stands alone:
#   src/epopenblas/CMakeLists.txt   ← this file
#   src/epopenblas/CHECKLIST.md     ← per-routine porting tracker
#   src/epopenblas/<TARGET_NAME>/CMakeLists.txt
#   src/epopenblas/<TARGET_NAME>/*.c
#
# Initial target: kind10. kind16 / multifloats can be added later by
# creating sibling target_<name>/ subdirs.

if(BASELINE_BUILD)
    message(STATUS "epopenblas: skipping under BASELINE_BUILD (no migrated target)")
    return()
endif()

set(_per_target_dir ${CMAKE_CURRENT_SOURCE_DIR}/${TARGET_NAME})
if(NOT IS_DIRECTORY ${_per_target_dir})
    message(STATUS "epopenblas: no kernels for target '${TARGET_NAME}' "
                   "at ${_per_target_dir} — skipping")
    return()
endif()

find_package(OpenMP COMPONENTS C)

set(EPOPENBLAS_MARCH "native" CACHE STRING
    "value passed to -march for epopenblas kernels (empty to disable)")

message(STATUS "epopenblas: building ${LIB_PREFIX}blas_epopenblas from ${_per_target_dir}")
add_subdirectory(${_per_target_dir} epopenblas_${TARGET_NAME})
