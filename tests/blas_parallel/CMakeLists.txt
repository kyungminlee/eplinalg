# Parallel BLAS overlay test suite.
#
# Two layers (only consistency lands in step 1):
#
#  1. Consistency / fuzz — `${LIB_PREFIX}blas_parallel` (overlay) vs
#     `${LIB_PREFIX}blas_migrated` (objcopy-renamed migrated archive)
#     at the working precision. Tolerance 10 ulp relative error
#     (target eps × 10, per doc/parallel-blas-20260513.md §4).
#
#  2. Bench (step 2) — geometric size sweep, JSON output.
#
# Both link the test exe against:
#   • the composite `${LIB_PREFIX}blas` (so `egemm_` resolves to the
#     overlay's C kernel) AND
#   • the renamed `${LIB_PREFIX}blas_migrated` (so `egemm_migrated_`
#     resolves to the Fortran migrated archive).
# Archive semantics ensure no symbol collision: the migrated archive's
# non-renamed members (eaxpy_, etc.) are only pulled in if referenced.

if(NOT HAVE_REAL10)
    message(STATUS "tests/blas_parallel: skipping — HAVE_REAL10 required")
    return()
endif()

# Shared utility module: env-driven seed + case count, random shape /
# parameter generation, error metric.
add_library(blas_parallel_fuzz_common STATIC
    common/fuzz_util.f90
)
fortran_module_layout(blas_parallel_fuzz_common)

function(add_consistency_test driver_src)
    get_filename_component(_name ${driver_src} NAME_WE)
    add_executable(${_name} ${driver_src})
    target_link_libraries(${_name} PRIVATE
        blas_parallel_fuzz_common
        ${LIB_PREFIX}blas          # composite — overlay wins for egemm_
        ${LIB_PREFIX}blas_migrated # renamed archive — provides egemm_migrated_
    )
    add_test(NAME blas_parallel_${_name} COMMAND ${_name})
    # Pass env defaults; user can override on the ctest invocation.
    set_tests_properties(blas_parallel_${_name} PROPERTIES
        ENVIRONMENT "BLAS_FUZZ_CASES=200")
endfunction()

add_consistency_test(consistency/fuzz_egemm.f90)
