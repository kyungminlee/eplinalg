# Hand-written parallel BLAS — kind16 (REAL(KIND=16) / __float128).
#
# Produces static archive `${LIB_PREFIX}blas_parallel` (= `qblas_parallel`
# for kind16). __float128 has no hardware support on x86-64; every
# add / multiply lowers to a libquadmath call (~hundreds of cycles).
# That makes the arithmetic the bottleneck, not memory, so cache
# blocking yields little but OpenMP scaling is excellent.

set(_kernels
    qgemm.c
    xgemm.c
)

set(_routines "")
foreach(_k ${_kernels})
    string(REGEX REPLACE "\\.c$" "" _r "${_k}")
    list(APPEND _routines "${_r}")
endforeach()
set(PARALLEL_BLAS_OVERLAY_ROUTINES "${_routines}" CACHE INTERNAL
    "BLAS routines provided by the parallel overlay for this target")

add_library(${LIB_PREFIX}blas_parallel STATIC ${_kernels})

target_compile_features(${LIB_PREFIX}blas_parallel PRIVATE c_std_99)
target_compile_options(${LIB_PREFIX}blas_parallel PRIVATE
    -O3
    -ffp-contract=fast
)
if(NOT "${PARALLEL_BLAS_MARCH}" STREQUAL "")
    target_compile_options(${LIB_PREFIX}blas_parallel PRIVATE
        -march=${PARALLEL_BLAS_MARCH})
endif()

# Quadmath provides __float128 arithmetic + the complex helpers
# (conjq, cabsq, ...). Link explicitly — libquadmath is a separate .so
# from libgcc.
target_link_libraries(${LIB_PREFIX}blas_parallel PUBLIC quadmath)

if(OpenMP_C_FOUND)
    target_link_libraries(${LIB_PREFIX}blas_parallel PUBLIC OpenMP::OpenMP_C)
endif()
