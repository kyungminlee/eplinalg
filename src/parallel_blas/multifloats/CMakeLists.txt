# Hand-written parallel BLAS — multifloats target (double-double).
#
# Type: upstream `float64x2` (2 IEEE doubles, ~32-digit precision).
# Each DD op is an error-free transform: ~8 scalar flops per add,
# ~12 per mul. Arithmetic-bound; OpenMP scales nearly linearly.
#
# Uses the multifloats C API (adddd/muldd/fmadd/cmuldd/etc.) from
# <multifloats/float64x2.h>. The multifloats library is provided by
# the upstream FetchContent target `multifloats`.

set(_kernels
    mgemm.cpp
    wgemm.cpp
)

set(_routines "")
foreach(_k ${_kernels})
    string(REGEX REPLACE "\\.cpp$" "" _r "${_k}")
    list(APPEND _routines "${_r}")
endforeach()
set(PARALLEL_BLAS_OVERLAY_ROUTINES "${_routines}" CACHE INTERNAL
    "BLAS routines provided by the parallel overlay for this target")

add_library(${LIB_PREFIX}blas_parallel STATIC ${_kernels})

# C++17 for std::aligned_alloc + multifloats's header requirements.
target_compile_features(${LIB_PREFIX}blas_parallel PRIVATE cxx_std_17)
target_compile_options(${LIB_PREFIX}blas_parallel PRIVATE
    -O3
    # multifloats forbids -ffast-math (breaks IEEE FMA semantics that
    # the error-free transforms rely on). Use -ffp-contract=on, not
    # =fast — fast lets the compiler fold across statements.
    -ffp-contract=on
)
if(NOT "${PARALLEL_BLAS_MARCH}" STREQUAL "")
    target_compile_options(${LIB_PREFIX}blas_parallel PRIVATE
        -march=${PARALLEL_BLAS_MARCH})
endif()

# multifloats C library provides adddd / muldd / fmadd / cmuldd / ...
target_link_libraries(${LIB_PREFIX}blas_parallel PUBLIC multifloats)

if(OpenMP_CXX_FOUND)
    target_link_libraries(${LIB_PREFIX}blas_parallel PUBLIC OpenMP::OpenMP_CXX)
elseif(OpenMP_C_FOUND)
    target_link_libraries(${LIB_PREFIX}blas_parallel PUBLIC OpenMP::OpenMP_C)
endif()
