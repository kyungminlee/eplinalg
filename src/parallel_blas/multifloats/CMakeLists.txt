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

# AVX2 SIMD double-double kernel for mgemm. 4-wide SoA — see
# mgemm_simd_kernel.h. Default ON for x86_64; the cmake check below
# disables it if -mavx2 isn't accepted by the compiler. Bench result
# on Raptor Lake: 3.8× over the scalar outer-product baseline single-
# threaded, 38× over migrated at OMP=4 s=1024.
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-mavx2 -mfma" _HAS_AVX2_FMA_FLAGS)
option(MBLAS_SIMD_DD
    "multifloats: AVX2 SIMD kernel for mgemm (4-wide SoA double-double)"
    ${_HAS_AVX2_FMA_FLAGS})
if(MBLAS_SIMD_DD)
    target_compile_definitions(${LIB_PREFIX}blas_parallel PRIVATE
        MBLAS_SIMD_DD=1)
    target_compile_options(${LIB_PREFIX}blas_parallel PRIVATE
        -mavx2 -mfma)
    # MR — register-tile rows for the SIMD micro-kernel.
    # Empirical sweep on Raptor Lake (i3-1315U) at OMP=1 s=1024:
    #   MR=1: 1.92 GFLOP/s            (baseline, 1 chain)
    #   MR=2: 3.49     (+82%, ILP win — two independent dd_mul chains)
    #   MR=3: 3.74     (+7%, sweet spot)
    #   MR=4: 3.78     (+1%, register pressure cancels parallelism gain)
    # Choose 3 by default — gains plateau there.
    set(MBLAS_SIMD_MR 3 CACHE STRING
        "multifloats: SIMD micro-kernel register-tile rows MR (1..4)")
    # NR_PAN: # of ymm-wide (4-col) panels per micro-kernel call.
    # 1 → NR=4 (current default). 2 → NR=8 (forces MR≤2 due to regs).
    set(MBLAS_SIMD_NR_PAN 1 CACHE STRING
        "multifloats: number of 4-lane B-panels per micro-kernel (1..2)")
    target_compile_definitions(${LIB_PREFIX}blas_parallel PRIVATE
        MGEMM_SIMD_MR=${MBLAS_SIMD_MR}
        MGEMM_SIMD_NR_PAN=${MBLAS_SIMD_NR_PAN})
endif()
