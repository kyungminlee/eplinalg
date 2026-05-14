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

# Opt-in: inline __multf3/__addtf3 via vendored soft-fp template
# headers under external/libgcc-softfp/. Saves the libgcc external
# call but balloons kernel code size (~700+1266 insns per qmul+qadd).
# Off by default; flip with -DQBLAS_INLINE_SOFTFP=ON at cmake time.
option(QBLAS_INLINE_SOFTFP
    "kind16: inline __multf3/__addtf3 from soft-fp headers" OFF)
option(QBLAS_USE_FMAQ
    "kind16: use fmaq() instead of separate __multf3 + __addtf3" OFF)
if(QBLAS_USE_FMAQ)
    target_compile_definitions(${LIB_PREFIX}blas_parallel PRIVATE
        QBLAS_USE_FMAQ=1)
endif()
if(QBLAS_INLINE_SOFTFP)
    # Soft-fp headers reference each other and the per-arch
    # `sfp-machine.h` via plain `#include "..."`. Add both the root
    # external dir AND the i386/64 subdir to the include path so the
    # nested includes resolve.
    set(_softfp_root
        ${CMAKE_SOURCE_DIR}/external/libgcc-softfp)
    if(NOT IS_DIRECTORY ${_softfp_root})
        message(FATAL_ERROR
            "QBLAS_INLINE_SOFTFP=ON but ${_softfp_root} is missing — "
            "ensure migrator stage copied external/libgcc-softfp/")
    endif()
    target_include_directories(${LIB_PREFIX}blas_parallel PRIVATE
        ${_softfp_root}
        ${_softfp_root}/config/i386/64)
    target_compile_definitions(${LIB_PREFIX}blas_parallel PRIVATE
        QBLAS_INLINE_SOFTFP=1)
    # Soft-fp emits ISO-C extensions that warn under strict modes.
    # Suppress; the headers are upstream LGPL code we don't modify.
    target_compile_options(${LIB_PREFIX}blas_parallel PRIVATE
        -Wno-error -w)
endif()

if(OpenMP_C_FOUND)
    target_link_libraries(${LIB_PREFIX}blas_parallel PUBLIC OpenMP::OpenMP_C)
endif()
