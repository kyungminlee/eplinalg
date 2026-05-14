# Hand-written parallel BLAS — kind10 (REAL(KIND=10) / 80-bit long double).
#
# Produces a STATIC archive `${LIB_PREFIX}blas_parallel` (= `eblas_parallel`
# for kind10) that supplies overlay implementations of selected BLAS
# entry points. Step 0 ships only a passthrough stub for egemm_ that
# forwards to the migrated archive's renamed symbol `egemm_migrated_`
# (produced by the objcopy rename in step 0b). This proves the
# Fortran-ABI C boundary end-to-end before any real kernel work.
#
# Linking model: see doc/parallel-blas-20260513.md §1. When
# PARALLEL_BLAS=ON, the public `${LIB_PREFIX}blas` target is rewired
# into an INTERFACE composite that links `_parallel` before `_serial`
# so overlay symbols win. This file builds the `_parallel` half only.

set(_kernels
    egemm.c
)

add_library(${LIB_PREFIX}blas_parallel STATIC ${_kernels})

target_compile_features(${LIB_PREFIX}blas_parallel PRIVATE c_std_99)
target_compile_options(${LIB_PREFIX}blas_parallel PRIVATE
    -O3
    -ffp-contract=fast
)
# -march=native is desirable but pinned per-host; gate so cross-builds
# and aarch64 CI don't break. Override with -DPARALLEL_BLAS_MARCH=...
set(PARALLEL_BLAS_MARCH "native" CACHE STRING
    "value passed to -march for parallel BLAS kernels (empty to disable)")
if(NOT "${PARALLEL_BLAS_MARCH}" STREQUAL "")
    target_compile_options(${LIB_PREFIX}blas_parallel PRIVATE
        -march=${PARALLEL_BLAS_MARCH})
endif()

# OpenMP comes in at step 4. Wire it now so the per-kernel #pragmas
# compile from the start.
if(OpenMP_C_FOUND)
    target_link_libraries(${LIB_PREFIX}blas_parallel PUBLIC OpenMP::OpenMP_C)
endif()

# Symbols this archive REFERS to (and the future composite resolves):
# the migrated archive renamed via objcopy in step 0b. We don't link
# anything here — the static archive can carry unresolved externals;
# they're satisfied at executable link time by the composite.
