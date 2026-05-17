# Parallel BLAS overlay test suite.
#
# Generates per-target fuzz_util + per-(target,routine) fuzz / bench
# drivers from fypp templates. Naming convention:
#   common/fuzz_util_body.fypp                  ← module body, parameterized
#   common/target_<name>/fuzz_util.fypp         ← per-target shim
#   consistency/fuzz_gemm_{real,complex}_body.fypp
#   consistency/target_<name>/fuzz_<prefix>gemm.fypp
#   bench/bench_gemm_real_body.fypp
#   bench/target_<name>/bench_<prefix>gemm.fypp
#
# add_fypp_source (defined in tests/CMakeLists.txt) drives fypp.

if(NOT (HAVE_REAL10 OR HAVE_REAL16))
    message(STATUS "tests/blas_parallel: skipping — no extended precision available")
    return()
endif()

set(_common_target_dir
    ${CMAKE_CURRENT_SOURCE_DIR}/common/target_${TARGET_NAME})
set(_consist_target_dir
    ${CMAKE_CURRENT_SOURCE_DIR}/consistency/target_${TARGET_NAME})
set(_bench_target_dir
    ${CMAKE_CURRENT_SOURCE_DIR}/bench/target_${TARGET_NAME})

if(NOT IS_DIRECTORY ${_common_target_dir})
    message(STATUS "tests/blas_parallel: no fuzz_util shim for target "
                   "'${TARGET_NAME}' — skipping")
    return()
endif()

# Generate fuzz_util.f90 for this target.
add_fypp_source(_fuzz_util_src
    ${_common_target_dir}/fuzz_util.fypp
    ${CMAKE_CURRENT_SOURCE_DIR}/common)

add_library(blas_parallel_fuzz_common STATIC ${_fuzz_util_src})
fortran_module_layout(blas_parallel_fuzz_common)
if(NEEDS_MULTIFLOATS AND TARGET multifloatsf)
    # multifloats target's fuzz_util uses real64x2/cmplx64x2 from the
    # upstream Fortran module.
    target_link_libraries(blas_parallel_fuzz_common PUBLIC
        $<BUILD_INTERFACE:multifloatsf>)
endif()

function(add_consistency_test shim_path)
    get_filename_component(_name ${shim_path} NAME_WE)
    add_fypp_source(_src ${shim_path}
        ${CMAKE_CURRENT_SOURCE_DIR}/consistency)
    add_executable(${_name} ${_src})
    target_link_libraries(${_name} PRIVATE
        blas_parallel_fuzz_common
        ${LIB_PREFIX}blas
        ${LIB_PREFIX}blas_migrated)
    add_test(NAME blas_parallel_${_name} COMMAND ${_name})
    set_tests_properties(blas_parallel_${_name} PROPERTIES
        ENVIRONMENT "BLAS_FUZZ_CASES=200")
endfunction()

function(add_bench_test shim_path)
    get_filename_component(_name ${shim_path} NAME_WE)
    add_fypp_source(_src ${shim_path}
        ${CMAKE_CURRENT_SOURCE_DIR}/bench)
    add_executable(${_name} ${_src})
    target_link_libraries(${_name} PRIVATE
        blas_parallel_fuzz_common
        ${LIB_PREFIX}blas
        ${LIB_PREFIX}blas_migrated)
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/bench_reports)
    add_test(NAME blas_parallel_${_name} COMMAND ${_name}
             WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bench_reports)
    set_tests_properties(blas_parallel_${_name} PROPERTIES LABELS "bench")
endfunction()

if(IS_DIRECTORY ${_consist_target_dir})
    file(GLOB _consist_shims CONFIGURE_DEPENDS
        ${_consist_target_dir}/fuzz_*.fypp)
    foreach(_s ${_consist_shims})
        add_consistency_test(${_s})
    endforeach()
endif()

if(IS_DIRECTORY ${_bench_target_dir})
    file(GLOB _bench_shims CONFIGURE_DEPENDS
        ${_bench_target_dir}/bench_*.fypp)
    foreach(_s ${_bench_shims})
        add_bench_test(${_s})
    endforeach()
endif()

# C perf harnesses — kernel-isolated timing without the Fortran bench's
# per-iter realloc/copy asymmetry. Optional, glob per-target.
set(_perf_target_dir
    ${CMAKE_CURRENT_SOURCE_DIR}/perf/target_${TARGET_NAME})
function(add_perf_test src_path)
    get_filename_component(_name ${src_path} NAME_WE)
    add_executable(${_name} ${src_path})
    target_link_libraries(${_name} PRIVATE
        ${LIB_PREFIX}blas
        ${LIB_PREFIX}blas_migrated)
    # Keep overlay and migrated entry points close in the linked binary.
    # Without --gc-sections, the full eblas_parallel archive's other
    # kernels get pulled in, pushing overlay's egemv ~100 KB away from
    # migrated's egemv. That gap costs ~10% on the timed function via
    # iTLB/page-walk overhead at the call site and would mask the actual
    # kernel comparison.
    target_compile_options(${_name} PRIVATE -ffunction-sections -fdata-sections)
    target_link_options(${_name} PRIVATE LINKER:--gc-sections)
    add_test(NAME blas_parallel_${_name} COMMAND ${_name})
    set_tests_properties(blas_parallel_${_name} PROPERTIES LABELS "perf")
endfunction()
if(IS_DIRECTORY ${_perf_target_dir})
    file(GLOB _perf_srcs CONFIGURE_DEPENDS
        ${_perf_target_dir}/perf_*.c)
    foreach(_s ${_perf_srcs})
        add_perf_test(${_s})
    endforeach()
endif()
