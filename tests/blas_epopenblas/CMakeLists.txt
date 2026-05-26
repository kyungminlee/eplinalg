# epopenblas overlay test suite.
#
# Mirrors tests/blas_parallel/ but compares the epopenblas archive
# against ${LIB_PREFIX}blas_epopenblas_migrated instead of
# ${LIB_PREFIX}blas_migrated. The Fortran-ABI shape of every routine
# is identical, so the consistency / bench *_body.fypp files are
# REUSED in place from tests/blas_parallel/ via fypp's -I include path.
#
# Per-target shims live under:
#   common/target_<name>/fuzz_util.fypp
#   consistency/target_<name>/fuzz_<routine>.fypp
#   bench/target_<name>/bench_<routine>.fypp
#
# Each bench routine is registered TWICE: once at OMP_NUM_THREADS=1 and
# once at OMP_NUM_THREADS=4 (CTest LABEL: epopenblas-bench-omp1 /
# epopenblas-bench-omp4). Fuzz tests carry the LABEL epopenblas-fuzz.

if(NOT (HAVE_REAL10 OR HAVE_REAL16))
    message(STATUS "tests/blas_epopenblas: skipping — no extended precision available")
    return()
endif()

set(_common_target_dir
    ${CMAKE_CURRENT_SOURCE_DIR}/common/target_${TARGET_NAME})
set(_consist_target_dir
    ${CMAKE_CURRENT_SOURCE_DIR}/consistency/target_${TARGET_NAME})
set(_bench_target_dir
    ${CMAKE_CURRENT_SOURCE_DIR}/bench/target_${TARGET_NAME})

if(NOT IS_DIRECTORY ${_common_target_dir})
    message(STATUS "tests/blas_epopenblas: no fuzz_util shim for target "
                   "'${TARGET_NAME}' — skipping")
    return()
endif()

# Body fypp files live under tests/blas_parallel/ — same Fortran ABI,
# nothing in them references "parallel" specifically, so we reuse them
# directly. The shim sets the fypp include path to its body dir.
set(_parallel_common_body
    ${CMAKE_CURRENT_SOURCE_DIR}/../blas_parallel/common)
set(_parallel_consist_body
    ${CMAKE_CURRENT_SOURCE_DIR}/../blas_parallel/consistency)
set(_parallel_bench_body
    ${CMAKE_CURRENT_SOURCE_DIR}/../blas_parallel/bench)

# fuzz_util shim for this target. The shim file itself lives here;
# the body it #:includes is the parallel suite's fuzz_util_body.fypp.
add_fypp_source(_fuzz_util_src
    ${_common_target_dir}/fuzz_util.fypp
    ${_parallel_common_body})

add_library(blas_epopenblas_fuzz_common STATIC ${_fuzz_util_src})
fortran_module_layout(blas_epopenblas_fuzz_common)
if(NEEDS_MULTIFLOATS AND TARGET multifloatsf)
    target_link_libraries(blas_epopenblas_fuzz_common PUBLIC
        $<BUILD_INTERFACE:multifloatsf>)
endif()

function(add_ep_consistency_test shim_path)
    get_filename_component(_name ${shim_path} NAME_WE)
    add_fypp_source(_src ${shim_path} ${_parallel_consist_body})
    add_executable(ep_${_name} ${_src})
    # eblas archive supplies lsame_/xerbla_ that the migrated L2 routines
    # call internally; listed last so the overlay/migrated archives'
    # routines are picked up first.
    target_link_libraries(ep_${_name} PRIVATE
        blas_epopenblas_fuzz_common
        ${LIB_PREFIX}blas_epopenblas
        ${LIB_PREFIX}blas_epopenblas_migrated
        blas)
    add_test(NAME blas_epopenblas_${_name} COMMAND ep_${_name})
    set_tests_properties(blas_epopenblas_${_name} PROPERTIES
        ENVIRONMENT "BLAS_FUZZ_CASES=200"
        LABELS "epopenblas-fuzz")
endfunction()

function(add_ep_bench_test shim_path)
    get_filename_component(_name ${shim_path} NAME_WE)
    add_fypp_source(_src ${shim_path} ${_parallel_bench_body})
    add_executable(ep_${_name} ${_src})
    target_link_libraries(ep_${_name} PRIVATE
        blas_epopenblas_fuzz_common
        ${LIB_PREFIX}blas_epopenblas
        ${LIB_PREFIX}blas_epopenblas_migrated
        blas)
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/reports/epopenblas)
    foreach(_threads 1 4)
        set(_test_name blas_epopenblas_${_name}_omp${_threads})
        set(_json ${CMAKE_BINARY_DIR}/reports/epopenblas/${_name}_omp${_threads}.json)
        add_test(NAME ${_test_name} COMMAND ep_${_name}
                 WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/reports/epopenblas)
        set_tests_properties(${_test_name} PROPERTIES
            ENVIRONMENT "OMP_NUM_THREADS=${_threads};BLAS_BENCH_OUT=${_json}"
            LABELS "epopenblas-bench;epopenblas-bench-omp${_threads}")
    endforeach()
endfunction()

if(IS_DIRECTORY ${_consist_target_dir})
    file(GLOB _consist_shims CONFIGURE_DEPENDS
        ${_consist_target_dir}/fuzz_*.fypp)
    foreach(_s ${_consist_shims})
        add_ep_consistency_test(${_s})
    endforeach()
endif()

if(IS_DIRECTORY ${_bench_target_dir})
    file(GLOB _bench_shims CONFIGURE_DEPENDS
        ${_bench_target_dir}/bench_*.fypp)
    foreach(_s ${_bench_shims})
        add_ep_bench_test(${_s})
    endforeach()
endif()

# ── C perf harnesses ─────────────────────────────────────────────────
# Kernel-isolated timing — no per-iter realloc/copy asymmetry, with
# -ffunction-sections / --gc-sections so overlay's and migrated's
# symbols stay within a few KB of each other (no iTLB-spread artifact).
#
# The perf_*.c sources live under tests/blas_parallel/perf/target_<TGT>/
# and reference `<routine>_` / `<routine>_migrated_` symbols. They are
# overlay-agnostic — same source, different archives linked. We reuse
# them here but only build the ones whose routine is in the epopenblas
# overlay set so the link succeeds without dragging in unported symbols.
set(_perf_target_dir
    ${CMAKE_CURRENT_SOURCE_DIR}/../blas_parallel/perf/target_${TARGET_NAME})
function(add_ep_perf_test src_path)
    get_filename_component(_name ${src_path} NAME_WE)  # perf_<routine>
    add_executable(ep_${_name} ${src_path})
    target_link_libraries(ep_${_name} PRIVATE
        ${LIB_PREFIX}blas_epopenblas
        ${LIB_PREFIX}blas_epopenblas_migrated
        blas)
    target_compile_options(ep_${_name} PRIVATE
        -ffunction-sections -fdata-sections)
    target_link_options(ep_${_name} PRIVATE LINKER:--gc-sections)
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/reports/epopenblas)
    foreach(_threads 1 4)
        set(_test_name blas_epopenblas_${_name}_omp${_threads})
        set(_json ${CMAKE_BINARY_DIR}/reports/epopenblas/${_name}_omp${_threads}.json)
        add_test(NAME ${_test_name} COMMAND ep_${_name}
                 WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/reports/epopenblas)
        set_tests_properties(${_test_name} PROPERTIES
            ENVIRONMENT "OMP_NUM_THREADS=${_threads};BLAS_PERF_JSON=${_json}"
            LABELS "epopenblas-perf;epopenblas-perf-omp${_threads}")
    endforeach()
endfunction()
if(IS_DIRECTORY ${_perf_target_dir} AND DEFINED EPOPENBLAS_OVERLAY_ROUTINES)
    foreach(_r ${EPOPENBLAS_OVERLAY_ROUTINES})
        set(_perf_src ${_perf_target_dir}/perf_${_r}.c)
        if(EXISTS ${_perf_src})
            add_ep_perf_test(${_perf_src})
        else()
            # Some routines have a .cpp harness (multifloats uses C++).
            set(_perf_src_cpp ${_perf_target_dir}/perf_${_r}.cpp)
            if(EXISTS ${_perf_src_cpp})
                add_ep_perf_test(${_perf_src_cpp})
            endif()
        endif()
    endforeach()
endif()
