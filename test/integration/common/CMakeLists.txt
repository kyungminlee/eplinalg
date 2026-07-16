# Per-target conversion module shared by every test framework.
#
# The conv module exposes a uniform API (q2t_r / t2q_r / q2t_c / t2q_c)
# that the templated test wrappers call to convert between the
# universal quad reference precision and whichever precision the
# target's migrated routines actually use. Each per-target conv module
# implements those four converters with whatever logic the target
# needs (identity for kind16, intrinsic cast for kind10, double-double
# split for multifloats).
#
# TARGET_NAME comes from target_config.cmake.

set(_conv_dir ${CMAKE_CURRENT_SOURCE_DIR}/target_${TARGET_NAME})
if(NOT IS_DIRECTORY ${_conv_dir})
    message(STATUS "tests/common: no conv module for target '${TARGET_NAME}' "
                   "at ${_conv_dir} — skipping")
    return()
endif()

add_library(target_conv STATIC ${_conv_dir}/target_conv.f90)
fortran_test_module_layout(target_conv)
if(NEEDS_MULTIFLOATS AND TARGET multifloatsf)
    target_link_libraries(target_conv PUBLIC
        $<BUILD_INTERFACE:multifloatsf>)
endif()
