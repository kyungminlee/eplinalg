# Top-level test directory — added by cmake/CMakeLists.txt when
# BUILD_TESTING=ON. Each library has its own subdir; we only descend
# into ones whose precision library was actually built (i.e., the
# library was migrated for this target) AND whose CMakeLists.txt
# exists in the staging copy.

# fypp is used to expand the per-framework target_*_body.fypp template
# into per-target target_*.f90 sources. Available via `uv run fypp`
# from the project venv.
find_program(FYPP fypp REQUIRED)

# Helper for per-target wrapper shims:
#   add_fypp_source(<out_var> <shim.fypp> <body_dir>)
# generates <build>/.../<basename>.f90 from <shim.fypp>, with -I
# pointing at <body_dir> so #:include resolves the shared body file.
function(add_fypp_source out_var shim_path body_dir)
    get_filename_component(_stem ${shim_path} NAME_WE)
    set(_out ${CMAKE_CURRENT_BINARY_DIR}/${_stem}.f90)
    # Depend on every .fypp in body_dir so edits to the shared body
    # (or any future #:include sibling) trigger regeneration.
    file(GLOB _body_fypps CONFIGURE_DEPENDS ${body_dir}/*.fypp)
    add_custom_command(
        OUTPUT ${_out}
        COMMAND ${FYPP} -I ${body_dir} ${shim_path} ${_out}
        DEPENDS ${shim_path} ${_body_fypps}
        COMMENT "fypp ${shim_path} -> ${_out}"
        VERBATIM
    )
    set(${out_var} ${_out} PARENT_SCOPE)
endfunction()

# Per-target conversion library shared by every framework's test wrappers.
add_subdirectory(common)

if(BASELINE_BUILD)
    # Baseline (kind4/kind8): no migrated archives exist; STAGED_LIBRARIES
    # is empty by design. Iterate the on-disk test framework dirs and
    # descend into each whose corresponding std archive built. blas/lapack
    # always; pblas/scalapack only when MPI was found.
    foreach(_lib blas lapack pblas scalapack)
        set(_libdir ${CMAKE_CURRENT_SOURCE_DIR}/${_lib})
        if(TARGET ${_lib} AND EXISTS ${_libdir}/CMakeLists.txt)
            add_subdirectory(${_lib})
        endif()
    endforeach()
else()
    # Non-baseline builds must carry a non-empty LIB_PREFIX; otherwise
    # the gate below would collapse to ``TARGET <lib>`` and accidentally
    # admit the std archives, running each test tree's non-baseline
    # branch against the wrong (standard-precision) target.
    if(LIB_PREFIX STREQUAL "")
        message(FATAL_ERROR
            "tests/CMakeLists.txt: LIB_PREFIX is empty in a "
            "non-baseline build — refusing to descend (would run "
            "migrated test branches against std archives). Check "
            "target_config.cmake and the active TARGET_NAME.")
    endif()
    foreach(_lib ${STAGED_LIBRARIES})
        set(_libdir ${CMAKE_CURRENT_SOURCE_DIR}/${_lib})
        if(TARGET ${LIB_PREFIX}${_lib}
                AND EXISTS ${_libdir}/CMakeLists.txt)
            add_subdirectory(${_lib})
        endif()
    endforeach()

endif()
