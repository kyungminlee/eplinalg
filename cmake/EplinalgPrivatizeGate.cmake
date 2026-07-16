# eplinalg — ep_ symbol-privatization audit gate (task 44, Phase 3)
#
# The migrator's privatize pass rewrites every symbol in
# ``privatize_ep_symbols.txt`` to its ``ep_`` twin at source level, so the
# extended-precision stack never shares a BLACS/PBLAS/ScaLAPACK-common
# symbol (and therefore never shares state) with an MKL linked into the
# same consumer. This gate is the build-time proof: after each privatized
# STATIC archive links, ``nm`` it and fail the build if any manifest name
# survives — defined (a rename the pass missed) or referenced (a caller
# that would bind to MKL's incompatible engine at consumer link time).
#
# include()'d from the staged top-level CMakeLists.txt. The actual check
# runs in script mode (EplinalgPrivatizeGateCheck.cmake) so it re-executes
# on every relink, not just at configure time.

# Attach the audit to <target>'s POST_BUILD step. No-op if <target> was
# not created in this configuration (mirrors link_if_present's tolerance).
function(eplinalg_privatize_gate target)
    if(NOT TARGET ${target})
        return()
    endif()
    set(_manifest "${CMAKE_SOURCE_DIR}/privatize_ep_symbols.txt")
    if(NOT EXISTS "${_manifest}")
        message(FATAL_ERROR
            "eplinalg_privatize_gate(${target}): manifest not found at "
            "${_manifest}. The staging step copies "
            "codegen/recipes/privatize_ep_symbols.txt next to the staged "
            "CMakeLists.txt; re-run `migrator stage`.")
    endif()
    set(_nm "${CMAKE_NM}")
    if(NOT _nm)
        set(_nm nm)
    endif()
    add_custom_command(TARGET ${target} POST_BUILD
        COMMAND ${CMAKE_COMMAND}
            -DARCHIVE=$<TARGET_FILE:${target}>
            -DMANIFEST=${_manifest}
            -DNM=${_nm}
            -DLABEL=${target}
            -P ${CMAKE_SOURCE_DIR}/EplinalgPrivatizeGateCheck.cmake
        COMMENT "privatize gate: auditing ${target} for pristine BLACS/PBLAS/ScaLAPACK symbols"
        VERBATIM)
endfunction()
