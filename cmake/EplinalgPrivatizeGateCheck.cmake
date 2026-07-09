# eplinalg — privatize-gate check body (task 44, Phase 3)
#
# Script mode (cmake -P). Inputs: ARCHIVE, MANIFEST, NM, LABEL.
# Audits one static archive against the privatization manifest:
#
#   1. No manifest name may be DEFINED — a definition means the migrator's
#      rename pass missed a symbol and the archive would collide with
#      MKL's identically-named engine (the task 44 BLACBUFF 48-vs-56 ABI
#      SIGSEGV).
#   2. No manifest name may be UNDEFINED (referenced) — a leftover
#      reference would bind to MKL's incompatible implementation at
#      consumer link time, defeating the hermetic cut.
#
# Whitelisted references: LAPACK's arithmetic-independent aux routines,
# which the extended stack intentionally resolves from the shared LAPACK
# archive (they carry no BLACS/PBLAS state and exist identically in every
# implementation). They are not manifest names today; the subtraction
# guards a future manifest regeneration under a different decoration.

# Script mode starts with no policies set; IN_LIST needs CMP0057.
cmake_minimum_required(VERSION 3.20)

foreach(_var ARCHIVE MANIFEST NM LABEL)
    if(NOT DEFINED ${_var})
        message(FATAL_ERROR "privatize gate: missing -D${_var}=")
    endif()
endforeach()

set(_whitelist
    lsame_ xerbla_ ilaenv_ ieeeck_ iparmq_
    chla_transtype_ ilaprec_ ilatrans_ ilauplo_)

file(STRINGS "${MANIFEST}" _raw)
set(_names "")
foreach(_line IN LISTS _raw)
    string(STRIP "${_line}" _line)
    if(_line AND NOT _line MATCHES "^#")
        list(APPEND _names "${_line}")
    endif()
endforeach()
list(REMOVE_DUPLICATES _names)
if(NOT _names)
    message(FATAL_ERROR "privatize gate: empty manifest ${MANIFEST}")
endif()

execute_process(
    COMMAND "${NM}" --format=posix "${ARCHIVE}"
    OUTPUT_VARIABLE _nm_out
    ERROR_VARIABLE _nm_err
    RESULT_VARIABLE _nm_rc)
if(NOT _nm_rc EQUAL 0)
    message(FATAL_ERROR
        "privatize gate: `${NM} --format=posix ${ARCHIVE}` failed "
        "(rc=${_nm_rc}): ${_nm_err}")
endif()

# POSIX format: one symbol per line, ``<name> <type> [<value> <size>]``;
# archive-member headers (``lib.a[member.o]:``) contain no space and are
# skipped by the match.
string(REPLACE ";" "" _nm_out "${_nm_out}")
string(REPLACE "\n" ";" _nm_lines "${_nm_out}")
set(_defined "")
set(_undefined "")
foreach(_line IN LISTS _nm_lines)
    if(_line MATCHES "^([^ ]+) ([A-Za-z])( |$)")
        if(CMAKE_MATCH_2 STREQUAL "U")
            list(APPEND _undefined "${CMAKE_MATCH_1}")
        else()
            list(APPEND _defined "${CMAKE_MATCH_1}")
        endif()
    endif()
endforeach()
list(REMOVE_DUPLICATES _defined)
list(REMOVE_DUPLICATES _undefined)

set(_bad_defs "")
set(_bad_refs "")
foreach(_n IN LISTS _names)
    if(_n IN_LIST _defined)
        list(APPEND _bad_defs "${_n}")
    endif()
    if(_n IN_LIST _undefined AND NOT _n IN_LIST _whitelist)
        list(APPEND _bad_refs "${_n}")
    endif()
endforeach()

if(_bad_defs OR _bad_refs)
    set(_msg "privatize gate FAILED for ${LABEL} (${ARCHIVE}):\n")
    if(_bad_defs)
        list(LENGTH _bad_defs _nd)
        list(JOIN _bad_defs "\n    " _defs_str)
        string(APPEND _msg
            "  ${_nd} pristine manifest symbol(s) DEFINED (rename pass "
            "missed them; would collide with MKL):\n    ${_defs_str}\n")
    endif()
    if(_bad_refs)
        list(LENGTH _bad_refs _nr)
        list(JOIN _bad_refs "\n    " _refs_str)
        string(APPEND _msg
            "  ${_nr} pristine manifest symbol(s) REFERENCED (would bind "
            "to MKL's incompatible engine at consumer link):\n    ${_refs_str}\n")
    endif()
    string(APPEND _msg
        "  The archive is NOT hermetic against an MKL-linked consumer. "
        "Re-stage after fixing the migrator's privatize pass, or "
        "regenerate the manifest (recipes/privatize_ep_symbols.txt) if "
        "the symbol surface legitimately changed.")
    message(FATAL_ERROR "${_msg}")
endif()

list(LENGTH _names _n_names)
list(LENGTH _defined _n_defined)
message(STATUS
    "privatize gate: ${LABEL} OK — ${_n_defined} defined symbols audited "
    "against ${_n_names} manifest names, 0 pristine defs, 0 pristine refs")
