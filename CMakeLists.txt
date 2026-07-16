# Repo-root project. The migrated library stacks are NOT built from
# here: `migrator stage <dir> --target <t>` materializes a
# self-contained staging tree (whose top-level CMakeLists.txt is
# cmake/CMakeLists.txt) and that tree is what gets configured and
# built — see doc/dev/build.md. This root project only carries the
# version and the Sphinx documentation target.
cmake_minimum_required(VERSION 3.20)

file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/VERSION" EPLINALG_VERSION LIMIT_COUNT 1)
project(eplinalg VERSION ${EPLINALG_VERSION} LANGUAGES NONE)

message(STATUS "eplinalg ${PROJECT_VERSION} — library stacks are built from a "
               "staged tree: `uv run python -m migrator stage <dir> --target <t>` "
               "then `cmake -S <dir> -B <dir>/build` (doc/dev/build.md)")

# ── Documentation (optional) ─────────────────────────────────────────
# `cmake --build build --target doc` renders doc/ with Sphinx when
# sphinx-build is available (pip install -r doc/requirements.txt).
find_program(SPHINX_BUILD sphinx-build)
if(SPHINX_BUILD)
    configure_file(doc/conf.py.in ${CMAKE_BINARY_DIR}/doc/conf.py @ONLY)
    add_custom_target(doc
        COMMAND ${SPHINX_BUILD} -b html
                -c ${CMAKE_BINARY_DIR}/doc
                ${CMAKE_CURRENT_SOURCE_DIR}/doc
                ${CMAKE_BINARY_DIR}/doc/html
        COMMENT "Building HTML documentation with Sphinx")
else()
    message(STATUS "sphinx-build not found — `doc` target disabled")
endif()
