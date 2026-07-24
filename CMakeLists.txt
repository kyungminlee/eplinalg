# Repo-root project. The migrated library stacks are NOT built from
# here: `migrator stage <dir> --target <t>` materializes a
# self-contained staging tree (whose top-level CMakeLists.txt is
# cmake/CMakeLists.txt) and that tree is what gets configured and
# built — see doc/dev/build.md. This root project only carries the
# version. Documentation is built by ./doc/build.sh (see doc/README.md),
# not through CMake.
cmake_minimum_required(VERSION 3.20)

file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/VERSION" EPLINALG_VERSION LIMIT_COUNT 1)
project(eplinalg VERSION ${EPLINALG_VERSION} LANGUAGES NONE)

message(STATUS "eplinalg ${PROJECT_VERSION} — library stacks are built from a "
               "staged tree: `uv run python -m migrator stage <dir> --target <t>` "
               "then `cmake -S <dir> -B <dir>/build` (doc/dev/build.md)")
