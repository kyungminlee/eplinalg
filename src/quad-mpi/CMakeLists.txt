# src/quad-mpi — first-party MPI bridge for the kind16
# (REAL(KIND=16) / __float128) arithmetic. Brought in via
# add_subdirectory from the parent build; NOT a standalone project.
#
# Custom MPI reduce ops (MPI_QQ_SUM / QQ_AMX / QQ_AMN and the XX_*
# complex analogues) registered on the *standard* MPI_REAL16 /
# MPI_COMPLEX32 datatypes. Intel MPI's builtin MPI_SUM/MAX/MIN have no
# 16-byte-real kernel and segfault on MPI_REDUCE at np>=2, so the
# migrated kind16 MUMPS routes every reduction through these instead —
# the quad analogue of the multifloats bridge, but plain C (no C++, no
# FetchContent, no custom datatypes).

# quad_mpi — C registration of the custom combine ops. fabsq lives
# in libquadmath.
add_library(quad_mpi STATIC "${CMAKE_CURRENT_SOURCE_DIR}/quad_mpi.c")
target_include_directories(quad_mpi PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
# ep_mpi_op_commute.h lives in the sibling mpiseq dir — the one runtime
# dir the migrator stages for every target (src/quad-mpi and src/mpiseq
# land as runtime/ siblings in the staged tree, so the relative path
# holds there too). PRIVATE: only quad_mpi.c reads it.
target_include_directories(quad_mpi PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../mpiseq)
# MPI::MPI_C is PRIVATE so consumers linking libmpiseq instead of
# the real MPI runtime (the _seq test variants) aren't forced to
# pull in the real MPI .so — mirrors multifloats_mpi's handling.
# (MPI_Op_c2f is a cast macro and MPI_Op_create is stubbed by
# libmpiseq, so quad_mpi links cleanly against mpiseq too.)
if(MPI_C_FOUND)
    target_link_libraries(quad_mpi PRIVATE MPI::MPI_C)
endif()
target_link_libraries(quad_mpi PRIVATE quadmath)

# quad_mpi_f — Fortran shim exposing the C-side op handles to MUMPS
# (which issues MPI_REDUCE from Fortran with these tokens directly).
# Linking quad_mpi PUBLIC pulls the C runtime that defines the *_f
# globals into any target that depends on quad_mpi_f.
if(TARGET quad_mpi)
    add_library(quad_mpi_f STATIC "${CMAKE_CURRENT_SOURCE_DIR}/quad_mpi_f.f90")
    fortran_module_layout(quad_mpi_f)
    target_link_libraries(quad_mpi_f PUBLIC quad_mpi)
endif()
