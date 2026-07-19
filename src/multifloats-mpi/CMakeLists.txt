# src/multifloats-mpi — first-party MPI bridge for the multifloats
# (real64x2 double-double) arithmetic. Brought in via add_subdirectory
# from the parent build; NOT a standalone project. Relies on the parent
# scope having already:
#   - found MPI (MPI_CXX_FOUND / MPI::MPI_CXX),
#   - created the `multifloats` interface target (the header-only C++ API),
#   - included FortranCompiler.cmake (provides fortran_module_layout).
# These all hold at the point the parent add_subdirectory()s us.

# multifloats_mpi — C++ MPI handle registration + bridge header
if(MPI_CXX_FOUND)
    add_library(multifloats_mpi STATIC
        "${CMAKE_CURRENT_SOURCE_DIR}/multifloats_mpi.cpp")
    target_include_directories(multifloats_mpi PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
    # mpiseq_dtype_tag.h / ep_mpi_op_commute.h live in the sibling
    # mpiseq dir — the one runtime dir the migrator stages for every
    # target (src/multifloats-mpi and src/mpiseq land as runtime/
    # siblings in the staged tree, so the relative path holds there
    # too). PRIVATE: only multifloats_mpi.cpp reads them; the installed
    # multifloats_bridge.h stays self-contained.
    target_include_directories(multifloats_mpi PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../mpiseq)
    # MPI::MPI_CXX is PRIVATE so consumers (mblacs, mlapack PUBLIC
    # chain, mumps test exes) aren't forced to link libmpicxx.so.
    # Each consumer that wants the real-MPI runtime adds MPI::* on
    # its own — keeps the libmpiseq link path (which substitutes
    # mpiseq for ALL of MPI/BLACS/PBLAS/scalapack) clean.
    # multifloats is the header-only C++ double-double core: needed only to
    # compile multifloats_mpi.cpp (symbols baked into libmultifloats_mpi.a),
    # so it is a BUILD_INTERFACE-only dep. Leaking multifloats:: into the
    # installed interface would make find_package(multifloats_mpi) — and every
    # precision package that find_dependency's it — fail, because the
    # FetchContent'd core is not an installed, find_package-able package.
    target_link_libraries(multifloats_mpi
        PRIVATE MPI::MPI_CXX
        PUBLIC  $<BUILD_INTERFACE:multifloats>)
    target_compile_features(multifloats_mpi PRIVATE cxx_std_17)
else()
    # No MPI — still provide the bridge header via an interface target
    add_library(multifloats_mpi INTERFACE)
    target_include_directories(multifloats_mpi INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
endif()

# multifloats_mpi_f — Fortran module exposing the C-side MPI_FLOAT64X2
# / MPI_MM_SUM / etc. handles to Fortran callers (MUMPS calls MPI
# from Fortran directly with these tokens; without the module they
# have no IMPLICIT type at compile time). Linking into multifloats_mpi
# ensures the C++ runtime that defines mf_mpi_*_f globals is pulled
# in transitively whenever a target depends on multifloats_mpi_f.
if(TARGET multifloats_mpi)
    add_library(multifloats_mpi_f STATIC
        "${CMAKE_CURRENT_SOURCE_DIR}/multifloats_mpi_f.f90")
    fortran_module_layout(multifloats_mpi_f)
    target_link_libraries(multifloats_mpi_f PUBLIC multifloats_mpi)
endif()
