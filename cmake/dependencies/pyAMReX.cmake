function(find_pyamrex)
    if(ImpactX_pyamrex_src)
        message(STATUS "Compiling local pyAMReX ...")
        message(STATUS "pyAMReX source path: ${ImpactX_pyamrex_src}")
        if(NOT IS_DIRECTORY ${ImpactX_pyamrex_src})
            message(FATAL_ERROR "Specified directory ImpactX_pyamrex_src='${ImpactX_pyamrex_src}' does not exist!")
        endif()
    elseif(ImpactX_pyamrex_internal)
        message(STATUS "Downloading pyAMReX ...")
        message(STATUS "pyAMReX repository: ${ImpactX_pyamrex_repo} (${ImpactX_pyamrex_branch})")
        include(FetchContent)
    endif()

    # transitive control for AMReX & pybind11 superbuild
    #   note: if we do superbuilds, we want the same AMReX commit for
    #           AMReX->ABLASTR->ImpactX and
    #           AMReX->pyAMReX->pyImpactX
    #   note: this is performed after we did the transitive logic control in
    #         ABLASTR.cmake
    set(pyAMReX_amrex_internal ${ImpactX_amrex_internal} CACHE BOOL
        "Download & build AMReX" FORCE)
    set(pyAMReX_pybind11_internal ${ImpactX_pybind11_internal} CACHE BOOL
        "Download & build AMReX" FORCE)

    if(ImpactX_amrex_src)
        set(pyAMReX_amrex_src ${ImpactX_amrex_src} CACHE PATH
            "Local path to AMReX source directory (preferred if set)" FORCE)
    elseif(ImpactX_amrex_internal)
        if(ImpactX_amrex_repo)
            set(pyAMReX_amrex_repo ${ImpactX_amrex_repo} CACHE STRING
                "Repository URI to pull and build AMReX from if(ImpactX_amrex_internal)" FORCE)
        endif()
        if(ImpactX_amrex_branch)
            set(pyAMReX_amrex_branch ${ImpactX_amrex_branch} CACHE STRING
                "Repository branch for ImpactX_amrex_repo if(ImpactX_amrex_internal)" FORCE)
        endif()
    endif()

    if(ImpactX_pyamrex_internal OR ImpactX_pyamrex_src)
        set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)

        if(ImpactX_pyamrex_src)
            add_subdirectory(${ImpactX_pyamrex_src} _deps/localpyamrex-build/)
        else()
            FetchContent_Declare(fetchedpyamrex
                GIT_REPOSITORY ${ImpactX_pyamrex_repo}
                GIT_TAG        ${ImpactX_pyamrex_branch}
                BUILD_IN_SOURCE 0
            )
            FetchContent_MakeAvailable(fetchedpyamrex)

            # advanced fetch options
            mark_as_advanced(FETCHCONTENT_BASE_DIR)
            mark_as_advanced(FETCHCONTENT_FULLY_DISCONNECTED)
            mark_as_advanced(FETCHCONTENT_QUIET)
            mark_as_advanced(FETCHCONTENT_SOURCE_DIR_FETCHEDpyamrex)
            mark_as_advanced(FETCHCONTENT_UPDATES_DISCONNECTED)
            mark_as_advanced(FETCHCONTENT_UPDATES_DISCONNECTED_FETCHEDpyamrex)
        endif()
    elseif(NOT ImpactX_pyamrex_internal)
        # TODO: MPI control
        find_package(pyAMReX 24.08 CONFIG REQUIRED)
        message(STATUS "pyAMReX: Found version '${pyAMReX_VERSION}'")
    endif()
endfunction()

# local source-tree
set(ImpactX_pyamrex_src ""
    CACHE PATH
    "Local path to pyAMReX source directory (preferred if set)")

# Git fetcher
option(ImpactX_pyamrex_internal "Download & build pyAMReX" ON)
set(ImpactX_pyamrex_repo "https://github.com/AMReX-Codes/pyamrex.git"
    CACHE STRING
    "Repository URI to pull and build pyamrex from if(ImpactX_pyamrex_internal)")
set(ImpactX_pyamrex_branch "abdf332e25bfeef2b4d613d7adbe93fb8cf3e2f7"
    CACHE STRING
    "Repository branch for ImpactX_pyamrex_repo if(ImpactX_pyamrex_internal)")

find_pyamrex()
