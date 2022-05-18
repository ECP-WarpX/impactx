function(find_pyamrex)
    if(ImpactX_pyamrex_src)
        message(STATUS "Compiling local pyAMReX ...")
        message(STATUS "pyAMReX source path: ${ImpactX_pyamrex_src}")
    elseif(ImpactX_pyamrex_internal)
        message(STATUS "Downloading pyAMReX ...")
        message(STATUS "pyAMReX repository: ${ImpactX_pyamrex_repo} (${ImpactX_pyamrex_branch})")
        include(FetchContent)
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
            FetchContent_GetProperties(fetchedpyamrex)

            if(NOT fetchedpyamrex_POPULATED)
                FetchContent_Populate(fetchedpyamrex)
                add_subdirectory(${fetchedpyamrex_SOURCE_DIR} ${fetchedpyamrex_BINARY_DIR})
            endif()

            # advanced fetch options
            mark_as_advanced(FETCHCONTENT_BASE_DIR)
            mark_as_advanced(FETCHCONTENT_FULLY_DISCONNECTED)
            mark_as_advanced(FETCHCONTENT_QUIET)
            mark_as_advanced(FETCHCONTENT_SOURCE_DIR_FETCHEDpyamrex)
            mark_as_advanced(FETCHCONTENT_UPDATES_DISCONNECTED)
            mark_as_advanced(FETCHCONTENT_UPDATES_DISCONNECTED_FETCHEDpyamrex)
        endif()
    else()
        find_package(pyAMReX 21.02 CONFIG REQUIRED)
        message(STATUS "pyAMReX: Found version '${pyamrex_VERSION}'")
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
set(ImpactX_pyamrex_branch "development"
    CACHE STRING
    "Repository branch for ImpactX_pyamrex_repo if(ImpactX_pyamrex_internal)")

find_pyamrex()
