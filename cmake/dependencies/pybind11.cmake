function(find_pybind11)
    if(ImpactX_pybind11_src)
        message(STATUS "Compiling local pybind11 ...")
        message(STATUS "pybind11 source path: ${ImpactX_pybind11_src}")
        if(NOT IS_DIRECTORY ${ImpactX_pybind11_src})
            message(FATAL_ERROR "Specified directory ImpactX_pybind11_src='${ImpactX_pybind11_src}' does not exist!")
        endif()
    elseif(ImpactX_pybind11_internal)
        message(STATUS "Downloading pybind11 ...")
        message(STATUS "pybind11 repository: ${ImpactX_pybind11_repo} (${ImpactX_pybind11_branch})")
        include(FetchContent)
    endif()
    if(ImpactX_pybind11_internal OR ImpactX_pybind11_src)
        set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)

        if(ImpactX_pybind11_src)
            add_subdirectory(${ImpactX_pybind11_src} _deps/localpybind11-build/)
        else()
            FetchContent_Declare(fetchedpybind11
                GIT_REPOSITORY ${ImpactX_pybind11_repo}
                GIT_TAG        ${ImpactX_pybind11_branch}
                BUILD_IN_SOURCE 0
            )
            FetchContent_MakeAvailable(fetchedpybind11)

            # advanced fetch options
            mark_as_advanced(FETCHCONTENT_BASE_DIR)
            mark_as_advanced(FETCHCONTENT_FULLY_DISCONNECTED)
            mark_as_advanced(FETCHCONTENT_QUIET)
            mark_as_advanced(FETCHCONTENT_SOURCE_DIR_FETCHEDpybind11)
            mark_as_advanced(FETCHCONTENT_UPDATES_DISCONNECTED)
            mark_as_advanced(FETCHCONTENT_UPDATES_DISCONNECTED_FETCHEDpybind11)
        endif()
    else()
        find_package(pybind11 2.12.0 CONFIG REQUIRED)
        message(STATUS "pybind11: Found version '${pybind11_VERSION}'")
    endif()
endfunction()

# local source-tree
set(ImpactX_pybind11_src ""
    CACHE PATH
    "Local path to pybind11 source directory (preferred if set)")

# Git fetcher
option(ImpactX_pybind11_internal "Download & build pybind11" ON)
set(ImpactX_pybind11_repo "https://github.com/pybind/pybind11.git"
    CACHE STRING
    "Repository URI to pull and build pybind11 from if(ImpactX_pybind11_internal)")
set(ImpactX_pybind11_branch "v2.12.0"
    CACHE STRING
    "Repository branch for ImpactX_pybind11_repo if(ImpactX_pybind11_internal)")

find_pybind11()
