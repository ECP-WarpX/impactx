macro(find_ablastr)
    if(ImpactX_ablastr_src)
        message(STATUS "Compiling local ABLASTR ...")
        message(STATUS "ABLASTR source path: ${ImpactX_ablastr_src}")
    elseif(ImpactX_ablastr_internal)
        message(STATUS "Downloading ABLASTR ...")
        message(STATUS "ABLASTR repository: ${ImpactX_ablastr_repo} (${ImpactX_ablastr_branch})")
        include(FetchContent)
    endif()

    if(ImpactX_amrex_src)
        set(WarpX_amrex_src ${ImpactX_amrex_src} CACHE PATH
            "Local path to AMReX source directory (preferred if set)" FORCE)
        list(APPEND CMAKE_MODULE_PATH "${WarpX_amrex_src}/Tools/CMake")
    endif()

    if(ImpactX_ablastr_internal OR ImpactX_ablastr_src)
        set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)

        set(WarpX_APP OFF CACHE BOOL "" FORCE)
        set(WarpX_LIB OFF CACHE BOOL "" FORCE)
        set(WarpX_QED OFF CACHE BOOL "" FORCE)
        set(WarpX_DIMS 3 CACHE INTERNAL "" FORCE)
        set(WarpX_COMPUTE ${ImpactX_COMPUTE} CACHE INTERNAL "" FORCE)
        set(WarpX_OPENPMD ${ImpactX_OPENPMD} CACHE INTERNAL "" FORCE)
        set(WarpX_PRECISION ${ImpactX_PRECISION} CACHE INTERNAL "" FORCE)
        set(WarpX_MPI_THREAD_MULTIPLE ${ImpactX_MPI_THREAD_MULTIPLE} CACHE INTERNAL "" FORCE)
        set(WarpX_IPO ${ImpactX_IPO} CACHE INTERNAL "" FORCE)

        if(ImpactX_ablastr_src)
            #list(APPEND CMAKE_MODULE_PATH "${WarpX_amrex_src}/Tools/CMake")
            if(ImpactX_COMPUTE STREQUAL CUDA)
                enable_language(CUDA)
                # ABLASTR 21.06+ supports CUDA_ARCHITECTURES
                #if(CMAKE_VERSION VERSION_LESS 3.20)
                #    include(ABLASTR_SetupCUDA)
                #endif()
            endif()
            add_subdirectory(${ImpactX_ablastr_src} _deps/localablastr-build/)
            # TODO: this is a bit hacky, check if we find a variable like
            #       fetchedamrex_SOURCE_DIR or FETCHCONTENT_SOURCE_DIR_FETCHEDAMREX
            #       that we could use for the named path instead
            list(APPEND CMAKE_MODULE_PATH "${FETCHCONTENT_BASE_DIR}/fetchedamrex-src/Tools/CMake")
        else()
            FetchContent_Declare(fetchedablastr
                GIT_REPOSITORY ${ImpactX_ablastr_repo}
                GIT_TAG        ${ImpactX_ablastr_branch}
                BUILD_IN_SOURCE 0
            )
            FetchContent_GetProperties(fetchedablastr)

            if(NOT fetchedablastr_POPULATED)
                FetchContent_Populate(fetchedablastr)
                #list(APPEND CMAKE_MODULE_PATH "${fetchedamrex_SOURCE_DIR}/Tools/CMake")
                if(ImpactX_COMPUTE STREQUAL CUDA)
                    enable_language(CUDA)
                    # ABLASTR 21.06+ supports CUDA_ARCHITECTURES
                    #if(CMAKE_VERSION VERSION_LESS 3.20)
                    #    include(ABLASTR_SetupCUDA)
                    #endif()
                endif()
                add_subdirectory(${fetchedablastr_SOURCE_DIR} ${fetchedablastr_BINARY_DIR})
                # TODO: this is a bit hacky, check if we find a variable like
                #       fetchedamrex_SOURCE_DIR or FETCHCONTENT_SOURCE_DIR_FETCHEDAMREX
                #       that we could use for the named path instead
                list(APPEND CMAKE_MODULE_PATH "${FETCHCONTENT_BASE_DIR}/fetchedamrex-src/Tools/CMake")
            endif()

            # advanced fetch options
            mark_as_advanced(FETCHCONTENT_BASE_DIR)
            mark_as_advanced(FETCHCONTENT_FULLY_DISCONNECTED)
            mark_as_advanced(FETCHCONTENT_QUIET)
            mark_as_advanced(FETCHCONTENT_SOURCE_DIR_FETCHEDABLASTR)
            mark_as_advanced(FETCHCONTENT_UPDATES_DISCONNECTED)
            mark_as_advanced(FETCHCONTENT_UPDATES_DISCONNECTED_FETCHEDABLASTR)
        endif()

        # ABLASTR options not relevant to most ImpactX users
        mark_as_advanced(AMREX_BUILD_DATETIME)

        message(STATUS "ABLASTR: Using version '${WarpX_VERSION}' (${WarpX_GIT_VERSION})")
    else()
        message(STATUS "Searching for pre-installed ABLASTR ...")
        message(FATAL_ERROR "Not yet supported!")
        set(COMPONENT_DIM 3D)
        set(COMPONENT_PRECISION ${ImpactX_PRECISION} P${ImpactX_PRECISION})

        find_package(ABLASTR 21.12 CONFIG REQUIRED COMPONENTS ${COMPONENT_DIM})
        message(STATUS "ABLASTR: Found version '${ABLASTR_VERSION}'")
    endif()
endmacro()

# local source-tree
set(ImpactX_amrex_src ""
    CACHE PATH
    "Local path to AMReX source directory (preferred if set)")
set(ImpactX_ablastr_src ""
    CACHE PATH
    "Local path to ABLASTR source directory (preferred if set)")

# Git fetcher
set(ImpactX_ablastr_repo "https://github.com/ECP-WarpX/WarpX.git"
    CACHE STRING
    "Repository URI to pull and build ABLASTR from if(ImpactX_ablastr_internal)")
set(ImpactX_ablastr_branch "2a24ba13c214b3e4345529deecf627d6d113004e"
    CACHE STRING
    "Repository branch for ImpactX_ablastr_repo if(ImpactX_ablastr_internal)")

find_ablastr()
