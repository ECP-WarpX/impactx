macro(find_ablastr)
    # if pyAMReX is external, AMReX must be as well
    if(DEFINED ImpactX_pyamrex_internal AND NOT ImpactX_pyamrex_internal)
        set(ImpactX_amrex_internal OFF CACHE BOOL
            "Download & build AMReX" FORCE)
    endif()

    if(ImpactX_ablastr_src)
        message(STATUS "Compiling local ABLASTR ...")
        message(STATUS "ABLASTR source path: ${ImpactX_ablastr_src}")
        if(NOT IS_DIRECTORY ${ImpactX_ablastr_src})
            message(FATAL_ERROR "Specified directory ImpactX_ablastr_src='${ImpactX_ablastr_src}' does not exist!")
        endif()
    elseif(ImpactX_ablastr_internal)
        message(STATUS "Downloading ABLASTR ...")
        message(STATUS "ABLASTR repository: ${ImpactX_ablastr_repo} (${ImpactX_ablastr_branch})")
        include(FetchContent)
    endif()

    # transitive control for AMReX superbuild
    set(WarpX_amrex_internal ${ImpactX_amrex_internal} CACHE BOOL
        "Download & build AMReX" FORCE)
    if(ImpactX_amrex_src)
        set(WarpX_amrex_src ${ImpactX_amrex_src} CACHE PATH
            "Local path to AMReX source directory (preferred if set)" FORCE)
        list(APPEND CMAKE_MODULE_PATH "${WarpX_amrex_src}/Tools/CMake")
    elseif(ImpactX_amrex_internal)
        if(ImpactX_amrex_repo)
            set(WarpX_amrex_repo ${ImpactX_amrex_repo} CACHE STRING
                "Repository URI to pull and build AMReX from if(ImpactX_amrex_internal)" FORCE)
        endif()
        if(ImpactX_amrex_branch)
            set(WarpX_amrex_branch ${ImpactX_amrex_branch} CACHE STRING
                "Repository branch for ImpactX_amrex_repo if(ImpactX_amrex_internal)" FORCE)
        endif()
    endif()

    # transitive control for openPMD superbuild
    if(ImpactX_openpmd_src)
        set(WarpX_openpmd_src ${ImpactX_openpmd_src} CACHE PATH
            "Local path to openPMD-api source directory (preferred if set)" FORCE)
    elseif(ImpactX_openpmd_internal)
        if(ImpactX_openpmd_repo)
            set(WarpX_openpmd_repo ${ImpactX_openpmd_repo} CACHE STRING
                "Repository URI to pull and build openPMD-api from if(ImpactX_openpmd_internal)" FORCE)
        endif()
        if(ImpactX_openpmd_branch)
            set(WarpX_openpmd_branch ${ImpactX_openpmd_branch} CACHE STRING
                "Repository branch for ImpactX_openpmd_repo if(ImpactX_openpmd_internal)" FORCE)
        endif()
    else()
        set(WarpX_openpmd_internal ${ImpactX_openpmd_internal} CACHE STRING
            "Download & build openPMD-api" FORCE)
    endif()

    # transitive control for ABLASTR superbuild
    if(ImpactX_ablastr_internal OR ImpactX_ablastr_src)
        set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)

        set(ABLASTR_FFT ${ImpactX_FFT} CACHE BOOL "" FORCE)

        set(WarpX_APP OFF CACHE BOOL "" FORCE)
        set(WarpX_LIB OFF CACHE BOOL "" FORCE)
        set(WarpX_QED OFF CACHE BOOL "" FORCE)
        set(WarpX_COMPUTE ${ImpactX_COMPUTE} CACHE INTERNAL "" FORCE)
        set(WarpX_DIMS 3 CACHE INTERNAL "" FORCE)
        set(WarpX_FFT ${ImpactX_FFT} CACHE BOOL "" FORCE)
        set(WarpX_OPENPMD ${ImpactX_OPENPMD} CACHE INTERNAL "" FORCE)
        set(WarpX_PRECISION ${ImpactX_PRECISION} CACHE INTERNAL "" FORCE)
        set(WarpX_MPI ${ImpactX_MPI} CACHE INTERNAL "" FORCE)
        set(WarpX_MPI_THREAD_MULTIPLE ${ImpactX_MPI_THREAD_MULTIPLE} CACHE INTERNAL "" FORCE)
        set(WarpX_IPO ${ImpactX_IPO} CACHE INTERNAL "" FORCE)

        # shared libs, i.e. for Python bindings, need relocatable code
        if(ImpactX_PYTHON OR BUILD_SHARED_LIBS)
            set(AMReX_PIC ON CACHE INTERNAL
                "Build AMReX with position independent code")
            set(ABLASTR_POSITION_INDEPENDENT_CODE ON CACHE INTERNAL
                "Build ABLASTR with position independent code")

            # WE NEED AMReX AS SHARED LIB, OTHERWISE WE CANNOT SHARE ITS GLOBALS
            # BETWEEN MULTIPLE PYTHON MODULES
            # TODO this is likely an export/symbol hiding issue that we could
            #      alleviate later on
            set(AMReX_BUILD_SHARED_LIBS ON CACHE BOOL "Build AMReX shared library" FORCE)
        endif()

        if(ImpactX_ablastr_src)
            #list(APPEND CMAKE_MODULE_PATH "${WarpX_amrex_src}/Tools/CMake")
            if(ImpactX_COMPUTE STREQUAL CUDA)
                enable_language(CUDA)
            endif()
            add_subdirectory(${ImpactX_ablastr_src} _deps/localablastr-build/)
            # TODO: this is a bit hacky, check if we find a variable like
            #       fetchedamrex_SOURCE_DIR or FETCHCONTENT_SOURCE_DIR_FETCHEDAMREX
            #       or AMReX_DIR or AMReX_MODULE_PATH that we could use for the named path instead
            list(APPEND CMAKE_MODULE_PATH "${FETCHCONTENT_BASE_DIR}/fetchedamrex-src/Tools/CMake")
        else()
            if(ImpactX_COMPUTE STREQUAL CUDA)
                enable_language(CUDA)
            endif()
            FetchContent_Declare(fetchedablastr
                GIT_REPOSITORY ${ImpactX_ablastr_repo}
                GIT_TAG        ${ImpactX_ablastr_branch}
                BUILD_IN_SOURCE 0
            )
            FetchContent_MakeAvailable(fetchedablastr)
            list(APPEND CMAKE_MODULE_PATH "${FETCHCONTENT_BASE_DIR}/fetchedamrex-src/Tools/CMake")

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
        # TODO: MPI & FFT control
        set(COMPONENT_DIM 3D)
        set(COMPONENT_PRECISION ${ImpactX_PRECISION} P${ImpactX_PRECISION})

        find_package(ABLASTR 24.08 CONFIG REQUIRED COMPONENTS ${COMPONENT_DIM})
        message(STATUS "ABLASTR: Found version '${ABLASTR_VERSION}'")
    endif()

    # AMReX CMake helper scripts
    list(APPEND CMAKE_MODULE_PATH "${AMReX_DIR}/AMReXCMakeModules")

    # transitive control for openPMD external
    if(NOT ImpactX_openpmd_src AND NOT ImpactX_openpmd_internal)
        if(ImpactX_MPI)
            set(COMPONENT_WMPI MPI)
        else()
            set(COMPONENT_WMPI NOMPI)
        endif()
        find_package(openPMD 0.15.2 CONFIG REQUIRED COMPONENTS ${COMPONENT_WMPI})
        message(STATUS "openPMD-api: Found version '${openPMD_VERSION}'")
    endif()
endmacro()

# local source-tree
set(ImpactX_amrex_src ""
    CACHE PATH
    "Local path to AMReX source directory (preferred if set)")
set(ImpactX_ablastr_src ""
    CACHE PATH
    "Local path to ABLASTR source directory (preferred if set)")
set(ImpactX_openpmd_src ""
    CACHE PATH
    "Local path to openPMD-api source directory (preferred if set)")

# Git fetcher
set(ImpactX_ablastr_repo "https://github.com/ECP-WarpX/WarpX.git"
    CACHE STRING
    "Repository URI to pull and build ABLASTR from if(ImpactX_ablastr_internal)")
set(ImpactX_ablastr_branch "0838941a6693df711769a90a8a989c9a920b0abe"
    CACHE STRING
    "Repository branch for ImpactX_ablastr_repo if(ImpactX_ablastr_internal)")

# AMReX is transitively pulled through ABLASTR
set(ImpactX_amrex_repo "https://github.com/AMReX-Codes/amrex.git"
    CACHE STRING
    "Repository URI to pull and build AMReX from if(ImpactX_amrex_internal)")
set(ImpactX_amrex_branch ""
    CACHE STRING
    "Repository branch for ImpactX_amrex_repo if(ImpactX_amrex_internal)")

# openPMD is transitively pulled through ABLASTR
set(ImpactX_openpmd_repo "https://github.com/openPMD/openPMD-api.git"
    CACHE STRING
    "Repository URI to pull and build openPMD-api from if(ImpactX_openpmd_internal)")
set(ImpactX_openpmd_branch "0.15.2"
    CACHE STRING
    "Repository branch for ImpactX_openPMD_repo if(ImpactX_openpmd_internal)")

find_ablastr()
