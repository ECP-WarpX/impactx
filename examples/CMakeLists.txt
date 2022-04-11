# Configuration ###############################################################
#
if(ImpactX_MPI)
    # OpenMPI root guard: https://github.com/open-mpi/ompi/issues/4451
    if("$ENV{USER}" STREQUAL "root")
        # calling even --help as root will abort and warn on stderr
        execute_process(COMMAND ${MPIEXEC_EXECUTABLE} --help
                ERROR_VARIABLE MPIEXEC_HELP_TEXT
                OUTPUT_STRIP_TRAILING_WHITESPACE)
        if(${MPIEXEC_HELP_TEXT} MATCHES "^.*allow-run-as-root.*$")
            set(MPI_ALLOW_ROOT --allow-run-as-root)
        endif()
    endif()
    set(MPI_TEST_EXE
        ${MPIEXEC_EXECUTABLE}
        ${MPI_ALLOW_ROOT}
        ${MPIEXEC_NUMPROC_FLAG} 2
    )
endif()


# FODO Cell ###################################################################
#
add_test(NAME FODO.run
         COMMAND $<TARGET_FILE:app> ${ImpactX_SOURCE_DIR}/examples/fodo/input_fodo.in
         WORKING_DIRECTORY ${ImpactX_RUNTIME_OUTPUT_DIRECTORY}
)
add_test(NAME FODO.analysis
         COMMAND ${ImpactX_SOURCE_DIR}/examples/fodo/analysis_fodo.py
         WORKING_DIRECTORY ${ImpactX_RUNTIME_OUTPUT_DIRECTORY}
)

set_tests_properties(FODO.analysis PROPERTIES DEPENDS "FODO.run")


# MPI-Parallel FODO Cell ######################################################
#
if(ImpactX_MPI)
    add_test(NAME FODO.MPI.run
            COMMAND ${MPI_TEST_EXE}
                $<TARGET_FILE:app> ${ImpactX_SOURCE_DIR}/examples/fodo/input_fodo.in
            WORKING_DIRECTORY ${ImpactX_RUNTIME_OUTPUT_DIRECTORY}
    )
    add_test(NAME FODO.MPI.analysis
            COMMAND ${ImpactX_SOURCE_DIR}/examples/fodo/analysis_fodo.py
            WORKING_DIRECTORY ${ImpactX_RUNTIME_OUTPUT_DIRECTORY}
    )

    set_tests_properties(FODO.MPI.run PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")
    set_tests_properties(FODO.MPI.analysis PROPERTIES DEPENDS "FODO.MPI.run")
endif()


# Chicane #####################################################################
#
add_test(NAME chicane.run
         COMMAND $<TARGET_FILE:app> ${ImpactX_SOURCE_DIR}/examples/chicane/input_chicane.in
         WORKING_DIRECTORY ${ImpactX_RUNTIME_OUTPUT_DIRECTORY}
)
add_test(NAME chicane.analysis
         COMMAND ${ImpactX_SOURCE_DIR}/examples/chicane/analysis_chicane.py
         WORKING_DIRECTORY ${ImpactX_RUNTIME_OUTPUT_DIRECTORY}
)

set_tests_properties(chicane.analysis PROPERTIES DEPENDS "chicane.run")


# Constant Focusing Channel ###################################################
#
add_test(NAME cfchannel.run
         COMMAND $<TARGET_FILE:app> ${ImpactX_SOURCE_DIR}/examples/cfchannel/input_cfchannel.in
         WORKING_DIRECTORY ${ImpactX_RUNTIME_OUTPUT_DIRECTORY}
)
add_test(NAME cfchannel.analysis
         COMMAND ${ImpactX_SOURCE_DIR}/examples/cfchannel/analysis_cfchannel.py
         WORKING_DIRECTORY ${ImpactX_RUNTIME_OUTPUT_DIRECTORY}
)

set_tests_properties(cfchannel.analysis PROPERTIES DEPENDS "cfchannel.run")


# Kurth Distribution Test ###################################################
#
add_test(NAME kurth.run
         COMMAND $<TARGET_FILE:app> ${ImpactX_SOURCE_DIR}/examples/kurth/input_kurth.in
         WORKING_DIRECTORY ${ImpactX_RUNTIME_OUTPUT_DIRECTORY}
)
add_test(NAME kurth.analysis
         COMMAND ${ImpactX_SOURCE_DIR}/examples/kurth/analysis_kurth.py
         WORKING_DIRECTORY ${ImpactX_RUNTIME_OUTPUT_DIRECTORY}
)

set_tests_properties(kurth.analysis PROPERTIES DEPENDS "kurth.run")
