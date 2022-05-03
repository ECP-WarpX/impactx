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


# 6D Gaussian Distribution Test ###################################################
#
add_test(NAME gaussian.run
         COMMAND $<TARGET_FILE:app> ${ImpactX_SOURCE_DIR}/examples/distgen/input_gaussian.in
         WORKING_DIRECTORY ${ImpactX_RUNTIME_OUTPUT_DIRECTORY}
)
add_test(NAME gaussian.analysis
         COMMAND ${ImpactX_SOURCE_DIR}/examples/distgen/analysis_gaussian.py
         WORKING_DIRECTORY ${ImpactX_RUNTIME_OUTPUT_DIRECTORY}
)

set_tests_properties(gaussian.analysis PROPERTIES DEPENDS "gaussian.run")


# K-V Distribution Test ############################################################
#
add_test(NAME kvdist.run
         COMMAND $<TARGET_FILE:app> ${ImpactX_SOURCE_DIR}/examples/distgen/input_kvdist.in
         WORKING_DIRECTORY ${ImpactX_RUNTIME_OUTPUT_DIRECTORY}
)
add_test(NAME kvdist.analysis
         COMMAND ${ImpactX_SOURCE_DIR}/examples/distgen/analysis_kvdist.py
         WORKING_DIRECTORY ${ImpactX_RUNTIME_OUTPUT_DIRECTORY}
)

set_tests_properties(kvdist.analysis PROPERTIES DEPENDS "kvdist.run")


# FODO Cell with RF ###################################################################
#
add_test(NAME FODO_RF.run
         COMMAND $<TARGET_FILE:app> ${ImpactX_SOURCE_DIR}/examples/fodo_rf/input_fodo_rf.in
         WORKING_DIRECTORY ${ImpactX_RUNTIME_OUTPUT_DIRECTORY}
)
add_test(NAME FODO_RF.analysis
         COMMAND ${ImpactX_SOURCE_DIR}/examples/fodo_rf/analysis_fodo_rf.py
         WORKING_DIRECTORY ${ImpactX_RUNTIME_OUTPUT_DIRECTORY}
)

set_tests_properties(FODO_RF.analysis PROPERTIES DEPENDS "FODO_RF.run")


# 4D Kurth Distribution Test ############################################################
#
add_test(NAME kurth4d.run
         COMMAND $<TARGET_FILE:app> ${ImpactX_SOURCE_DIR}/examples/distgen/input_kurth4d.in
         WORKING_DIRECTORY ${ImpactX_RUNTIME_OUTPUT_DIRECTORY}
)
add_test(NAME kurth4d.analysis
         COMMAND ${ImpactX_SOURCE_DIR}/examples/distgen/analysis_kurth4d.py
         WORKING_DIRECTORY ${ImpactX_RUNTIME_OUTPUT_DIRECTORY}
)

set_tests_properties(kurth4d.analysis PROPERTIES DEPENDS "kurth4d.run")


# Semi-Gaussian Distribution Test ############################################################
#
add_test(NAME semigaussian.run
         COMMAND $<TARGET_FILE:app> ${ImpactX_SOURCE_DIR}/examples/distgen/input_semigaussian.in
         WORKING_DIRECTORY ${ImpactX_RUNTIME_OUTPUT_DIRECTORY}
)
add_test(NAME semigaussian.analysis
         COMMAND ${ImpactX_SOURCE_DIR}/examples/distgen/analysis_semigaussian.py
         WORKING_DIRECTORY ${ImpactX_RUNTIME_OUTPUT_DIRECTORY}
)

set_tests_properties(semigaussian.analysis PROPERTIES DEPENDS "semigaussian.run")

# Chain of Multipoles Test ############################################################
#
add_test(NAME multipole.run
         COMMAND $<TARGET_FILE:app> ${ImpactX_SOURCE_DIR}/examples/multipole/input_multipole.in
         WORKING_DIRECTORY ${ImpactX_RUNTIME_OUTPUT_DIRECTORY}
)
add_test(NAME multipole.analysis
         COMMAND ${ImpactX_SOURCE_DIR}/examples/multipole/analysis_multipole.py
         WORKING_DIRECTORY ${ImpactX_RUNTIME_OUTPUT_DIRECTORY}
)

set_tests_properties(multipole.analysis PROPERTIES DEPENDS "multipole.run")
