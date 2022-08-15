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

function(impactx_test_set_pythonpath test_name)
    if(WIN32)
        string(REPLACE ";" "\\;" WIN_PYTHONPATH "$ENV{PYTHONPATH}")
        string(REGEX REPLACE "/" "\\\\" WIN_PYTHON_OUTPUT_DIRECTORY ${CMAKE_PYTHON_OUTPUT_DIRECTORY})
        set_tests_properties(${test_name}
            PROPERTIES ENVIRONMENT "PYTHONPATH=${WIN_PYTHON_OUTPUT_DIRECTORY}\;${WIN_PYTHONPATH}"
        )
    else()
        set_tests_properties(${test_name}
            PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_PYTHON_OUTPUT_DIRECTORY}:$ENV{PYTHONPATH}"
        )
    endif()
endfunction()


# FODO Cell ###################################################################
#
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/FODO)
add_test(NAME FODO.run
         COMMAND $<TARGET_FILE:app> ${ImpactX_SOURCE_DIR}/examples/fodo/input_fodo.in
         WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/FODO
)
add_test(NAME FODO.analysis
         COMMAND ${ImpactX_SOURCE_DIR}/examples/fodo/analysis_fodo.py
         WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/FODO
)
add_test(NAME FODO.plot
         COMMAND ${ImpactX_SOURCE_DIR}/examples/fodo/plot_fodo.py --save-png
         WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/FODO
)

set_tests_properties(FODO.analysis PROPERTIES DEPENDS "FODO.run")
set_tests_properties(FODO.plot PROPERTIES DEPENDS "FODO.run")


# MPI-Parallel FODO Cell ######################################################
#
if(ImpactX_MPI)
    file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/FODO.MPI)
    add_test(NAME FODO.MPI.run
             COMMAND ${MPI_TEST_EXE}
                 $<TARGET_FILE:app> ${ImpactX_SOURCE_DIR}/examples/fodo/input_fodo.in
             WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/FODO.MPI
    )
    add_test(NAME FODO.MPI.analysis
             COMMAND ${ImpactX_SOURCE_DIR}/examples/fodo/analysis_fodo.py
             WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/FODO.MPI
    )
    add_test(NAME FODO.MPI.plot
             COMMAND ${ImpactX_SOURCE_DIR}/examples/fodo/plot_fodo.py --save-png
             WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/FODO.MPI
    )

    set_tests_properties(FODO.MPI.run PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")
    set_tests_properties(FODO.MPI.analysis PROPERTIES DEPENDS "FODO.MPI.run")
    set_tests_properties(FODO.MPI.plot PROPERTIES DEPENDS "FODO.MPI.run")
endif()


# Python: FODO Cell ###########################################################
#
if(ImpactX_PYTHON)
    file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/FODO.py)
    add_test(NAME FODO.py.run
             COMMAND ${Python_EXECUTABLE}
                 ${ImpactX_SOURCE_DIR}/examples/fodo/run_fodo.py
             WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/FODO.py
    )
    add_test(NAME FODO.py.analysis
             COMMAND ${ImpactX_SOURCE_DIR}/examples/fodo/analysis_fodo.py
             WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/FODO.py
    )
    add_test(NAME FODO.py.plot
             COMMAND ${ImpactX_SOURCE_DIR}/examples/fodo/plot_fodo.py --save-png
             WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/FODO.py
    )

    set_tests_properties(FODO.py.analysis PROPERTIES DEPENDS "FODO.py.run")
    set_tests_properties(FODO.py.plot PROPERTIES DEPENDS "FODO.py.run")
    impactx_test_set_pythonpath(FODO.py.run)
endif()


# Python: MPI-parallel FODO Cell ##############################################
#
if(ImpactX_PYTHON AND ImpactX_MPI)
    file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/FODO.py.MPI)
    add_test(NAME FODO.py.MPI.run
             COMMAND ${MPI_TEST_EXE} ${Python_EXECUTABLE}
                 ${ImpactX_SOURCE_DIR}/examples/fodo/run_fodo.py
             WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/FODO.py.MPI
    )
    add_test(NAME FODO.py.MPI.analysis
             COMMAND ${ImpactX_SOURCE_DIR}/examples/fodo/analysis_fodo.py
             WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/FODO.py.MPI
    )
    add_test(NAME FODO.py.MPI.plot
             COMMAND ${ImpactX_SOURCE_DIR}/examples/fodo/plot_fodo.py --save-png
             WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/FODO.py.MPI
    )

    set_tests_properties(FODO.py.MPI.run PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")
    set_tests_properties(FODO.py.MPI.analysis PROPERTIES DEPENDS "FODO.py.MPI.run")
    set_tests_properties(FODO.py.MPI.plot PROPERTIES DEPENDS "FODO.py.MPI.run")
    impactx_test_set_pythonpath(FODO.py.MPI.run)
endif()


# Chicane #####################################################################
#
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/chicane)
add_test(NAME chicane.run
         COMMAND $<TARGET_FILE:app> ${ImpactX_SOURCE_DIR}/examples/chicane/input_chicane.in
         WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/chicane
)
add_test(NAME chicane.analysis
         COMMAND ${ImpactX_SOURCE_DIR}/examples/chicane/analysis_chicane.py
         WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/chicane
)
add_test(NAME chicane.plot
         COMMAND ${ImpactX_SOURCE_DIR}/examples/chicane/plot_chicane.py --save-png
         WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/chicane
)

set_tests_properties(chicane.analysis PROPERTIES DEPENDS "chicane.run")
set_tests_properties(chicane.plot PROPERTIES DEPENDS "chicane.run")


# Python: Chicane #####################################################################
#
if(ImpactX_PYTHON)
    file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/chicane.py)
    add_test(NAME chicane.py.run
            COMMAND ${Python_EXECUTABLE}
                ${ImpactX_SOURCE_DIR}/examples/chicane/run_chicane.py
            WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/chicane.py
    )
    add_test(NAME chicane.py.analysis
            COMMAND ${ImpactX_SOURCE_DIR}/examples/chicane/analysis_chicane.py
            WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/chicane.py
    )
    add_test(NAME chicane.py.plot
            COMMAND ${ImpactX_SOURCE_DIR}/examples/chicane/plot_chicane.py --save-png
            WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/chicane.py
    )

    set_tests_properties(chicane.py.analysis PROPERTIES DEPENDS "chicane.py.run")
    set_tests_properties(chicane.py.plot PROPERTIES DEPENDS "chicane.py.run")
    impactx_test_set_pythonpath(chicane.py.run)
endif()


# Constant Focusing Channel ###################################################
#
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/cfchannel)
add_test(NAME cfchannel.run
         COMMAND $<TARGET_FILE:app> ${ImpactX_SOURCE_DIR}/examples/cfchannel/input_cfchannel.in
         WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/cfchannel
)
add_test(NAME cfchannel.analysis
         COMMAND ${ImpactX_SOURCE_DIR}/examples/cfchannel/analysis_cfchannel.py
         WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/cfchannel
)

set_tests_properties(cfchannel.analysis PROPERTIES DEPENDS "cfchannel.run")


# Kurth Distribution Test ###################################################
#
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/kurth)
add_test(NAME kurth.run
         COMMAND $<TARGET_FILE:app> ${ImpactX_SOURCE_DIR}/examples/kurth/input_kurth.in
         WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/kurth
)
add_test(NAME kurth.analysis
         COMMAND ${ImpactX_SOURCE_DIR}/examples/kurth/analysis_kurth.py
         WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/kurth
)

set_tests_properties(kurth.analysis PROPERTIES DEPENDS "kurth.run")


# 6D Gaussian Distribution Test ###################################################
#
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gaussian)
add_test(NAME gaussian.run
         COMMAND $<TARGET_FILE:app> ${ImpactX_SOURCE_DIR}/examples/distgen/input_gaussian.in
         WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gaussian
)
add_test(NAME gaussian.analysis
         COMMAND ${ImpactX_SOURCE_DIR}/examples/distgen/analysis_gaussian.py
         WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gaussian
)

set_tests_properties(gaussian.analysis PROPERTIES DEPENDS "gaussian.run")


# K-V Distribution Test ############################################################
#
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/kvdist)
add_test(NAME kvdist.run
         COMMAND $<TARGET_FILE:app> ${ImpactX_SOURCE_DIR}/examples/distgen/input_kvdist.in
         WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/kvdist
)
add_test(NAME kvdist.analysis
         COMMAND ${ImpactX_SOURCE_DIR}/examples/distgen/analysis_kvdist.py
         WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/kvdist
)

set_tests_properties(kvdist.analysis PROPERTIES DEPENDS "kvdist.run")


# FODO Cell with RF ###################################################################
#
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/FODO_RF)
add_test(NAME FODO_RF.run
         COMMAND $<TARGET_FILE:app> ${ImpactX_SOURCE_DIR}/examples/fodo_rf/input_fodo_rf.in
         WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/FODO_RF
)
add_test(NAME FODO_RF.analysis
         COMMAND ${ImpactX_SOURCE_DIR}/examples/fodo_rf/analysis_fodo_rf.py
         WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/FODO_RF
)

set_tests_properties(FODO_RF.analysis PROPERTIES DEPENDS "FODO_RF.run")


# 4D Kurth Distribution Test ############################################################
#
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/kurth4d)
add_test(NAME kurth4d.run
         COMMAND $<TARGET_FILE:app> ${ImpactX_SOURCE_DIR}/examples/distgen/input_kurth4d.in
         WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/kurth4d
)
add_test(NAME kurth4d.analysis
         COMMAND ${ImpactX_SOURCE_DIR}/examples/distgen/analysis_kurth4d.py
         WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/kurth4d
)

set_tests_properties(kurth4d.analysis PROPERTIES DEPENDS "kurth4d.run")


# Semi-Gaussian Distribution Test ############################################################
#
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/semigaussian)
add_test(NAME semigaussian.run
         COMMAND $<TARGET_FILE:app> ${ImpactX_SOURCE_DIR}/examples/distgen/input_semigaussian.in
         WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/semigaussian
)
add_test(NAME semigaussian.analysis
         COMMAND ${ImpactX_SOURCE_DIR}/examples/distgen/analysis_semigaussian.py
         WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/semigaussian
)

set_tests_properties(semigaussian.analysis PROPERTIES DEPENDS "semigaussian.run")


# Chain of Multipoles Test ############################################################
#
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/multipole)
add_test(NAME multipole.run
         COMMAND $<TARGET_FILE:app> ${ImpactX_SOURCE_DIR}/examples/multipole/input_multipole.in
         WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/multipole
)
add_test(NAME multipole.analysis
         COMMAND ${ImpactX_SOURCE_DIR}/examples/multipole/analysis_multipole.py
         WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/multipole
)

set_tests_properties(multipole.analysis PROPERTIES DEPENDS "multipole.run")


# IOTA Nonlinear Focusing Channel Test ############################################################
#
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/iotalens)
add_test(NAME iotalens.run
         COMMAND $<TARGET_FILE:app> ${ImpactX_SOURCE_DIR}/examples/iota_lens/input_iotalens.in
         WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/iotalens
)
add_test(NAME iotalens.analysis
         COMMAND ${ImpactX_SOURCE_DIR}/examples/iota_lens/analysis_iotalens.py
         WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/iotalens
)

set_tests_properties(iotalens.analysis PROPERTIES DEPENDS "iotalens.run")


# IOTA Linear Lattice Test ############################################################
#
if(ImpactX_MPI)
    file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/iotalattice.MPI)
    add_test(NAME iotalattice.MPI.run
             COMMAND ${MPI_TEST_EXE}
                 $<TARGET_FILE:app> ${ImpactX_SOURCE_DIR}/examples/iota_lattice/input_iotalattice.in
             WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/iotalattice.MPI
    )
    add_test(NAME iotalattice.MPI.analysis
             COMMAND ${ImpactX_SOURCE_DIR}/examples/iota_lattice/analysis_iotalattice.py
             WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/iotalattice.MPI
    )

    set_tests_properties(iotalattice.MPI.analysis PROPERTIES DEPENDS "iotalattice.MPI.run")
endif()


# Python: IOTA Linear Lattice Test ############################################
#
if(ImpactX_PYTHON AND ImpactX_MPI)
    file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/iotalattice.py.MPI)
    add_test(NAME iotalattice.py.MPI.run
             COMMAND ${MPI_TEST_EXE} ${Python_EXECUTABLE}
                 ${ImpactX_SOURCE_DIR}/examples/iota_lattice/run_iotalattice.py
             WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/iotalattice.py.MPI
    )
    add_test(NAME iotalattice.py.MPI.analysis
             COMMAND ${ImpactX_SOURCE_DIR}/examples/iota_lattice/analysis_iotalattice.py
             WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/iotalattice.py.MPI
    )

    set_tests_properties(iotalattice.py.MPI.analysis PROPERTIES DEPENDS "iotalattice.py.MPI.run")
    impactx_test_set_pythonpath(iotalattice.py.MPI.run)
endif()
