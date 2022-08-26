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


function(add_impactx_test name input is_mpi is_python analysis_script plot_script)
    # make a unique run directory
    file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${name})

    # test run
    set(THIS_WORKING_DIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${name})
    set(THIS_MPI_TEST_EXE)
    if(is_mpi)
        set(THIS_MPI_TEST_EXE ${MPI_TEST_EXE})
    endif()
    set(THIS_Python_EXE)
    if(is_python)
        set(THIS_Python_EXE ${Python_EXECUTABLE})
    endif()
    if(is_python)
        add_test(NAME ${name}.run
                 COMMAND ${THIS_MPI_TEST_EXE} ${THIS_Python_EXE} ${ImpactX_SOURCE_DIR}/${input}
                 WORKING_DIRECTORY ${THIS_WORKING_DIR}
        )
        impactx_test_set_pythonpath(${name}.run)
    else()
        add_test(NAME ${name}.run
                 COMMAND ${THIS_MPI_TEST_EXE} $<TARGET_FILE:app> ${ImpactX_SOURCE_DIR}/${input}
                 WORKING_DIRECTORY ${THIS_WORKING_DIR}
        )
    endif()
    if(is_mpi)
        set_tests_properties(${name}.run PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")
    else()
        # TODO: Change to 2 after OpenMP support was added #195
        set_tests_properties(${name}.run PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")
    endif()

    # analysis and plots
    set(THIS_Python_SCRIPT_EXE)
    if(is_python OR WIN32)
        set(THIS_Python_SCRIPT_EXE ${Python_EXECUTABLE})
    endif()
    if(analysis_script)
        add_test(NAME ${name}.analysis
                 COMMAND ${THIS_Python_SCRIPT_EXE} ${ImpactX_SOURCE_DIR}/${analysis_script}
                 WORKING_DIRECTORY ${THIS_WORKING_DIR}
        )
        set_tests_properties(${name}.analysis PROPERTIES DEPENDS "${name}.run")
    endif()
    if(plot_script)
        add_test(NAME ${name}.plot
                 COMMAND ${THIS_Python_SCRIPT_EXE} ${ImpactX_SOURCE_DIR}/${plot_script}
                 WORKING_DIRECTORY ${THIS_WORKING_DIR}
        )
        set_tests_properties(${name}.plot PROPERTIES DEPENDS "${name}.run")
    endif()
endfunction()


# FODO Cell ###################################################################
#
add_impactx_test(FODO
    examples/fodo/input_fodo.in
      OFF  # ImpactX MPI-parallel
      OFF  # ImpactX Python interface
    examples/fodo/analysis_fodo.py
    examples/fodo/plot_fodo.py --save-png
)

# MPI-Parallel FODO Cell ######################################################
#
add_impactx_test(FODO.MPI
    examples/fodo/input_fodo.in
      ON   # ImpactX MPI-parallel
      OFF  # ImpactX Python interface
    examples/fodo/analysis_fodo.py
    examples/fodo/plot_fodo.py --save-png
)

# Python: FODO Cell ###########################################################
#
add_impactx_test(FODO.py
    examples/fodo/run_fodo.py
      OFF  # ImpactX MPI-parallel
      ON   # ImpactX Python interface
    examples/fodo/analysis_fodo.py
    examples/fodo/plot_fodo.py --save-png
)

# Python: MPI-parallel FODO Cell ##############################################
#
add_impactx_test(FODO.py.MPI
    examples/fodo/run_fodo.py
      ON  # ImpactX MPI-parallel
      ON   # ImpactX Python interface
    examples/fodo/analysis_fodo.py
    examples/fodo/plot_fodo.py --save-png
)

# Chicane #####################################################################
#
add_impactx_test(chicane
    examples/chicane/input_chicane.in
      OFF  # ImpactX MPI-parallel
      OFF  # ImpactX Python interface
    examples/chicane/analysis_chicane.py
    examples/chicane/plot_chicane.py --save-png
)

# Python: Chicane #####################################################################
#
add_impactx_test(chicane.py
    examples/chicane/run_chicane.py
      OFF  # ImpactX MPI-parallel
      ON   # ImpactX Python interface
    examples/chicane/analysis_chicane.py
    examples/chicane/plot_chicane.py --save-png
)

# Constant Focusing Channel ###################################################
#
add_impactx_test(cfchannel
    examples/cfchannel/input_cfchannel.in
      OFF  # ImpactX MPI-parallel
      OFF  # ImpactX Python interface
    examples/cfchannel/analysis_cfchannel.py
    OFF  # not plotting script yet
)

# Kurth Distribution Test ###################################################
#
add_impactx_test(kurth
    examples/kurth/input_kurth.in
      OFF  # ImpactX MPI-parallel
      OFF  # ImpactX Python interface
    examples/kurth/analysis_kurth.py
    OFF  # not plotting script yet
)

# Python: Kurth Distribution Test ##################################################
#
add_impactx_test(kurth.py
    examples/kurth/run_kurth.py
      OFF  # ImpactX MPI-parallel
      ON  # ImpactX Python interface
    examples/kurth/analysis_kurth.py
    OFF  # not plotting script yet
)

# 6D Gaussian Distribution Test ###################################################
#
add_impactx_test(gaussian
    examples/distgen/input_gaussian.in
      OFF  # ImpactX MPI-parallel
      OFF  # ImpactX Python interface
    examples/distgen/analysis_gaussian.py
    OFF  # not plotting script yet
)

# K-V Distribution Test ############################################################
#
add_impactx_test(kvdist
    examples/distgen/input_kvdist.in
      OFF  # ImpactX MPI-parallel
      OFF  # ImpactX Python interface
    examples/distgen/analysis_kvdist.py
    OFF  # not plotting script yet
)

# FODO Cell with RF ###################################################################
#
add_impactx_test(FODO_RF
    examples/fodo_rf/input_fodo_rf.in
      OFF  # ImpactX MPI-parallel
      OFF  # ImpactX Python interface
    examples/fodo_rf/analysis_fodo_rf.py
    OFF  # not plotting script yet
)

# 4D Kurth Distribution Test ############################################################
#
add_impactx_test(kurth4d
    examples/distgen/input_kurth4d.in
      OFF  # ImpactX MPI-parallel
      OFF  # ImpactX Python interface
    examples/distgen/analysis_kurth4d.py
    OFF  # not plotting script yet
)

# Semi-Gaussian Distribution Test ############################################################
#
add_impactx_test(semigaussian
    examples/distgen/input_semigaussian.in
      OFF  # ImpactX MPI-parallel
      OFF  # ImpactX Python interface
    examples/distgen/analysis_semigaussian.py
    OFF  # not plotting script yet
)

# Chain of Multipoles Test ############################################################
#
add_impactx_test(multipole
    examples/multipole/input_multipole.in
      OFF  # ImpactX MPI-parallel
      OFF  # ImpactX Python interface
    examples/multipole/analysis_multipole.py
    OFF  # not plotting script yet
)

# IOTA Nonlinear Focusing Channel Test ############################################################
#
add_impactx_test(iotalens
    examples/iota_lens/input_iotalens.in
      OFF  # ImpactX MPI-parallel
      OFF  # ImpactX Python interface
    examples/iota_lens/analysis_iotalens.py
    OFF  # not plotting script yet
)

# Python: IOTA Nonlinear Focusing Channel Test #####################################
#
add_impactx_test(iotalens.py
    examples/iota_lens/run_iotalens.py
      OFF  # ImpactX MPI-parallel
      ON  # ImpactX Python interface
    examples/iota_lens/analysis_iotalens.py
    OFF  # not plotting script yet
)

# IOTA Linear Lattice Test ############################################################
#
add_impactx_test(iotalattice.MPI
    examples/iota_lattice/input_iotalattice.in
      ON   # ImpactX MPI-parallel
      OFF  # ImpactX Python interface
    examples/iota_lattice/analysis_iotalattice.py
    OFF  # not plotting script yet
)


# Python: IOTA Linear Lattice Test ############################################
#
add_impactx_test(iotalattice.py.MPI
    examples/iota_lattice/run_iotalattice.py
      ON   # ImpactX MPI-parallel
      ON   # ImpactX Python interface
    examples/iota_lattice/analysis_iotalattice.py
    OFF  # not plotting script yet
)
