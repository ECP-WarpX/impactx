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


function(add_impactx_test name input is_mpi is_python analysis_script plot_script)
    # cannot run Python tests w/o Python support
    if(is_python AND NOT ImpactX_PYTHON)
        return()
    endif()

    # cannot run MPI tests w/o MPI support
    if(is_mpi AND NOT ImpactX_MPI)
        return()
    endif()

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
                 COMMAND
                    ${THIS_MPI_TEST_EXE} $<TARGET_FILE:app> ${ImpactX_SOURCE_DIR}/${input}
                        amrex.abort_on_unused_inputs=1
                        impactx.always_warn_immediately=1
                        impactx.abort_on_warning_threshold=low
                 WORKING_DIRECTORY ${THIS_WORKING_DIR}
        )
    endif()
    if(is_mpi)
        set_property(TEST ${name}.run APPEND PROPERTY ENVIRONMENT "OMP_NUM_THREADS=1")
    else()
        set_property(TEST ${name}.run APPEND PROPERTY ENVIRONMENT "OMP_NUM_THREADS=2")
    endif()

    # make HDF5 I/O more robust on various filesystems
    set_property(TEST ${name}.run APPEND PROPERTY ENVIRONMENT "HDF5_USE_FILE_LOCKING=FALSE")

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
        set_property(TEST ${name}.analysis APPEND PROPERTY DEPENDS "${name}.run")
    endif()
    if(plot_script)
        add_test(NAME ${name}.plot
                 COMMAND ${THIS_Python_SCRIPT_EXE} ${ImpactX_SOURCE_DIR}/${plot_script} --save-png
                 WORKING_DIRECTORY ${THIS_WORKING_DIR}
        )
        set_property(TEST ${name}.plot APPEND PROPERTY DEPENDS "${name}.run")
    endif()
endfunction()


# FODO Cell ###################################################################
#
add_impactx_test(FODO
    examples/fodo/input_fodo.in
      OFF  # ImpactX MPI-parallel
      OFF  # ImpactX Python interface
    examples/fodo/analysis_fodo.py
    examples/fodo/plot_fodo.py
)

# MPI-Parallel FODO Cell ######################################################
#
add_impactx_test(FODO.MPI
    examples/fodo/input_fodo.in
      ON   # ImpactX MPI-parallel
      OFF  # ImpactX Python interface
    examples/fodo/analysis_fodo.py
    examples/fodo/plot_fodo.py
)

# Python: FODO Cell ###########################################################
#
add_impactx_test(FODO.py
    examples/fodo/run_fodo.py
      OFF  # ImpactX MPI-parallel
      ON   # ImpactX Python interface
    examples/fodo/analysis_fodo.py
    examples/fodo/plot_fodo.py
)

# Python: FODO Cell w/ programmed element for the Drifts ######################
#
add_impactx_test(FODO.programmable.py
    examples/fodo/run_fodo_programmable.py
      OFF  # ImpactX MPI-parallel
      ON   # ImpactX Python interface
    examples/fodo/analysis_fodo.py
    examples/fodo/plot_fodo.py
)

# Python MADX: FODO Cell ######################################################
#
# copy MAD-X lattice file
file(COPY ${ImpactX_SOURCE_DIR}/examples/fodo/fodo.madx
     DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/FODO.MADX.py)

add_impactx_test(FODO.MADX.py
    examples/fodo/run_fodo_madx.py
      OFF  # ImpactX MPI-parallel
      ON   # ImpactX Python interface
    examples/fodo/analysis_fodo.py
    examples/fodo/plot_fodo.py
)

# Python: MPI-parallel FODO Cell ##############################################
#
add_impactx_test(FODO.py.MPI
    examples/fodo/run_fodo.py
      ON  # ImpactX MPI-parallel
      ON   # ImpactX Python interface
    examples/fodo/analysis_fodo.py
    examples/fodo/plot_fodo.py
)

# Chicane #####################################################################
#
add_impactx_test(chicane
    examples/chicane/input_chicane.in
      OFF  # ImpactX MPI-parallel
      OFF  # ImpactX Python interface
    examples/chicane/analysis_chicane.py
    examples/chicane/plot_chicane.py
)

# Python: Chicane #############################################################
#
add_impactx_test(chicane.py
    examples/chicane/run_chicane.py
      OFF  # ImpactX MPI-parallel
      ON   # ImpactX Python interface
    examples/chicane/analysis_chicane.py
    examples/chicane/plot_chicane.py
)

# Python MADX: Chicane ########################################################
#
# copy MAD-X lattice file
file(COPY ${ImpactX_SOURCE_DIR}/examples/chicane/chicane.madx DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/chicane.MADX.py)

add_impactx_test(chicane.MADX.py
    examples/chicane/run_chicane_madx.py
      OFF  # ImpactX MPI-parallel
      ON   # ImpactX Python interface
    examples/chicane/analysis_chicane.py
    examples/chicane/plot_chicane.py
)

# Constant Focusing Channel ###################################################
#
add_impactx_test(cfchannel
    examples/cfchannel/input_cfchannel.in
      OFF  # ImpactX MPI-parallel
      OFF  # ImpactX Python interface
    examples/cfchannel/analysis_cfchannel.py
    OFF  # no plot script yet
)

# Constant Focusing Channel with Space Charge #################################
#
add_impactx_test(cfchannel_spacecharge
    examples/cfchannel/input_cfchannel_10nC.in
      OFF  # ImpactX MPI-parallel
      OFF  # ImpactX Python interface
    examples/cfchannel/analysis_cfchannel_10nC.py
    OFF  # no plot script yet
)

# Python: Constant Focusing Channel with Space Charge #########################
#
add_impactx_test(cfchannel_spacecharge.py
    examples/cfchannel/run_cfchannel_10nC.py
      OFF  # ImpactX MPI-parallel
      ON   # ImpactX Python interface
    examples/cfchannel/analysis_cfchannel_10nC.py
    OFF  # no plot script yet
)

# Python: Constant Focusing Channel ###########################################
#
add_impactx_test(cfchannel.py
    examples/cfchannel/run_cfchannel.py
      OFF  # ImpactX MPI-parallel
      ON   # ImpactX Python interface
    examples/cfchannel/analysis_cfchannel.py
    OFF  # no plot script yet
)

# 6D Gaussian Distribution Test ###############################################
#
add_impactx_test(gaussian
    examples/distgen/input_gaussian.in
      OFF  # ImpactX MPI-parallel
      OFF  # ImpactX Python interface
    examples/distgen/analysis_gaussian.py
    OFF  # no plot script yet
)

# K-V Distribution Test #######################################################
#
add_impactx_test(kvdist
    examples/distgen/input_kvdist.in
      OFF  # ImpactX MPI-parallel
      OFF  # ImpactX Python interface
    examples/distgen/analysis_kvdist.py
    OFF  # no plot script yet
)

# FODO Cell with RF ###########################################################
#
add_impactx_test(FODO_RF
    examples/fodo_rf/input_fodo_rf.in
      OFF  # ImpactX MPI-parallel
      OFF  # ImpactX Python interface
    examples/fodo_rf/analysis_fodo_rf.py
    OFF  # no plot script yet
)

# Python: FODO Cell with RF ###################################################
#
add_impactx_test(FODO_RF.py
    examples/fodo_rf/run_fodo_rf.py
      OFF  # ImpactX MPI-parallel
      ON  # ImpactX Python interface
    examples/fodo_rf/analysis_fodo_rf.py
    OFF  # no plot script yet
)

# 4D Kurth Distribution Test ##################################################
#
add_impactx_test(kurth4d
    examples/distgen/input_kurth4d.in
      OFF  # ImpactX MPI-parallel
      OFF  # ImpactX Python interface
    examples/distgen/analysis_kurth4d.py
    OFF  # no plot script yet
)

# Semi-Gaussian Distribution Test #############################################
#
add_impactx_test(semigaussian
    examples/distgen/input_semigaussian.in
      OFF  # ImpactX MPI-parallel
      OFF  # ImpactX Python interface
    examples/distgen/analysis_semigaussian.py
    OFF  # no plot script yet
)

# Chain of Multipoles Test ####################################################
#
add_impactx_test(multipole
    examples/multipole/input_multipole.in
      OFF  # ImpactX MPI-parallel
      OFF  # ImpactX Python interface
    examples/multipole/analysis_multipole.py
    OFF  # no plot script yet
)

# Expanding Beam Test #########################################################
#
add_impactx_test(expanding_beam
    examples/expanding_beam/input_expanding.in
      OFF  # ImpactX MPI-parallel
      OFF  # ImpactX Python interface
    examples/expanding_beam/analysis_expanding.py
    OFF  # no plot script yet
)

# Python: Expanding Beam Test #################################################
#
add_impactx_test(expanding_beam.py
    examples/expanding_beam/run_expanding.py
      OFF  # ImpactX MPI-parallel
      ON   # ImpactX Python interface
    examples/expanding_beam/analysis_expanding.py
    OFF  # no plot script yet
)

# Python: Chain of Multipoles Test ############################################
#
add_impactx_test(multipole.py
    examples/multipole/run_multipole.py
      OFF  # ImpactX MPI-parallel
      ON  # ImpactX Python interface
    examples/multipole/analysis_multipole.py
    OFF  # not plotting script yet
)

# IOTA Nonlinear Focusing Channel Test ########################################
#
add_impactx_test(iotalens
    examples/iota_lens/input_iotalens.in
      OFF  # ImpactX MPI-parallel
      OFF  # ImpactX Python interface
    examples/iota_lens/analysis_iotalens.py
    OFF  # no plot script yet
)

# Python: IOTA Nonlinear Focusing Channel Test ################################
#
add_impactx_test(iotalens.py
    examples/iota_lens/run_iotalens.py
      OFF  # ImpactX MPI-parallel
      ON  # ImpactX Python interface
    examples/iota_lens/analysis_iotalens.py
    OFF  # not plotting script yet
)

# IOTA Linear Lattice Test ####################################################
#
add_impactx_test(iotalattice.MPI
    examples/iota_lattice/input_iotalattice.in
      ON   # ImpactX MPI-parallel
      OFF  # ImpactX Python interface
    examples/iota_lattice/analysis_iotalattice.py
    OFF  # no plot script yet
)

# Python: IOTA Linear Lattice Test ############################################
#
add_impactx_test(iotalattice.py.MPI
    examples/iota_lattice/run_iotalattice.py
      ON   # ImpactX MPI-parallel
      ON   # ImpactX Python interface
    examples/iota_lattice/analysis_iotalattice.py
    OFF  # no plot script yet
)

# Kurth Beam in a Periodic Channel Test ###################################################
#
# w/o space charge
add_impactx_test(kurth_periodic
    examples/kurth/input_kurth_periodic.in
      OFF  # ImpactX MPI-parallel
      OFF  # ImpactX Python interface
    examples/kurth/analysis_kurth_periodic.py
    OFF  # no plot script yet
)
add_impactx_test(kurth_periodic.py
    examples/kurth/run_kurth_periodic.py
      OFF  # ImpactX MPI-parallel
      ON   # ImpactX Python interface
    examples/kurth/analysis_kurth_periodic.py
    OFF  # no plot script yet
)

# w/ space charge
add_impactx_test(kurth_10nC_periodic
    examples/kurth/input_kurth_10nC_periodic.in
      OFF  # ImpactX MPI-parallel
      OFF  # ImpactX Python interface
    examples/kurth/analysis_kurth_10nC_periodic.py
    OFF  # no plot script yet
)
add_impactx_test(kurth_10nC_periodic.py
    examples/kurth/run_kurth_10nC_periodic.py
      OFF  # ImpactX MPI-parallel
      ON   # ImpactX Python interface
    examples/kurth/analysis_kurth_10nC_periodic.py
    OFF  # no plot script yet
)

# Acceleration by RF Cavities #################################################
#
add_impactx_test(rfcavity.MPI
    examples/rfcavity/input_rfcavity.in
      ON   # ImpactX MPI-parallel
      OFF  # ImpactX Python interface
    examples/rfcavity/analysis_rfcavity.py
    OFF  # no plot script yet
)

# Python: Acceleration by RF Cavities #########################################
#
add_impactx_test(rfcavity.py.MPI
    examples/rfcavity/run_rfcavity.py
      ON   # ImpactX MPI-parallel
      ON   # ImpactX Python interface
    examples/rfcavity/analysis_rfcavity.py
    OFF  # no plot script yet
)

# Ideal, Hard-Edge Solenoid ###################################################
#
# w/o space charge
add_impactx_test(solenoid
    examples/solenoid/input_solenoid.in
      OFF  # ImpactX MPI-parallel
      OFF  # ImpactX Python interface
    examples/solenoid/analysis_solenoid.py
    OFF  # no plot script yet
)
add_impactx_test(solenoid.py
    examples/solenoid/run_solenoid.py
      OFF  # ImpactX MPI-parallel
      ON   # ImpactX Python interface
    examples/solenoid/analysis_solenoid.py
    OFF  # no plot script yet
)

# copy MAD-X lattice file
file(COPY ${ImpactX_SOURCE_DIR}/examples/solenoid/solenoid.madx
     DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/solenoid.MADX.py)
add_impactx_test(solenoid.MADX.py
    examples/solenoid/run_solenoid_madx.py
      OFF  # ImpactX MPI-parallel
      ON   # ImpactX Python interface
    examples/solenoid/analysis_solenoid.py
    OFF  # no plot script yet
)

# Soft-Edge Solenoid ###################################################
#
# w/o space charge
add_impactx_test(solenoid_softedge
    examples/solenoid_softedge/input_solenoid_softedge.in
      OFF  # ImpactX MPI-parallel
      OFF  # ImpactX Python interface
    examples/solenoid_softedge/analysis_solenoid_softedge.py
    OFF  # no plot script yet
)
add_impactx_test(solenoid_softedge.py
    examples/solenoid_softedge/run_solenoid_softedge.py
      OFF  # ImpactX MPI-parallel
      ON   # ImpactX Python interface
    examples/solenoid_softedge/analysis_solenoid_softedge.py
    OFF  # no plot script yet
)

# Soft-Edge Quadrupole ###################################################
#
# w/o space charge
add_impactx_test(quadrupole_softedge
    examples/quadrupole_softedge/input_quadrupole_softedge.in
      OFF  # ImpactX MPI-parallel
      OFF  # ImpactX Python interface
    examples/quadrupole_softedge/analysis_quadrupole_softedge.py
    OFF  # no plot script yet
)
add_impactx_test(quadrupole_softedge.py
    examples/quadrupole_softedge/run_quadrupole_softedge.py
      OFF  # ImpactX MPI-parallel
      ON   # ImpactX Python interface
    examples/quadrupole_softedge/analysis_quadrupole_softedge.py
    OFF  # no plot script yet
)
