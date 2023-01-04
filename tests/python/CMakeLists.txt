# add pytest tests
#
set(pytest_name pytest.ImpactX)
set(pytest_rundir ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${pytest_name})

#   make a unique run directory
file(MAKE_DIRECTORY ${pytest_rundir})

#   copy input files
file(COPY ${ImpactX_SOURCE_DIR}/examples
     DESTINATION ${pytest_rundir})

#   run
add_test(NAME ${pytest_name}
    COMMAND ${Python_EXECUTABLE} -m pytest -s -vvvv
        ${ImpactX_SOURCE_DIR}/tests/python
    WORKING_DIRECTORY ${pytest_rundir}
)

#   set PYTHONPATH and PATH (for .dll files)
impactx_test_set_pythonpath(${pytest_name})
