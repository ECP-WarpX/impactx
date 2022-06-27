add_test(NAME pytest
    COMMAND ${Python_EXECUTABLE} -m pytest
        ${ImpactX_SOURCE_DIR}/tests/python
    WORKING_DIRECTORY
        ${CMAKE_PYTHON_OUTPUT_DIRECTORY}
)
