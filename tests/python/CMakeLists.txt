# copy input files
file(COPY ${ImpactX_SOURCE_DIR}/examples DESTINATION ${CMAKE_PYTHON_OUTPUT_DIRECTORY})

# add pytest tests
add_test(NAME pytest.ImpactX
    COMMAND ${Python_EXECUTABLE} -m pytest -s -vvvv
        ${ImpactX_SOURCE_DIR}/tests/python
    WORKING_DIRECTORY
        ${CMAKE_PYTHON_OUTPUT_DIRECTORY}
)
