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
