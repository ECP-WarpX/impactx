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
