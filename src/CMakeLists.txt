###############################################################################
# These are the files with the central application logic.
# We build this into a library internally, so we can call this both as an
# application as well as a Python module.
#
target_sources(ImpactX
  PRIVATE
    ImpactX.cpp
)

add_subdirectory(initialization)
add_subdirectory(particles)

###############################################################################
# This is the entry point of the application.
# This fail primarily contains main(), so we can run ImpactX as an application
# on the command line.
#
if(ImpactX_APP)
    target_sources(app PRIVATE main.cpp)
endif()
