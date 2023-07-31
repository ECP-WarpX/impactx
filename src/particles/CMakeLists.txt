target_sources(lib
  PRIVATE
    ChargeDeposition.cpp
    ImpactXParticleContainer.cpp
    Push.cpp
)

add_subdirectory(diagnostics)
add_subdirectory(elements)
add_subdirectory(spacecharge)
add_subdirectory(transformation)
