target_sources(lib
  PRIVATE
    ChargeDeposition.cpp
    CollectLost.cpp
    ImpactXParticleContainer.cpp
    Push.cpp
)

add_subdirectory(diagnostics)
add_subdirectory(elements)
add_subdirectory(spacecharge)
add_subdirectory(transformation)
add_subdirectory(wakefields)
