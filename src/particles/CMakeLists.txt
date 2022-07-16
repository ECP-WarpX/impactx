target_sources(ImpactX
  PRIVATE
    ChargeDeposition.cpp
    ImpactXParticleContainer.cpp
    Push.cpp
)

add_subdirectory(diagnostics)
add_subdirectory(spacecharge)
add_subdirectory(transformation)
