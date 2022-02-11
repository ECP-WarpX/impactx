target_sources(ImpactX
  PRIVATE
    ImpactXParticleContainer.cpp
    Push.cpp
)

add_subdirectory(transformation)
add_subdirectory(diagnostics)
