###############################################################################
# These are the files equivalent to the primary C++ implementation - but here
# we define how they will appear in our Python module (aka Python bindings).
#
target_sources(pyImpactX
  PRIVATE
    distribution.cpp
    elements.cpp
    ImpactX.cpp
    ImpactXParticleContainer.cpp
    ReferenceParticle.cpp
    transformation.cpp
    WakeConvolution.cpp
)
