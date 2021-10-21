target_sources(ImpactX
  PRIVATE
    ImpactX.cpp
)

if(ImpactX_APP)
    target_sources(app PRIVATE main.cpp)
endif()
