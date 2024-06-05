#ifndef WAKEPUSH_H
#define WAKEPUSH_H

#include <AMReX_REAL.H>
#include <vector>
#include "../ImpactXParticleContainer.H" 

namespace impactx::wakepush
{
    //Function to update particle momentum using wake force kick
    void WakePush(ImpactXParticleContainer & pc, const std::vector<double>& convoluted_wakefield, amrex::ParticleReal const slice_ds);
}

#endif // WAKEPUSH_H
