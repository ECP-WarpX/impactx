#ifndef WAKEPUSH_H
#define WAKEPUSH_H

#include <AMReX_REAL.H>
#include <vector>
#include "../ImpactXParticleContainer.H"

namespace impactx::wakepush
{
    //Function to update particle momentum using wake force kick
    void WakePush(ImpactXParticleContainer & pc, std::vector<double>& convoluted_wakefield, amrex::ParticleReal const slice_ds, double bin_size, double bin_min, int padding_factor);
}

#endif // WAKEPUSH_H
