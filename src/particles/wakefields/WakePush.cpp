#include "WakePush.H"

#include <ablastr/particles/NodalFieldGather.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_REAL.H>
#include <AMReX_SPACE.H>

namespace impactx::wakepush
{
    void WakePush(ImpactXParticleContainer & pc, const std::vector<double>& convoluted_wakefield, amrex::ParticleReal const slice_ds)
    {
        BL_PROFILE("impactx::wakepush::WakePush");

        using namespace amrex::literals;

        //Loop over refinement levels
        int const nLevel = pc.finestLevel();
        for (int lev = 0; lev <= nLevel; ++lev)
        {
            //Loop over all particle boxes
            using ParIt = ImpactXParticleContainer::iterator;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif

            for (ParIt pti(pc, lev); pti.isValid(); ++pti)
            {
                const int np = pti.numParticles();

                //Physical constants and reference quantities
                amrex::ParticleReal const mc_SI = pc.GetRefParticle().mass * (ablastr::constant::SI::c);
                amrex::ParticleReal const pz_ref_SI = pc.GetRefParticle().beta_gamma() * mc_SI;
                amrex::ParticleReal const gamma = pc.GetRefParticle().gamma();
                amrex::ParticleReal const inv_gamma2 = 1.0_prt / (gamma * gamma);

                amrex::ParticleReal const dt = slice_ds / pc.GetRefParticle().beta() / c0_SI;

                //Access data from StructOfArrays (soa)
                auto& soa_real = pti.GetStructOfArrays().GetRealData();

                amrex::ParticleReal* const AMREX_RESTRICT part_x = soa_real[RealSoA::x].dataPtr();
                amrex::ParticleReal* const AMREX_RESTRICT part_y = soa_real[RealSoA::y].dataPtr();
                amrex::ParticleReal* const AMREX_RESTRICT part_z = soa_real[RealSoA::z].dataPtr(); //Note: Currently for a fixed t

                amrex::ParticleReal* const AMREX_RESTRICT part_px = soa_real[RealSoA::px].dataPtr();
                amrex::ParticleReal* const AMREX_RESTRICT part_py = soa_real[RealSoA::py].dataPtr();
                amrex::ParticleReal* const AMREX_RESTRICT part_pz = soa_real[RealSoA::pz].dataPtr(); //Note: Currently for a fixed t

                //Obtain constants for force normalization
                amrex::ParticleReal const push_consts = 1.0 / ((ablastr::constant::SI::c) * pz_ref_SI);

                //Gather particles and push momentum
                amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i)
                {
                    //Access SoA Real data
                    amrex::ParticleReal & AMREX_RESTRICT x = part_x[i];
                    amrex::ParticleReal & AMREX_RESTRICT y = part_y[i];
                    amrex::ParticleReal & AMREX_RESTRICT z = part_z[i];

                    amrex::ParticleReal & AMREX_RESTRICT px = part_px[i];
                    amrex::ParticleReal & AMREX_RESTRICT py = part_py[i];
                    amrex::ParticleReal & AMREX_RESTRICT pz = part_pz[i];

                    //Update longitudinal momentum with the convoluted wakefield force
                    int idx = static_cast<int>((z - slice_ds) / slice_ds * convoluted_wakefield.size()); //Find index position along z
                    idx = std::min(std::max(idx, 0), static_cast<int>(convoluted_wakefield.size()) - 1); //Ensure index in range of wakefield force vector
                    amrex::ParticleReal const F_L = convoluted_wakefield[idx];

                    //Update longitudinal momentum
                    pz -=  push_consts * slice_ds * F_L;

                    //Other dimensions (x, y) remain unchanged
                    // px = px;
                    // py = py;
                });
            } //End loop over all particle boxes
        } //End mesh-refinement level loop
    }
}
