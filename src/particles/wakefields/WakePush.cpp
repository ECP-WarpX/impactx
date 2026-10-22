/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Alex Bojanich, Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "WakePush.H"
#include "InterpolateLinear1D.H"

#include <ablastr/particles/NodalFieldGather.H>

#include <AMReX_BLassert.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_REAL.H>
#include <AMReX_SPACE.H>


namespace impactx::particles::wakefields
{
    void WakePush (
        ImpactXParticleContainer & pc,
        amrex::Gpu::DeviceVector<amrex::Real> const & convolved_wakefield,
        amrex::ParticleReal slice_ds,
        amrex::Real bin_size,
        amrex::Real bin_min
    )
    {
        BL_PROFILE("impactx::particles::wakefields::WakePush")

        using namespace amrex::literals;

        int const num_values = static_cast<int>(convolved_wakefield.size());
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(num_values >= 2,
            "WakePush: the wake needs at least two values to interpolate");

        amrex::Real const inv_bin_size = 1_rt / bin_size;

        // Loop over refinement levels
        int const nLevel = pc.finestLevel();
        for (int lev = 0; lev <= nLevel; ++lev)
        {
            // Loop over all particle boxes
            using ParIt = ImpactXParticleContainer::iterator;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (ParIt pti(pc, lev); pti.isValid(); ++pti)
            {
                const int np = pti.numParticles();

                // Physical constants and reference quantities
                amrex::ParticleReal const mc_SI = pc.GetRefParticle().mass * (ablastr::constant::SI::c);
                amrex::ParticleReal const pz_ref_SI = pc.GetRefParticle().beta_gamma() * mc_SI;

                // Access data from StructOfArrays (soa)
                auto& soa_real = pti.GetStructOfArrays().GetRealData();

                amrex::ParticleReal* const AMREX_RESTRICT part_t = soa_real[RealSoA::t].dataPtr();
                amrex::ParticleReal* const AMREX_RESTRICT part_pt = soa_real[RealSoA::pt].dataPtr();

                // Obtain constants for force normalization
                amrex::ParticleReal const push_consts = 1_prt / ((ablastr::constant::SI::c) * pz_ref_SI);

                // Gather particles and push momentum
                const amrex::Real* wakefield_ptr = convolved_wakefield.data();
                amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i)
                {
                    // Access SoA Real data
                    amrex::ParticleReal const & AMREX_RESTRICT t = part_t[i];
                    amrex::ParticleReal & AMREX_RESTRICT pt = part_pt[i];

#if (defined(AMREX_DEBUG) || defined(DEBUG)) && !defined(AMREX_USE_GPU)
                    // Position in units of bin_size: the wake value k is at the bin edge v = k
                    amrex::Real const v = (t - bin_min) * inv_bin_size;
                    if (v < 0_rt || v > static_cast<amrex::Real>(num_values - 1))
                    {
                        std::cerr << "Warning: Position out of range for wakefield: " << v << std::endl;
                    }
#endif

                    // Update longitudinal momentum: interpolate the wake between the bin edges
                    amrex::ParticleReal const F_L = InterpolateLinear1D(
                        t, wakefield_ptr, num_values, bin_min, inv_bin_size);

#if (defined(AMREX_DEBUG) || defined(DEBUG)) && !defined(AMREX_USE_GPU)
                    // Check if the force (convolution) values are within a reasonable range
                    if (!std::isfinite(F_L))
                    {
                        // Handle unexpected values: log warning and skip momentum update
                        std::cerr << "Warning: Invalid or out-of-range values detected." << std::endl;
                    }
#endif

                    pt -= push_consts * slice_ds * F_L;
                    // Other dimensions (x, y) remain unchanged
                });
            } // End loop over all particle boxes
        } // End mesh-refinement level loop
    }
} // namespace impactx::particles::wakefields
