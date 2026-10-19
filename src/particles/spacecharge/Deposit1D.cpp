/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Ji Qiang
 * License: BSD-3-Clause-LBNL
 */
#include "Deposit1D.H"

#include "particles/wakefields/InterpolateLinear1D.H"

#include <ablastr/constant.H>

#include <AMReX_BLassert.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuParallelReduce.H>
#include <AMReX_REAL.H>


namespace impactx::particles::spacecharge
{

    amrex::Gpu::DeviceVector<amrex::Real>
    Deposit1D (
        ImpactXParticleContainer & pc,
        amrex::Real bin_min,
        amrex::Real bin_max,
        int num_bins
    )
    {
        BL_PROFILE("impactx::spacecharge::Deposit1D");

        using namespace amrex::literals;

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(num_bins >= 2,
            "Deposit1D: needs at least two grid points");

        // num_bins equidistant grid points, the first on bin_min and the last on bin_max
        amrex::Real const bin_size = (bin_max - bin_min) / (num_bins - 1);
        amrex::Real const inv_bin_size = 1.0_rt / bin_size;

        // charge density per unit particle weight
        amrex::Real const charge_density_per_weight = ablastr::constant::SI::q_e * inv_bin_size;

        amrex::Gpu::DeviceVector<amrex::Real> charge_distribution(num_bins, 0.0_rt);
        amrex::Real * const dptr_charge = charge_distribution.data();

        // Loop over each grid level
        int const nlevs = pc.finestLevel();
        for (int lev = 0; lev <= nlevs; ++lev)
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (impactx::ParIterSoA pti(pc, lev); pti.isValid(); ++pti)
            {
                auto & soa = pti.GetStructOfArrays();
                long const np = pti.numParticles();

                amrex::ParticleReal const * const AMREX_RESTRICT part_w =
                    soa.GetRealData(RealSoA::w).dataPtr();
                amrex::ParticleReal const * const AMREX_RESTRICT part_z =
                    soa.GetRealData(RealSoA::z).dataPtr();

                // amrex::For: different particles add to the same grid points
                amrex::For(np, [=] AMREX_GPU_DEVICE (long i)
                {
                    // linear (cloud-in-cell) weights, the same as InterpolateLinear1D
                    int k = 0;
                    amrex::Real w = 0.0_rt;
                    wakefields::LinearWeights1D(part_z[i], num_bins, bin_min, inv_bin_size, k, w);

                    amrex::Real const charge_density = charge_density_per_weight * part_w[i];
                    amrex::HostDevice::Atomic::Add(&dptr_charge[k], (1.0_rt - w) * charge_density);
                    amrex::HostDevice::Atomic::Add(&dptr_charge[k + 1], w * charge_density);
                });
            }
        }

        // Sum up all partial charge densities to each MPI process to calculate
        // the global charge slope.
        amrex::ParallelAllReduce::Sum(
            charge_distribution,
            amrex::ParallelDescriptor::Communicator()
        );

        return charge_distribution;
    }

    void DensityAndSlope1D (
        amrex::Gpu::DeviceVector<amrex::Real> const & charge_distribution,
        amrex::Gpu::DeviceVector<amrex::Real> & density,
        amrex::Gpu::DeviceVector<amrex::Real> & slope,
        amrex::Real bin_size,
        bool GetNumberDensity
    )
    {
        using namespace amrex::literals;

        int const num_points = static_cast<int>(charge_distribution.size());
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(num_points >= 2,
            "DensityAndSlope1D: needs at least two grid points");
        density.resize(num_points);
        slope.resize(num_points);

        // central difference, optionally of the number density
        amrex::Real const slope_factor = GetNumberDensity
            ? 0.5_rt / (bin_size * ablastr::constant::SI::q_e)
            : 0.5_rt / bin_size;

        amrex::Real const * const dptr_charge = charge_distribution.data();
        amrex::Real * const dptr_density = density.data();
        amrex::Real * const dptr_slope = slope.data();
        amrex::ParallelFor(num_points, [=] AMREX_GPU_DEVICE(int k)
        {
            // the first and the last grid point own only half a bin
            bool const is_end = k == 0 || k == num_points - 1;
            dptr_density[k] = is_end ? 2_rt * dptr_charge[k] : dptr_charge[k];

            // zero charge density outside of the grid points, so the slopes next to the bunch
            // ends include its steps there
            amrex::Real const charge_below = (k > 0) ? dptr_charge[k - 1] : 0_rt;
            amrex::Real const charge_above = (k + 1 < num_points) ? dptr_charge[k + 1] : 0_rt;
            dptr_slope[k] = (charge_above - charge_below) * slope_factor;
        });
    }

}  // namespace impactx::particles::spacecharge
