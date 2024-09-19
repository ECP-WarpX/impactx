/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Alex Bojanich, Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "ChargeBinning.H"
#include "particles/ImpactXParticleContainer.H"

#include <cmath>


namespace impactx::particles::wakefields
{
    void DepositCharge1D (
        impactx::ImpactXParticleContainer& myspc,
        amrex::Gpu::DeviceVector<amrex::Real> & charge_distribution,
        amrex::Real bin_min,
        amrex::Real bin_size,
        bool is_unity_particle_weight
    )
    {
        int const num_bins = charge_distribution.size();
        amrex::Real * const dptr_data = charge_distribution.data();

        // Loop over each grid level
        int const nlevs = myspc.finestLevel();
        for (int lev = 0; lev <= nlevs; ++lev)
        {
            // TODO We could make this more efficient on CPU, see:
            // https://github.com/ECP-WarpX/WarpX/pull/5161/files#r1735068373
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            {
                // Loop over particles at the current grid level
                for (impactx::ParIterSoA pti(myspc, lev); pti.isValid(); ++pti)
                {
                    auto& soa = pti.GetStructOfArrays();  // Access data directly from StructOfArrays (soa)

                    // Number of particles
                    long const np = pti.numParticles();

                    // Access particle weights and momenta
                    amrex::ParticleReal* const AMREX_RESTRICT d_w = soa.GetRealData(impactx::RealSoA::w).dataPtr();

                    // Access particle positions
                    amrex::ParticleReal* const AMREX_RESTRICT pos_z = soa.GetRealData(impactx::RealSoA::z).dataPtr();

                    // Parallel loop over particles
                    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(int i)
                    {
                        // Access particle z-position directly
                        amrex::ParticleReal const z = pos_z[i];  // (Macro)Particle longitudinal position at i
                        auto const w = amrex::Real(d_w[i]);  // (Macro)Particle weight at i

                        /*
                        Weight w is given in [number of electrons]:
                        For w = 1 --> macroparticle = 1 electron making up a macroparticle
                        For w > 1 --> macroparticle > 1 electrons making up a macroparticle
                        */

                        // Calculate bin index based on z-position
                        int const bin = int(amrex::Math::floor((z - bin_min) / bin_size));  // Round to nearest bin index integer
                        if (bin < 0 || bin >= num_bins) { return; }  // Discard if bin index is out of range

                        // Divide charge by bin size to get binned charge density
                        amrex::Real add_value = ablastr::constant::SI::q_e / bin_size;

                        // Calculate the charge contribution of the macro particle
                        if (!is_unity_particle_weight)
                        {
                            add_value *= w;
                        }

                        // Add to histogram bin
                        amrex::HostDevice::Atomic::Add(&dptr_data[bin], add_value);
                    });
                }
            }
        }
    }

    void DerivativeCharge1D (
        amrex::Gpu::DeviceVector<amrex::Real> const & charge_distribution,
        amrex::Gpu::DeviceVector<amrex::Real> & slopes,
        amrex::Real bin_size,
        bool GetNumberDensity
    )
    {
        int const num_bins = charge_distribution.size();
        amrex::Real const * const dptr_charge_distribution = charge_distribution.data();
        amrex::Real * const dptr_slopes = slopes.data();

        amrex::ParallelFor(num_bins - 1, [=] AMREX_GPU_DEVICE(int i)
        {
            // Compute the charge density derivative
            amrex::Real const charge_derivative = (dptr_charge_distribution[i + 1] - dptr_charge_distribution[i]) / bin_size;

            // If GetNumberDensity = True, convert charge density derivative to number density derivative for CSR convolution
            if (GetNumberDensity)
            {
                dptr_slopes[i] = charge_derivative / ablastr::constant::SI::q_e;
            }
            else
            {
                dptr_slopes[i] = charge_derivative;
            }
        });
    }

    void MeanTransversePosition (
        impactx::ImpactXParticleContainer& myspc,
        amrex::Gpu::DeviceVector<amrex::Real> & mean_x,
        amrex::Gpu::DeviceVector<amrex::Real> & mean_y,
        amrex::Real bin_min,
        amrex::Real bin_size,
        bool is_unity_particle_weight
    )
    {
        using namespace amrex::literals;

        int const num_bins = mean_x.size();
        amrex::Real* dptr_mean_x = mean_x.data();
        amrex::Real* dptr_mean_y = mean_y.data();

        // Declare arrays for sums of positions and weights
        auto sum_x = amrex::Gpu::DeviceVector<amrex::Real>(num_bins, 0.0_rt);
        auto sum_y = amrex::Gpu::DeviceVector<amrex::Real>(num_bins, 0.0_rt);
        auto sum_w = amrex::Gpu::DeviceVector<amrex::Real>(num_bins, 0.0_rt);

        amrex::Real* const sum_w_ptr = sum_w.data();
        amrex::Real* const sum_x_ptr = sum_x.data();
        amrex::Real* const sum_y_ptr = sum_y.data();

        int const nlevs = myspc.finestLevel();
        for (int lev = 0; lev <= nlevs; ++lev)
        {
            // TODO We could make this more efficient on CPU, see:
            // https://github.com/ECP-WarpX/WarpX/pull/5161/files#r1735068373
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            {
                for (impactx::ParIterSoA pti(myspc, lev); pti.isValid(); ++pti)
                {
                    auto& soa = pti.GetStructOfArrays();
                    long const np = pti.numParticles();

                    amrex::Real* const AMREX_RESTRICT pos_x = soa.GetRealData(impactx::RealSoA::x).dataPtr();
                    amrex::Real* const AMREX_RESTRICT pos_y = soa.GetRealData(impactx::RealSoA::y).dataPtr();
                    amrex::Real* const AMREX_RESTRICT pos_z = soa.GetRealData(impactx::RealSoA::z).dataPtr();
                    amrex::Real* const AMREX_RESTRICT d_w = soa.GetRealData(impactx::RealSoA::w).dataPtr();

                    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(int i)
                    {
                        amrex::Real const w = d_w[i];
                        amrex::Real const x = pos_x[i];
                        amrex::Real const y = pos_y[i];
                        amrex::Real const z = pos_z[i];

                        int const bin = int(amrex::Math::floor((z - bin_min) / bin_size));
                        if (bin < 0 || bin >= num_bins) { return; }

                        amrex::Real const weight = is_unity_particle_weight ? 1.0_rt : w;  // Check is macroparticle made up of 1 or more particles

                        amrex::HostDevice::Atomic::Add(&sum_w_ptr[bin], weight);      // Deposit the number of particles composing macroparticle
                        amrex::HostDevice::Atomic::Add(&sum_x_ptr[bin], x * weight);  // Deposit x position multiplied by number of particles at this position
                        amrex::HostDevice::Atomic::Add(&sum_y_ptr[bin], y * weight);
                    });
                }
            }
        }

        amrex::ParallelFor(num_bins, [=] AMREX_GPU_DEVICE(int i)
        {
            if (sum_w_ptr[i] > 0) // Ensure number of particles in a bin is >= 1 before taking the mean
            {
                dptr_mean_x[i] = sum_x_ptr[i] / sum_w_ptr[i];
                dptr_mean_y[i] = sum_y_ptr[i] / sum_w_ptr[i];
            }
            else
            {
                dptr_mean_x[i] = 0.0;
                dptr_mean_y[i] = 0.0;
            }
        });
    }
}
