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
#include <vector>

using namespace amrex;


namespace impactx::particles::wakefields
{
    void DepositCharge1D (
        impactx::ImpactXParticleContainer& myspc,
        Real* dptr_data,
        int num_bins,
        Real bin_min,
        Real bin_size,
        bool is_unity_particle_weight
    )
    {   //Access and change data directly using '&' to pass by reference

        //Determine the number of grid levels in the simulation
        int const nlevs = std::max(0, myspc.finestLevel() + 1);

        //Loop over each grid level
        for (int lev = 0; lev < nlevs; ++lev)
        {
            //OpenMP parallelization if enabled
            #ifdef AMREX_USE_OMP
            #pragma omp parallel if (Gpu::notInLaunchRegion())
            #endif
            {
                //Loop over particles at the current grid level
                for (impactx::ParIterSoA pti(myspc, lev); pti.isValid(); ++pti)
                {
                    auto& soa = pti.GetStructOfArrays(); //Access data directly from StructOfArrays (soa)

                    //Number of particles
                    long const np = pti.numParticles();

                    //Access particle weights and momenta
                    amrex::ParticleReal* const AMREX_RESTRICT d_w = soa.GetRealData(impactx::RealSoA::w).dataPtr();

                    //Access particle positions
                    amrex::ParticleReal* const AMREX_RESTRICT pos_z = soa.GetRealData(impactx::RealSoA::z).dataPtr();

                    //Parallel loop over particles
                    ParallelFor(np, [=] AMREX_GPU_DEVICE(int i)
                    {
                        //Access particle z-position directly
                        ParticleReal z = pos_z[i]; //(Macro)Particle longitudinal position at i
                        auto const w = (Real)d_w[i]; //(Macro)Particle weight at i

                        /*
                        Weight w is given in [number of electrons]:
                        For w = 1 --> macroparticle = 1 electron making up a macroparticle
                        For w > 1 --> macroparticle > 1 electrons making up a macroparticle
                        */

                        // Calculate bin index based on z-position
                        int const bin = int(Math::floor((z - bin_min) / bin_size)); //Round to nearest bin index integer
                        if (bin < 0 || bin >= num_bins) { return; } //Discard if bin index is out of range

                        // Divide charge by bin size to get binned charge density
                        Real add_value = ablastr::constant::SI::q_e / bin_size;

                        // Calculate the charge contribution of the macro particle
                        if (!is_unity_particle_weight)
                        {
                            add_value *= w;
                        }

                        // Add to histogram bin
                        HostDevice::Atomic::Add(&dptr_data[bin], add_value);
                    });
                }
            }
        }
    }

    void DerivativeCharge1D (
        amrex::Real* charge_distribution,
        amrex::Real* slopes,
        int num_bins,
        Real bin_size,
        bool GetNumberDensity
    )
    {

        for (int i = 0; i < num_bins - 1; ++i)
        {
            //Compute the charge density derivative
            amrex::Real charge_derivative = (charge_distribution[i + 1] - charge_distribution[i]) / bin_size;

            //If GetNumberDensity = True, convert charge density derivative to number density derivative for CSR convolution
            if (GetNumberDensity)
            {
                slopes[i] = charge_derivative / ablastr::constant::SI::q_e;
            }
            else
            {
                slopes[i] = charge_derivative;
            }
        }
    }

    void MeanTransversePosition (
        impactx::ImpactXParticleContainer& myspc,
        amrex::Real* mean_x,
        amrex::Real* mean_y,
        int num_bins,
        amrex::Real bin_min,
        amrex::Real bin_size,
        bool is_unity_particle_weight
    )
    {
        int const nlevs = std::max(0, myspc.finestLevel() + 1);

        // Declare arrays for sums of positions and weights
        std::vector<amrex::Real> sum_x(num_bins, 0.0);
        std::vector<amrex::Real> sum_y(num_bins, 0.0);
        std::vector<amrex::Real> sum_w(num_bins, 0.0);

        amrex::Real* sum_w_ptr = sum_w.data();
        amrex::Real* sum_x_ptr = sum_x.data();
        amrex::Real* sum_y_ptr = sum_y.data();

        for (int lev = 0; lev < nlevs; ++lev)
        {
            #ifdef AMREX_USE_OMP
            #pragma omp parallel if (Gpu::notInLaunchRegion())
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

                    ParallelFor(np, [=] AMREX_GPU_DEVICE(int i)
                    {
                        amrex::Real w = d_w[i];
                        amrex::Real x = pos_x[i];
                        amrex::Real y = pos_y[i];
                        amrex::Real z = pos_z[i];

                        int const bin = int(Math::floor((z - bin_min) / bin_size));
                        if (bin < 0 || bin >= num_bins) { return; }

                        amrex::Real weight = is_unity_particle_weight ? 1.0 : w; // Check is macroparticle made up of 1 or more particles

                        sum_w_ptr[bin] += weight; // Deposit the number of particles composing macroparticle
                        sum_x_ptr[bin] += x * weight; // Deposit x position multiplied by number of particles at this position
                        sum_y_ptr[bin] += y * weight;

                        //HostDevice::Atomic::Add(&sum_w[bin], weight); // Deposit the number of particles composing macroparticle
                        //HostDevice::Atomic::Add(&sum_x[bin], x * weight); // Deposit x position multiplied by number of particles at this position
                        //HostDevice::Atomic::Add(&sum_y[bin], y * weight);
                    });
                }
            }
        }

        for (int i = 0; i < num_bins; ++i)
        {
            if (sum_w[i] > 0) // Ensure number of particles in a bin is >= 1 before taking the mean
            {
                mean_x[i] = sum_x_ptr[i] / sum_w_ptr[i];
                mean_y[i] = sum_y_ptr[i] / sum_w_ptr[i];
            }
            else
            {
                mean_x[i] = 0.0;
                mean_y[i] = 0.0;
            }
        }
    }
}
