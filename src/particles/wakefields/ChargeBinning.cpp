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
#include "../ImpactXParticleContainer.H" //Includes all necessary AMReX headers

#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace amrex;

//Beam charge distribution function

void DepositCharge1D(impactx::ImpactXParticleContainer& myspc, Real* dptr_data, int num_bins, Real bin_min, Real bin_size, bool is_unity_particle_weight)
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
                //amrex::ParticleReal* const AMREX_RESTRICT pos_x = soa.GetRealData(impactx::RealSoA::x).dataPtr();
                //amrex::ParticleReal* const AMREX_RESTRICT pos_y = soa.GetRealData(impactx::RealSoA::y).dataPtr();
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

                    //Calculate bin index based on z-position
                    int const bin = int(Math::floor((z - bin_min) / bin_size)); //Round to nearest bin index integer
                    if (bin < 0 || bin >= num_bins) { return; } //Discard if bin index is out of range

                    //Calculate the charge contribution of the macro particle
                    Real charge_contribution = w * ablastr::constant::SI::q_e;

                    //Divide charge by bin size to get binned charge density
                    Real charge_density = charge_contribution / bin_size;
                    Real unity_density = ablastr::constant::SI::q_e / bin_size;

                    //Add to histogram bin
                    if (is_unity_particle_weight)
                    {
                        HostDevice::Atomic::Add(&dptr_data[bin], unity_density);  //Unity weight: Add elementary charge of 1 particle
                    }
                    else
                    {
                        HostDevice::Atomic::Add(&dptr_data[bin], charge_density);  //Non-unity weight: Add calculated charge from >1 particle
                    }
                });
            }
        }
    }
}

void DerivativeCharge1D(double* charge_distribution, double* slopes, int num_bins, Real bin_size, bool GetNumberDensity)
{

    for (int i = 0; i < num_bins - 1; ++i)
    {
        //Compute the charge density derivative
        double charge_derivative = (charge_distribution[i + 1] - charge_distribution[i]) / bin_size;

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
