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

#include "diagnostics/ReducedBeamCharacteristics.H"
#include "particles/wakefields/ChargeBinning.H"

#include <AMReX_REAL.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_GpuParallelReduce.H>
#include <AMReX_ParmParse.H>

#include <cmath>


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

        // Set parameters for charge deposition
        bool const is_unity_particle_weight = false;

        // num_bins bins of equal size tile [bin_min, bin_max]
        amrex::Real const bin_size = (bin_max - bin_min) / num_bins;
        // Allocate memory for the charge profile, with one more bin for the particles at bin_max
        amrex::Gpu::DeviceVector<amrex::Real> charge_distribution(num_bins + 1, 0.0);
        // Call charge deposition function
        impactx::particles::wakefields::DepositCharge1D(pc, charge_distribution, bin_min, bin_size, is_unity_particle_weight);

        // The particles at bin_max belong to the last bin
        amrex::Real * const dptr_charge = charge_distribution.data();
        amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE (int)
        {
            dptr_charge[num_bins - 1] += dptr_charge[num_bins];
        });
        charge_distribution.resize(num_bins);

        // Sum up all partial charge histograms to each MPI process to calculate
        // the global charge slope.
        amrex::ParallelAllReduce::Sum(
            charge_distribution,
            amrex::ParallelDescriptor::Communicator()
        );

        return charge_distribution;
    }

}  // namespace impactx::particles::spacecharge
