/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Alex Bojanich, Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef CHARGE_BINNING_H
#define CHARGE_BINNING_H

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_GpuContainers.H>
#include <AMReX_REAL.H>


namespace impactx::particles::wakefields
{
    /** Function to calculate charge density profile
     *
     * @param[in] myspc the particle species to deposit along s
     * @param[out] charge_distribution the output array (1D)
     * @param[in] bin_min lower end of the beam in s
     * @param[in] bin_size size of the beam in s divided by num_bins
     * @param[in] is_unity_particle_weight ignore the particle weighting per macro particle,
     *                                     otherwise treat each particle as one physical particle
     */
    void DepositCharge1D (
        impactx::ImpactXParticleContainer& myspc,
        amrex::Gpu::DeviceVector<amrex::Real> & charge_distribution,
        amrex::Real bin_min,
        amrex::Real bin_size,
        bool is_unity_particle_weight = false
    );

    /** Function to calculate the slope of the number (or charge) density
     *
     * @param[in] charge_distribution deposited charge in 1D along s
     * @param[out] slopes derivative of charge_distribution along s with len(charge_distribution)-1
     * @param[in] bin_size size of the beam in s divided by num_bins
     * @param[in] GetNumberDensity number density if true, otherwise charge density
     */
    void DerivativeCharge1D (
        amrex::Gpu::DeviceVector<amrex::Real> const & charge_distribution,
        amrex::Gpu::DeviceVector<amrex::Real> & slopes,
        amrex::Real bin_size,
        bool GetNumberDensity = true
    );

    /** Function to calculate the mean transverse positions (x and y) at each bin along s
     *
     * Bin k covers [bin_min + k * bin_size, bin_min + (k + 1) * bin_size), as in DepositCharge1D.
     * The means are over the particles of all MPI ranks: every rank must call this function.
     * Bins without particles get a mean of zero.
     *
     * @param[in] myspc the particle species to deposit along s
     * @param[in,out] mean_x the output array for mean x positions (1D), its size sets the number
     *                       of bins
     * @param[out] mean_y the output array for mean y positions (1D), the same size as mean_x
     * @param[in] bin_min lower end of the beam in s
     * @param[in] bin_size width of a bin in s
     * @param[in] is_unity_particle_weight ignore the particle weighting per macro particle,
     *                                     otherwise treat each particle as one physical particle
     */
    void MeanTransversePosition (
        impactx::ImpactXParticleContainer& myspc,
        amrex::Gpu::DeviceVector<amrex::Real> & mean_x,
        amrex::Gpu::DeviceVector<amrex::Real> & mean_y,
        amrex::Real bin_min,
        amrex::Real bin_size,
        bool is_unity_particle_weight = false
    );

} // namespace impactx::particles::wakefields

#endif // CHARGE_BINNING_H
