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

    /** Function to calculate the charge density and its slope at the bin centers
     *
     * The bins tile the bunch (see Deposit1D): value i of charge_distribution is the charge
     * density at the center of bin i, at bin_min + (i + 1/2) * bin_size. Both outputs have one
     * more value on each side, so that linear interpolation with InterpolateLinear1D covers the
     * whole bunch: value k belongs to bin_min + (k - 1/2) * bin_size.
     *
     * The slope is the central difference, with zero charge density outside of the bins, so it
     * includes the steps at the bunch ends. The density outside of the bins is extrapolated
     * linearly from the two outermost bins, limited to [0, outermost bin], so that the line
     * density of a flat-top bunch stays flat up to its ends.
     *
     * @param[in] charge_distribution deposited charge density in 1D along s, at least one bin
     * @param[out] density charge density at the bin centers, len(charge_distribution)+2 values
     * @param[out] slope derivative of the number (or charge) density at the bin centers,
     *                   len(charge_distribution)+2 values
     * @param[in] bin_size size of one bin in s
     * @param[in] GetNumberDensity number density if true, otherwise charge density
     */
    void DensityAndSlope1D (
        amrex::Gpu::DeviceVector<amrex::Real> const & charge_distribution,
        amrex::Gpu::DeviceVector<amrex::Real> & density,
        amrex::Gpu::DeviceVector<amrex::Real> & slope,
        amrex::Real bin_size,
        bool GetNumberDensity = true
    );

    /** Function to calculate the mean transverse positions (x and y) at each bin along s
     *
     * @param[in] myspc the particle species to deposit along s
     * @param[out] mean_x the output array for mean x positions (1D)
     * @param[out] mean_y the output array for mean y positions (1D)
     * @param[in] bin_min lower end of the beam in s
     * @param[in] bin_size size of the beam in s divided by num_bins
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
