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

#include <ablastr/constant.H>

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>
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

    /** Interpolate the charge density profile and its slope to a longitudinal position
     *
     * DepositCharge1D bins with nearest-grid-point weighting, so value i of charge_distribution
     * belongs to the center of bin i, at bin_min + (i + 1/2) * bin_size, and the charge density
     * is zero outside of the bins. The slope is the central difference at the bin centers.
     * Both are interpolated linearly between the bin centers.
     *
     * @param[in] z longitudinal position
     * @param[in] charge_distribution deposited charge density in 1D along s
     * @param[in] num_values number of values in charge_distribution
     * @param[in] bin_min lower end of the beam in s
     * @param[in] bin_size size of one bin in s
     * @param[out] density charge density at z
     * @param[out] slope derivative of the number (or charge) density at z
     * @param[in] GetNumberDensity number density if true, otherwise charge density
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void InterpolateCharge1D (
        amrex::Real z,
        amrex::Real const * charge_distribution,
        int num_values,
        amrex::Real bin_min,
        amrex::Real bin_size,
        amrex::Real & density,
        amrex::Real & slope,
        bool GetNumberDensity = true
    )
    {
        using namespace amrex::literals;

        // could also be calculated on the host upfront
        amrex::Real const inv_bin_size = 1.0_rt / bin_size;

        // position relative to the bin centers: center j at or below z, weight w in [0, 1)
        amrex::Real const v = (z - bin_min) * inv_bin_size - 0.5_rt;
        int const j = static_cast<int>(amrex::Math::floor(v));
        amrex::Real const w = v - static_cast<amrex::Real>(j);

        // charge density at the bin centers j-1, j, j+1, j+2
        amrex::Real rho[4];
        for (int k = 0; k < 4; ++k) {
            int const bin = j - 1 + k;
            rho[k] = (bin >= 0 && bin < num_values) ? charge_distribution[bin] : 0.0_rt;
        }

        density = (1.0_rt - w) * rho[1] + w * rho[2];

        // central differences at the bin centers j and j+1
        amrex::Real const half_inv_bin_size = 0.5_rt * inv_bin_size;
        amrex::Real const slope_j = (rho[2] - rho[0]) * half_inv_bin_size;
        amrex::Real const slope_jp1 = (rho[3] - rho[1]) * half_inv_bin_size;
        slope = (1.0_rt - w) * slope_j + w * slope_jp1;

        if (GetNumberDensity)
        {
            constexpr amrex::Real inv_q_e = 1.0_rt / ablastr::constant::SI::q_e;
            slope *= inv_q_e;
        }
    }

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
