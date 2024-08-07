/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
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

namespace impactx::particles::wakefields
{
    /** Function to calculate charge density profile
     *
     * @param[in] myspc the particle species to deposit along s
     * @param[out] dptr_data the output array (1D)
     * @param[in] num_bins number of bins in dptr_data
     * @param[in] bin_min lower end of the beam in s
     * @param[in] bin_size size of the beam in s divided by num_bins
     * @param[in] is_unity_particle_weight ignore the particle weighting per macro particle,
     *                                     otherwise treat each particle as one physical particle
     */
    void DepositCharge1D (
        impactx::ImpactXParticleContainer& myspc,
        amrex::Real* dptr_data,
        int num_bins,
        amrex::Real bin_min,
        amrex::Real bin_size,
        bool is_unity_particle_weight = false
    );

    /** Function to calculate the slope of the number (or charge) density
     *
     * @param[in] charge_distribution deposited charge in 1D along s
     * @param[out] slopes derivative of charge_distribution along s
     * @param[in] num_bins number of bins in charge_distribution and in len(slopes)-1
     * @param[in] bin_size size of the beam in s divided by num_bins
     * @param[in] GetNumberDensity number density if true, otherwise charge density
     */
    void DerivativeCharge1D (
        amrex::Real* charge_distribution,
        amrex::Real* slopes,
        int num_bins,
        amrex::Real bin_size,
        bool GetNumberDensity = true
    );

} // namespace impactx::particles::wakefields

#endif // CHARGE_BINNING_H
