/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Alex Bojanich, Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef WAKEPUSH_H
#define WAKEPUSH_H

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_REAL.H>

#include <vector>


namespace impactx::particles::wakefields
{
    /** Function to update particle momentum using wake force kick
     *
     * @param[inout] pc
     * @param[in] convoluted_wakefield
     * @param[in] slice_ds
     * @param[in] bin_size
     * @param[in] bin_min
     * @param[in] padding_factor
     */
    void WakePush (
        ImpactXParticleContainer & pc,
        amrex::Gpu::DeviceVector<amrex::Real> const & convoluted_wakefield,
        amrex::ParticleReal slice_ds,
        amrex::Real bin_size,
        amrex::Real bin_min,
        int padding_factor
    );

} // namespace impactx::particles::wakefields

#endif // WAKEPUSH_H
