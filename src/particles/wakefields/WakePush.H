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

#include <AMReX_GpuContainers.H>
#include <AMReX_REAL.H>


namespace impactx::particles::wakefields
{
    /** Function to update particle momentum using wake force kick
     *
     * @param[inout] pc the particle species container
     * @param[in] convolved_wakefield wake functions in [V*pc/mm]
     * @param[in] slice_ds slice spacing along s
     * @param[in] bin_size size of the beam in s divided by num_bins
     * @param[in] bin_min lower end of the beam in s
     */
    void WakePush (
        ImpactXParticleContainer & pc,
        amrex::Gpu::DeviceVector<amrex::Real> const & convolved_wakefield,
        amrex::ParticleReal slice_ds,
        amrex::Real bin_size,
        amrex::Real bin_min
    );

} // namespace impactx::particles::wakefields

#endif // WAKEPUSH_H
