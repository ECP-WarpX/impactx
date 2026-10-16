/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Ji Qiang
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACT_DEPOSIT1D_PUSH_H
#define IMPACT_DEPOSIT1D_PUSH_H

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_REAL.H>


namespace impactx::particles::spacecharge
{
    /** Deposit charge using longitudinal binning for the calculation of space charge
     *
     *  This function performs particle binning for particles in the beam along the
     *  longitudinal direction, for the purposes of computing space charge in
     *  so-called 2.5D models. The num_bins bins have equal size and tile [bin_min, bin_max].
     *
     * @param[in] pc container of the particles that deposited rho
     * @param[in] bin_min location of lower endpoint for binning
     * @param[in] bin_max location of upper endpoint for binning
     * @param[in] num_bins number of longitudinal bins
     * @return beam_profile deposited charge density, one value per bin, summed over all MPI ranks
     */
    amrex::Gpu::DeviceVector<amrex::Real>
    Deposit1D (
        ImpactXParticleContainer & pc,
        amrex::Real bin_min,
        amrex::Real bin_max,
        int num_bins
    );

}  // namespace impactx::particles::spacecharge

#endif // IMPACT_DEPOSIT1D_H
