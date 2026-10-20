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

#include <AMReX_GpuContainers.H>
#include <AMReX_REAL.H>


namespace impactx::particles::spacecharge
{
    /** Deposit charge using longitudinal binning for the calculation of space charge
     *
     *  This function deposits the particles of the beam along the longitudinal direction, for
     *  the purposes of computing space charge in so-called 2.5D models. The num_bins grid points
     *  are equidistant, at bin_min + k * bin_size with bin_size = (bin_max - bin_min) /
     *  (num_bins - 1), so the first and the last grid point lie on bin_min and bin_max.
     *  Each particle is deposited with linear (cloud-in-cell) weights onto its two neighboring
     *  grid points, the same weights as the linear interpolation back to the particles.
     *  A particle deposits its charge magnitude |q| w, with the charge q of the reference
     *  particle and the weight w as the number of real particles it represents.
     *
     * @param[in] pc container of the particles that deposited rho
     * @param[in] bin_min location of lower endpoint for binning
     * @param[in] bin_max location of upper endpoint for binning
     * @param[in] num_bins number of longitudinal grid points, at least two
     * @return beam_profile deposited charge magnitude per grid point divided by bin_size, summed
     *                      over all MPI ranks
     */
    amrex::Gpu::DeviceVector<amrex::Real>
    Deposit1D (
        ImpactXParticleContainer & pc,
        amrex::Real bin_min,
        amrex::Real bin_max,
        int num_bins
    );

    /** Charge density and its slope at the grid points of Deposit1D
     *
     * The first and the last grid point lie on the ends of the bunch and own only half a
     * bin, so their deposited charge density is doubled for the density, which keeps the line
     * density of a flat-top bunch flat up to its ends. The slope is the central difference of
     * the deposited charge density, with zero charge density outside of the grid points.
     * Interpolated linearly to the particles (InterpolateLinear1D), this slope conserves the
     * total longitudinal momentum exactly.
     *
     * @param[in] charge_distribution deposited charge density from Deposit1D, at least two values
     * @param[out] density charge density at the grid points, len(charge_distribution) values
     * @param[out] slope derivative of the charge density at the grid points,
     *                   len(charge_distribution) values
     * @param[in] bin_size distance between two grid points
     */
    void DensityAndSlope1D (
        amrex::Gpu::DeviceVector<amrex::Real> const & charge_distribution,
        amrex::Gpu::DeviceVector<amrex::Real> & density,
        amrex::Gpu::DeviceVector<amrex::Real> & slope,
        amrex::Real bin_size
    );

}  // namespace impactx::particles::spacecharge

#endif // IMPACT_DEPOSIT1D_H
