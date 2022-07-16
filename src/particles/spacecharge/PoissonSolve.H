/* Copyright 2022 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_POISSONSOLVE_H
#define IMPACTX_POISSONSOLVE_H

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_MultiFab.H>

#include <unordered_map>


namespace impactx::spacecharge
{
    /** Calculate the electric potential from charge density
     *
     * This resets the values in phi to zero and then calculates the space
     * charge potential phi.
     *
     * @param[in] pc container of the particles that deposited rho
     * @param[in] rho charge per level
     * @param[inout] phi scalar potential per level
     */
    void PoissonSolve (
        ImpactXParticleContainer const & pc,
        std::unordered_map<int, amrex::MultiFab> & rho,
        std::unordered_map<int, amrex::MultiFab> & phi
    );

} // namespace impactx

#endif // IMPACTX_POISSONSOLVE_H
