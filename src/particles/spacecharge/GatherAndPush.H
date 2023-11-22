/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Remi Lehe
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_GATHER_AND_PUSH_H
#define IMPACTX_GATHER_AND_PUSH_H

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#include <unordered_map>
#include <string>


namespace impactx::spacecharge
{
    /** Gather force fields and push particles in x,y,z
     *
     * This gathers the space charge field with respect to particle position
     * and shape. The momentum of all particles is then pushed using a common
     * time step given by the reference particle speed and ds slice. The
     * position push is done in the lattice elements and not here.
     *
     * @param[inout] pc container of the particles that deposited rho
     * @param[in] space_charge_field space charge force component in x,y,z per level
     * @param[in] geom geometry object
     * @param[in] slice_ds segment length in meters
     */
    void GatherAndPush (
        ImpactXParticleContainer & pc,
        std::unordered_map<int, std::unordered_map<std::string, amrex::MultiFab> > const & space_charge_field,
        const amrex::Vector<amrex::Geometry>& geom,
        amrex::ParticleReal slice_ds
    );

} // namespace impactx

#endif // IMPACTX_GATHER_AND_PUSH_H
