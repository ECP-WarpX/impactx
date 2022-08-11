/* Copyright 2022 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Marco Garten, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_FORCEFROMSELFFIELDS_H
#define IMPACTX_FORCEFROMSELFFIELDS_H

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#include <unordered_map>


namespace impactx::spacecharge
{
    /** Calculate the space charge force field from the electric potential
     *
     * This resets the values in scf_<component> to zero and then calculates the space
     * charge force field.
     *
     * @param[inout] space_charge_force space charge force component in x,y,z per level
     * @param[in] phi scalar potential per level
     * @param[in] geom geometry object
     */
    void ForceFromSelfFields (
        std::unordered_map<int, std::unordered_map<std::string, amrex::MultiFab> > & space_charge_force,
        std::unordered_map<int, amrex::MultiFab> const & phi,
        const amrex::Vector<amrex::Geometry>& geom
    );

} // namespace impactx

#endif // IMPACTX_FORCEFROMSELFFIELDS_H
