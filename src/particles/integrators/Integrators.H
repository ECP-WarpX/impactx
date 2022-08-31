/* Copyright 2022 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_INTEGRATORS_H
#define IMPACTX_INTEGRATORS_H

#include "particles/ImpactXParticleContainer.H"

#include <string>


namespace impactx::integrators
{
    /** Comments here
     *
     */

    void symp4_integrate (RefPart & refpart,
                           amrex::ParticleReal const zin,
                           amrex::ParticleReal const zout,
                           int const nsteps, std::string element_type);

    void map1 (amrex::ParticleReal const tau, RefPart & refpart,
                           amrex::ParticleReal & zeval, std::string element_type);

    void map2 (amrex::ParticleReal const tau, RefPart & refpart,
                           amrex::ParticleReal & zeval, std::string element_type);


} // namespace impactx::integrators

#endif // IMPACTX_INTEGRATORS_H
