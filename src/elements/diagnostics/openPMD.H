/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_OPENPMD_H
#define IMPACTX_OPENPMD_H

#include "particles/ReferenceParticle.H"

#ifdef ImpactX_USE_OPENPMD
#   include <openPMD/openPMD.hpp>
#endif

#include <string>
#include <utility>


namespace impactx::elements::diagnostics::detail
{
#ifdef ImpactX_USE_OPENPMD
    /** Unclutter a real_names to openPMD record
     *
     * @TODO move to ABLASTR
     *
     * @param fullName name as in real_names variable
     * @return pair of openPMD record and component name
     */
    std::pair< std::string, std::string >
    name2openPMD (const std::string& fullName);


    //! TODO@ move to ablastr
    openPMD::RecordComponent get_component_record (
        openPMD::ParticleSpecies & species,
        std::string comp_name
    );

    /** Momentum of the reference particle in SI
     *
     * Particle momenta at fixed t are stored in units of it: it is their openPMD unitSI.
     *
     * @param ref reference particle
     * @return beta*gamma*m*c of the reference particle, in kg*m/s
     */
    double
    reference_momentum_SI (RefPart const & ref);

    /** Does the reference particle move along the z axis of the lab frame?
     *
     * Particles at fixed t are stored in the frame of the reference particle, whose z axis
     * points along its momentum. Only then is it also the lab frame, e.g., of WarpX.
     *
     * @param ref reference particle
     * @return true if its momentum points along +z, up to round-off
     */
    bool
    reference_along_z (RefPart const & ref);
#endif
} // namespace impactx::elements::diagnostics::detail

#endif // IMPACTX_OPENPMD_H
