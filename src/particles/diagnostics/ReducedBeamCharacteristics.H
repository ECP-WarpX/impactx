/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Marco Garten, Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_REDUCED_BEAM_CHARACTERISTICS_H
#define IMPACTX_REDUCED_BEAM_CHARACTERISTICS_H

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuComplex.H>

#include <cmath>

namespace impactx::diagnostics
{

    /** Compute momenta of the beam distribution and print them out
      */
    std::unordered_map<std::string, amrex::ParticleReal> reduced_beam_characteristics (ImpactXParticleContainer const & pc);

} // namespace impactx

#endif // IMPACTX_REDUCED_BEAM_CHARACTERISTICS
