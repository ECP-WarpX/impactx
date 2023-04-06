/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Marco Garten, Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_BEAM_RELEVANT_H
#define IMPACTX_BEAM_RELEVANT_H

#include "particles/ImpactXParticleContainer.H"
#include "ImpactX.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuComplex.H>

#include <cmath>

using namespace amrex;

namespace impactx::diagnostics
{

    /** Compute momenta of the beam distribution and print them out
      */
    void compute_beam_relevant (ImpactXParticleContainer const & pc, int const step);

} // namespace impactx

#endif // BEAM_RELEVANT_H
