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

#include <AMReX_REAL.H>
#include <vector>
#include "particles/ImpactXParticleContainer.H"

namespace impactx::particles::wakefields
{
    //Function to update particle momentum using wake force kick
    void WakePush(ImpactXParticleContainer & pc, const std::vector<amrex::Real>& convoluted_wakefield, amrex::ParticleReal const slice_ds, amrex::Real bin_size, amrex::Real bin_min, int padding_factor);
}

#endif // WAKEPUSH_H
