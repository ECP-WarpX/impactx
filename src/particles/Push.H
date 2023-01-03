/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_PUSH_H
#define IMPACTX_PUSH_H

#include "elements/All.H"
#include "particles/ImpactXParticleContainer.H"

#include <list>


namespace impactx
{
    /** Push particles
     *
     * @param pc container of the particles to push
     * @param element_variant a single element to push the particles through
     */
    void Push (ImpactXParticleContainer & pc,
               KnownElements const & element_variant);

} // namespace impactx

#endif // IMPACTX_PUSH_H
