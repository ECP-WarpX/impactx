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
     * @param[inout] pc container of the particles to push
     * @param[inout] element_variant a single element to push the particles through
     * @param[in] step global step for diagnostics
     * @param[in] cycle for periodic lattices, this is the cycle or turn
     */
    void Push (
        ImpactXParticleContainer & pc,
        KnownElements & element_variant,
        int step,
        int cycle
    );

} // namespace impactx

#endif // IMPACTX_PUSH_H
