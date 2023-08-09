/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_COLLECT_LOST_H
#define IMPACTX_COLLECT_LOST_H

#include "particles/ImpactXParticleContainer.H"


namespace impactx
{
    /** ...
     *
     * ...
     *
     * @param source
     * @param dest
     */
    void collect_lost_particles (ImpactXParticleContainer& source, ImpactXParticleContainer& dest);

} // namespace impactx

#endif // IMPACTX_COLLECT_LOST_H
