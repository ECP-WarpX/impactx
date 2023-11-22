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
    /** Move lost particles into a separate container
     *
     * If particles are marked as lost, by setting their id to negative, we
     * will move them to another particle container, store their position when
     * lost and stop pushing them in the beamline.
     *
     * @param source the beam particle container that might loose particles
     */
    void collect_lost_particles (ImpactXParticleContainer& source);

} // namespace impactx

#endif // IMPACTX_COLLECT_LOST_H
