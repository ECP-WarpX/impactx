/* Copyright 2021 Axel Huebl
 *
 * This file is part of ImpactX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_PUSH_H
#define IMPACTX_PUSH_H

#include "particles/ImpactXParticleContainer.H"


namespace impactx
{
    /** Push particles
     *
     * @param pc container of the particles to push
     */
    void Push (ImpactXParticleContainer & pc);

} // namespace impactx

#endif // IMPACTX_PUSH_H
