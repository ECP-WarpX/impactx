/* Copyright 2021 Axel Huebl
 *
 * This file is part of ImpactX.
 *
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
     * @param lattice elements to push the particles through
     */
    void Push (ImpactXParticleContainer & pc,
               std::list<KnownElements> const & lattice);

} // namespace impactx

#endif // IMPACTX_PUSH_H
