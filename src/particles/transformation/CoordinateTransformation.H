/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_COORDINATE_TRANSFORMATION_H
#define IMPACTX_COORDINATE_TRANSFORMATION_H

#include "particles/ImpactXParticleContainer.H"


namespace impactx::transformation
{
    /** Transform the coordinate system of all particles
     *
     * @param pc container of the particles to push
     * @param direction the transformation (to fixed s or to fixed t)
     */
    void CoordinateTransformation (ImpactXParticleContainer & pc,
                                   CoordSystem direction);

} // namespace impactx::transformation

#endif // IMPACTX_COORDINATE_TRANSFORMATION_H
