/* Copyright 2022 The Regents of the University of California, through Lawrence
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


namespace impactx
{
namespace transformation
{
    /** Direction of the \see CoordinateTransformation
     */
    enum class Direction
    {
        T2Z,
        Z2T
    };

    /** Transform the coordinate system of all particles
     *
     * @param pc container of the particles to push
     * @param direction the direction (Z to T or vice versa)
     */
    void CoordinateTransformation (ImpactXParticleContainer &pc,
                                   Direction const & direction);

} // namespace transformation
} // namespace impactx

#endif // IMPACTX_COORDINATE_TRANSFORMATION_H
