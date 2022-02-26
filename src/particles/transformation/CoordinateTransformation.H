/* Copyright 2021 Axel Huebl, Chad Mitchell
 *
 * This file is part of ImpactX.
 *
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
    void CoordinateTransformation(ImpactXParticleContainer &pc,
                                  Direction const & direction);

} // namespace transformation
} // namespace impactx

#endif // IMPACTX_COORDINATE_TRANSFORMATION_H
