/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_PARTICLES_WAKEFIELDS_INTERPOLATE_LINEAR_1D_H
#define IMPACTX_PARTICLES_WAKEFIELDS_INTERPOLATE_LINEAR_1D_H

#include <AMReX_Algorithm.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>


namespace impactx::particles::wakefields
{
    /** Interpolate equidistant 1D values linearly to a position
     *
     * Value k of values belongs to the position x_first + k / inv_dx. Positions outside of the
     * values get the nearest end value.
     *
     * @param[in] x position to interpolate to
     * @param[in] values equidistant values, at least two
     * @param[in] num_values number of values
     * @param[in] x_first position of the first value
     * @param[in] inv_dx inverse of the distance between two values
     * @return value at x
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real InterpolateLinear1D (
        amrex::Real x,
        amrex::Real const * values,
        int num_values,
        amrex::Real x_first,
        amrex::Real inv_dx
    )
    {
        using namespace amrex::literals;

        // between the values k and k+1, at the weight w in [0, 1] from value k
        amrex::Real const v = (x - x_first) * inv_dx;
        int const k = amrex::Clamp(static_cast<int>(amrex::Math::floor(v)), 0, num_values - 2);
        amrex::Real const w = amrex::Clamp(v - static_cast<amrex::Real>(k), 0_rt, 1_rt);

        return (1_rt - w) * values[k] + w * values[k + 1];
    }

} // namespace impactx::particles::wakefields

#endif // IMPACTX_PARTICLES_WAKEFIELDS_INTERPOLATE_LINEAR_1D_H
