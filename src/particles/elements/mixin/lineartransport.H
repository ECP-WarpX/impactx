/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_MIXIN_LINEAR_TRANSPORT_H
#define IMPACTX_ELEMENTS_MIXIN_LINEAR_TRANSPORT_H

#include "particles/ImpactXParticleContainer.H"

#include <ablastr/constant.H>

#include <AMReX_Math.H>
#include <AMReX_Extension.H>
#include <AMReX_REAL.H>
#include <AMReX_SmallMatrix.H>


namespace impactx::elements
{
    /** This is a helper class for lattice elements that can be expressed as linear transport maps.
     */
    struct LinearTransport
    {
        /** ...
         */
        LinearTransport (
        )
        {
        }

        //LinearTransport () = default;
        LinearTransport (LinearTransport const &) = default;
        LinearTransport& operator= (LinearTransport const &) = default;
        LinearTransport (LinearTransport&&) = default;
        LinearTransport& operator= (LinearTransport&& rhs) = default;

        ~LinearTransport () = default;

        // 6x6 linear transport map
        using Map6x6 = amrex::SmallMatrix<amrex::ParticleReal, 6, 6, amrex::Order::F, 1>;
        Map6x6 m_transport_map; ///< linearized map
    };

} // namespace impactx::elements

#endif // IMPACTX_ELEMENTS_MIXIN_LINEAR_TRANSPORT_H
