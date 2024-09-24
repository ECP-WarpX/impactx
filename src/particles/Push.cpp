/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "Push.H"

#include <AMReX_BLProfiler.H>

#include <variant>


namespace impactx
{
    void Push (
        ImpactXParticleContainer & pc,
        KnownElements & element_variant,
        int step,
        int cycle
    )
    {
        // here we just access the element by its respective type
        std::visit([&pc, step, cycle](auto&& element)
        {
            BL_PROFILE("impactx::Push");

            // push reference particle & all particles
            element(pc, step, cycle);
        }, element_variant);
    }

} // namespace impactx
