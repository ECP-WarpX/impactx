/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "Push.H"
#include "elements/mixin/gpuelement.H"

#include <AMReX_BLProfiler.H>

#include <variant>


namespace impactx
{
    void push (
        ImpactXParticleContainer & pc,
        elements::ElementHandle & element_variant,
        int step,
        int period
    )
    {
        // here we just access the element by its respective type
        elements::visit([&pc, step, period](auto&& element)
        {
            BL_PROFILE("impactx::push");

            // push reference particle & all particles
            element(pc, step, period);
        }, element_variant);
    }

    void push (
        RefPart & ref,
        elements::ElementHandle & element_variant
    )
    {
        // here we just access the element by its respective type
        elements::visit([&ref](auto&& element)
        {
            // push reference particle in global coordinates
            BL_PROFILE("impactx::push::RefPart");
            element(ref);
        }, element_variant);
    }

} // namespace impactx
