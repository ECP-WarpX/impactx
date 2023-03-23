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

#include <string>
#include <variant>


namespace impactx
{
    void Push (ImpactXParticleContainer & pc,
               KnownElements & element_variant,
               int step)
    {
        // here we just access the element by its respective type
        std::visit([&pc, step](auto&& element)
        {
            // performance profiling per element
            std::string element_name;
            element_name = element.name;
            std::string const profile_name = "impactx::Push::" + element_name;
            BL_PROFILE("impactx::Push");
            BL_PROFILE(profile_name);

            // push reference particle & all particles
            element(pc, step);
        }, element_variant);
    }

} // namespace impactx
