/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Eric G. Stern, Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "PolygonAperture.H"

#include <stdexcept>
#include <string>

std::string
impactx::elements::PolygonApertureCore::action_name (Action const & action)
{
    return amrex::getEnumNameString(action);
}

IMPACTX_PUSH_INSTANTIATE(impactx::elements::PolygonAperturePhysics)
