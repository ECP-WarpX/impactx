/* Copyright 2021 Axel Huebl, Ji Qiang
 *
 * This file is part of ImpactX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_LIST_H
#define IMPACTX_ELEMENTS_LIST_H

#include "Drift.H"
#include "Sbend.H"
#include "Quad.H"

#include <variant>


namespace impactx
{
    using KnownElements = std::variant<Drift, Sbend, Quad>;

} // namespace impactx

#endif // IMPACTX_ELEMENTS_LIST_H
