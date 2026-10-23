/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Alex Bojanich, Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_PARTICLES_WAKEFIELDS_SHORTRANGEWAKEELEMENT_H_
#define IMPACTX_PARTICLES_WAKEFIELDS_SHORTRANGEWAKEELEMENT_H_

#include "elements/helper/Handle.H"
#include "elements/mixin/shortrangewake.H"
#include "elements/mixin/traits.H"

#include <optional>


namespace impactx::particles::wakefields
{
    /** The user-specified short-range structure wake of an element
     *
     * Short-range structure wakes are applied only in elements that carry one,
     * @see elements::mixin::ShortRangeWake.
     *
     * @param[in] element_variant the lattice element to check
     * @returns a view on the element's wake function data, or nothing if the element does not
     *          contribute a short-range wake
     */
    template <typename T_Element>
    std::optional<elements::mixin::ShortRangeWakeView>
    ShortRangeWakeElement (
        T_Element const & element_variant
    )
    {
        return elements::visit(
            [](auto const & element) -> std::optional<elements::mixin::ShortRangeWakeView> {
                if constexpr (elements::mixin::has_short_range_wake_v<decltype(element)>) {
                    if (element.has_short_range_wake()) {
                        return element.short_range_wake();
                    }
                }
                return std::nullopt;
            },
            element_variant
        );
    }

} // namespace impactx::particles::wakefields

#endif // IMPACTX_PARTICLES_WAKEFIELDS_SHORTRANGEWAKEELEMENT_H_
