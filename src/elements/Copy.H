/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_COPY_H
#define IMPACTX_ELEMENTS_COPY_H

#include "Handle.H"

#include <memory>
#include <type_traits>
#include <utility>
#include <variant>


namespace impactx::elements
{
namespace copy_detail
{
    template<typename T_Element, typename = void>
    struct has_member_copy : std::false_type {};

    template<typename T_Element>
    struct has_member_copy<
        T_Element,
        std::void_t<decltype(std::declval<T_Element const &>().copy())>
    > : std::true_type {};

} // namespace copy_detail

    /** A distinct element with the same configuration.
     *
     * Most elements are pure configuration, so a value copy is exactly right -- including
     * the ones holding coefficient arrays, since those deep-copy.
     *
     * An element that also holds state which two elements must not share provides its own
     * ``copy()``, and that is used instead. @see diagnostics::BeamMonitor, whose copy must
     * not inherit an already-open output series.
     */
    template<typename T_Element>
    [[nodiscard]] T_Element copy_of (T_Element const & element)
    {
        if constexpr (copy_detail::has_member_copy<T_Element>::value)
        {
            return element.copy();
        }
        else
        {
            return element;
        }
    }

    /** A handle owning a distinct element with the same configuration
     *
     * This is the operation behind the user-facing ``copy()``: the result is a separate
     * physical element, so putting it in a lattice adds an element rather than another
     * occurrence of an existing one.
     */
    [[nodiscard]] inline ElementHandle copy (ElementHandle const & handle)
    {
        require_valid(handle, "copy");
        return std::visit(
            [](auto const & el) {
                using element_type = std::decay_t<decltype(*el)>;
                return ElementHandle{std::make_shared<element_type>(copy_of(*el))};
            },
            handle
        );
    }

} // namespace impactx::elements

#endif // IMPACTX_ELEMENTS_COPY_H
