/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_HANDLE_H
#define IMPACTX_ELEMENTS_HANDLE_H

#include "All.H"

#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <variant>


namespace impactx::elements
{
    /** Turn a variant of values into a variant of shared pointers to those values */
    template<typename T_Variant>
    struct SharedVariant;

    template<typename... T_Elements>
    struct SharedVariant<std::variant<T_Elements...>>
    {
        using type = std::variant<std::shared_ptr<T_Elements>...>;
    };

    /** Shared ownership of one concrete lattice element.
     *
     * Each alternative is an ordinary @c std::shared_ptr to a concrete element type, so
     * the pointee keeps its exact type: a visitor still sees a @c Quad, not a base class,
     * and no downcast is involved.
     *
     * Sharing is what lets the same element appear at several lattice positions, be held
     * by more than one lattice, and stay alive while something outside the lattice -- a
     * Python variable, a tracking hook -- still refers to it.
     *
     * A default-constructed handle holds a null pointer in its first alternative. That is
     * never a valid lattice entry; reject it at every boundary that stores or visits one,
     * @see valid.
     */
    using ElementHandle = SharedVariant<KnownElements>::type;

    /** True if @c T is a way of denoting one lattice element
     *
     * Either an owning @see ElementHandle or a plain @c KnownElements value. Accessors are
     * written once against this, so a call site does not care which it holds.
     */
    template<typename T>
    inline constexpr bool is_element_repr_v =
        std::is_same_v<std::remove_cvref_t<T>, KnownElements> ||
        std::is_same_v<std::remove_cvref_t<T>, ElementHandle>;

    /** True if this handle owns an element */
    [[nodiscard]] inline bool valid (ElementHandle const & handle) noexcept
    {
        return std::visit([](auto const & el) { return static_cast<bool>(el); }, handle);
    }

    /** Throw unless this handle owns an element
     *
     * @param handle the handle to check
     * @param context what was attempted, used in the error message
     */
    inline void require_valid (ElementHandle const & handle, char const * context)
    {
        if (!valid(handle))
        {
            throw std::runtime_error(
                std::string("Lattice: ") + context + " with an empty element handle!");
        }
    }

    /** Build a handle owning a freshly constructed element
     *
     * @tparam T_Element the concrete element type
     * @param args arguments forwarded to the element constructor
     */
    template<typename T_Element, typename... T_Args>
    [[nodiscard]] ElementHandle make_handle (T_Args&&... args)
    {
        return ElementHandle{std::make_shared<T_Element>(std::forward<T_Args>(args)...)};
    }

    /** Build a handle owning a copy of an element value */
    [[nodiscard]] inline ElementHandle make_handle (KnownElements const & element)
    {
        return std::visit(
            [](auto const & el) {
                return ElementHandle{std::make_shared<std::decay_t<decltype(el)>>(el)};
            },
            element
        );
    }

    /** Apply a visitor to the element a handle owns
     *
     * The visitor receives the concrete element by reference, exactly as it would for a
     * @c KnownElements value, so the same visitors serve both representations.
     *
     * @param visitor the callable to apply
     * @param handle a non-null element handle
     */
    template<typename T_Visitor>
    decltype(auto) visit (T_Visitor && visitor, ElementHandle & handle)
    {
        require_valid(handle, "visit");
        return std::visit(
            [&](auto & el) -> decltype(auto) { return std::forward<T_Visitor>(visitor)(*el); },
            handle
        );
    }

    /** @see visit
     *
     * A const handle protects the element it owns. Dereferencing the shared pointer alone
     * would not say so: it hands out a mutable element however const the pointer is.
     */
    template<typename T_Visitor>
    decltype(auto) visit (T_Visitor && visitor, ElementHandle const & handle)
    {
        require_valid(handle, "visit");
        return std::visit(
            [&](auto const & el) -> decltype(auto) {
                return std::forward<T_Visitor>(visitor)(std::as_const(*el));
            },
            handle
        );
    }

    /** Apply a visitor to an element value
     *
     * Provided so that call sites can be written once against either representation.
     */
    template<typename T_Visitor>
    decltype(auto) visit (T_Visitor && visitor, KnownElements & element)
    {
        return std::visit(std::forward<T_Visitor>(visitor), element);
    }

    /** @see visit */
    template<typename T_Visitor>
    decltype(auto) visit (T_Visitor && visitor, KnownElements const & element)
    {
        return std::visit(std::forward<T_Visitor>(visitor), element);
    }

    /** Access the element if it is of the requested type, else nullptr
     *
     * The counterpart of @c std::get_if for handles: @c std::holds_alternative and
     * @c std::get do not work here, because the alternatives are shared pointers rather
     * than the element types themselves.
     *
     * @tparam T_Element the concrete element type to look for
     */
    template<typename T_Element>
    [[nodiscard]] T_Element * get_if (ElementHandle & handle) noexcept
    {
        auto * const held = std::get_if<std::shared_ptr<T_Element>>(&handle);
        return held ? held->get() : nullptr;
    }

    /** @see get_if */
    template<typename T_Element>
    [[nodiscard]] T_Element const * get_if (ElementHandle const & handle) noexcept
    {
        auto const * const held = std::get_if<std::shared_ptr<T_Element>>(&handle);
        return held ? held->get() : nullptr;
    }

    /** The address of the owned element, as an identity key
     *
     * Two handles denote the same physical element exactly when this matches.
     */
    [[nodiscard]] inline void const * address (ElementHandle const & handle) noexcept
    {
        return std::visit(
            [](auto const & el) { return static_cast<void const *>(el.get()); },
            handle
        );
    }

    /** Copy the owned element into a fresh, independently owned handle
     *
     * A plain value copy of the concrete element. Use it where a value copy is known to be
     * right -- splitting a thick beam-optic element, say.
     *
     * It is not the user-facing ``copy()``, which additionally resets state that must not
     * be shared between two elements, such as a beam monitor's open output series.
     */
    [[nodiscard]] inline ElementHandle copy_element (ElementHandle const & handle)
    {
        require_valid(handle, "copy");
        return std::visit(
            [](auto const & el) {
                return ElementHandle{std::make_shared<std::decay_t<decltype(*el)>>(*el)};
            },
            handle
        );
    }

} // namespace impactx::elements

#endif // IMPACTX_ELEMENTS_HANDLE_H
