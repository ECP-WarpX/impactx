/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_ELEMENT_H
#define IMPACTX_ELEMENTS_ELEMENT_H

#include "mixin/hostonly.H"
#include "mixin/named.H"

#include <optional>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>


namespace impactx::elements
{
    /** Common host-only metadata carried by every public lattice element.
     *
     * This is the single home for element state that describes the element to a human or
     * to tooling, as opposed to state that takes part in the physics: names today, and
     * descriptions, provenance, import annotations, tags or visualization data later.
     *
     * Adding a field here costs nothing in a kernel and requires no edit to any physics
     * implementation, because @see Element keeps this subobject out of the execution
     * snapshot that is copied into a launch closure.
     */
    struct ElementMetadata : mixin::Named
    {
        using mixin::Named::Named;
    };

/** Implementation details of the element envelope.
 *
 * Deliberately not called ``detail``: ``impactx::detail`` already exists and is used
 * unqualified from translation units that also do ``using namespace elements``
 * (e.g. initialization/InitElement.cpp), where a nested ``elements::detail`` would make
 * every such reference ambiguous.
 */
namespace envelope_detail
{
    /** True if @c T_Physics can be constructed from the whole argument pack */
    template<typename T_Physics, typename... T_Args>
    inline constexpr bool physics_from_all_v = std::is_constructible_v<T_Physics, T_Args...>;

    template<typename T_Physics, typename T_Tuple, typename T_Idx>
    struct physics_from_prefix : std::false_type
    {
    };

    template<typename T_Physics, typename... T_Args, std::size_t... T_I>
    struct physics_from_prefix<T_Physics, std::tuple<T_Args...>, std::index_sequence<T_I...>>
      : std::bool_constant<
            std::is_constructible_v<T_Physics, std::tuple_element_t<T_I, std::tuple<T_Args...>>...>
        >
    {
    };

    /** True if the pack is "physics arguments, then a trailing element name"
     *
     * The all-physics form is checked first and wins, so that an element whose physics
     * genuinely accepts a trailing string keeps today's overload resolution.
     *
     * This is a class template rather than one `&&` expression on purpose. `&&` in a
     * variable-template initializer short-circuits the *value*, not the instantiation of
     * the operands' types, so an empty pack would still form `sizeof...(T_Args) - 1` and
     * ask for a `std::make_index_sequence` of `SIZE_MAX`. Splitting the empty case into
     * the primary template keeps that expression from ever being formed.
     */
    template<typename T_Physics, typename... T_Args>
    struct physics_with_trailing_name : std::false_type
    {
    };

    template<typename T_Physics, typename T_First, typename... T_Rest>
    struct physics_with_trailing_name<T_Physics, T_First, T_Rest...>
      : std::bool_constant<
            !physics_from_all_v<T_Physics, T_First, T_Rest...> &&
            physics_from_prefix<
                T_Physics,
                std::tuple<T_First, T_Rest...>,
                std::make_index_sequence<sizeof...(T_Rest)>
            >::value &&
            std::is_convertible_v<
                std::tuple_element_t<sizeof...(T_Rest), std::tuple<T_First, T_Rest...>>,
                std::optional<std::string>
            >
        >
    {
    };

    template<typename T_Physics, typename... T_Args>
    inline constexpr bool physics_with_trailing_name_v =
        physics_with_trailing_name<T_Physics, T_Args...>::value;

    /** True if the pack is a single argument of (a type derived from) @c T_Self
     *
     * Used to keep the variadic constructors out of copy and move construction. Without
     * this, @c Element e2(e1) with a non-const lvalue @c e1 prefers the variadic template
     * over the copy constructor, because @c Element& converts to @c T_Physics const&
     * (derived-to-base) and so is an exact match for the template.
     *
     * The failure mode is not a compile error, which is why this guard is easy to lose
     * and important to keep: the variadic would copy-construct only the *physics* base
     * and default-construct the metadata, so the copy would silently come out unnamed.
     * See tests/python/test_element_metadata.py for the behavioral regression test.
     */
    template<typename T_Self, typename... T_Args>
    inline constexpr bool is_self_like_v =
        sizeof...(T_Args) == 1 &&
        (std::is_base_of_v<T_Self, std::remove_cvref_t<T_Args>> && ...);

} // namespace envelope_detail

    /** A public lattice element: element-specific physics plus common host metadata.
     *
     * ImpactX copies an element by value into the closure of a particle-push kernel. Any
     * host-only state reachable from that copy is either dead weight or, for owning state,
     * an outright hazard. This envelope makes the split explicit and generic:
     *
     * - @c T_Physics holds configuration, reference-particle and transfer-map algorithms,
     *   and the particle equations. It is what becomes an execution snapshot.
     * - @see ElementMetadata holds everything that only ever exists on the host.
     *
     * Because both are public bases, the scalar accessors and mixin properties that
     * existed before this split remain available directly on the public element, and no
     * per-element forwarding constructor or accessor is needed. Tracking and preparation
     * code nevertheless unwraps through @see physics explicitly, so the boundary is
     * visible at every call site rather than implied.
     *
     * @tparam T_Physics the element-specific physics implementation
     */
    template<typename T_Physics>
    class Element
      : public T_Physics,
        public ElementMetadata
    {
    public:
        using physics_type = T_Physics;

        /** Default-construct physics and metadata */
        Element () = default;

        /** Construct from physics arguments only; the element stays unnamed
         *
         * All arguments are forwarded to the physics constructor.
         */
        template<
            typename... T_Args,
            std::enable_if_t<
                envelope_detail::physics_from_all_v<T_Physics, T_Args...> &&
                !envelope_detail::is_self_like_v<Element, T_Args...>,
                int
            > = 0
        >
        explicit (sizeof...(T_Args) == 1)
        Element (T_Args&&... args)
          : T_Physics(std::forward<T_Args>(args)...)
        {
        }

        /** Construct from physics arguments plus a trailing element name
         *
         * The leading arguments are forwarded to the physics constructor and the last one
         * becomes the element name.
         */
        template<
            typename... T_Args,
            std::enable_if_t<
                envelope_detail::physics_with_trailing_name_v<T_Physics, T_Args...> &&
                !envelope_detail::is_self_like_v<Element, T_Args...>,
                int
            > = 0
        >
        explicit (sizeof...(T_Args) == 1)
        Element (T_Args&&... args)
          : Element(
                std::forward_as_tuple(std::forward<T_Args>(args)...),
                std::make_index_sequence<sizeof...(T_Args) - 1>{}
            )
        {
        }

        /** Access the physics subobject, i.e. everything that may enter a kernel */
        [[nodiscard]] T_Physics & physics () noexcept
        {
            return static_cast<T_Physics &>(*this);
        }

        /** Access the physics subobject, i.e. everything that may enter a kernel */
        [[nodiscard]] T_Physics const & physics () const noexcept
        {
            return static_cast<T_Physics const &>(*this);
        }

    private:
        /** Peel the trailing name off the argument pack
         *
         * @param args physics arguments, followed by the element name
         */
        template<typename T_Tuple, std::size_t... T_I>
        Element (T_Tuple&& args, std::index_sequence<T_I...>)
          : T_Physics(std::get<T_I>(std::forward<T_Tuple>(args))...),
            ElementMetadata(
                std::get<std::tuple_size_v<std::remove_reference_t<T_Tuple>> - 1>(
                    std::forward<T_Tuple>(args)
                )
            )
        {
        }
    };

    /** True if @c T is a public element envelope */
    template<typename T>
    inline constexpr bool is_element_envelope_v = false;

    template<typename T_Physics>
    inline constexpr bool is_element_envelope_v<Element<T_Physics>> = true;

    /** Unwrap a public element down to its physics.
     *
     * Use this at every point where an element crosses from host book-keeping into
     * reference tracking, transfer-map evaluation or pusher preparation, so that the
     * boundary is visible in the code rather than implied by an inherited member.
     *
     * Elements that carry no common metadata, and are therefore not wrapped
     * (@see diagnostics::BeamMonitor), are passed through unchanged.
     */
    template<typename T_Element>
    [[nodiscard]] decltype(auto) physics_of (T_Element && element) noexcept
    {
        if constexpr (is_element_envelope_v<std::remove_cvref_t<T_Element>>)
        {
            return std::forward<T_Element>(element).physics();
        }
        else
        {
            return static_cast<T_Element &&>(element);
        }
    }

} // namespace impactx::elements

#endif // IMPACTX_ELEMENTS_ELEMENT_H
