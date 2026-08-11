/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_MIXIN_GPUELEMENT_H
#define IMPACTX_ELEMENTS_MIXIN_GPUELEMENT_H

#include "hostonly.H"
#include "named.H"

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
     * to tooling, as opposed to state that enters a kernel: names today, and
     * descriptions, provenance, import annotations, tags or visualization data later.
     *
     * Adding a field here costs nothing in a kernel and requires no edit to any GPU
     * implementation, because @see Element keeps this subobject out of the execution
     * pusher that is copied into a launch closure.
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
    /** True if @c T_GPU can be constructed from the whole argument pack */
    template<typename T_GPU, typename... T_Args>
    inline constexpr bool gpu_from_all_v = std::is_constructible_v<T_GPU, T_Args...>;

    template<typename T_GPU, typename T_Tuple, typename T_Idx>
    struct gpu_from_prefix : std::false_type
    {
    };

    template<typename T_GPU, typename... T_Args, std::size_t... T_I>
    struct gpu_from_prefix<T_GPU, std::tuple<T_Args...>, std::index_sequence<T_I...>>
      : std::bool_constant<
            std::is_constructible_v<T_GPU, std::tuple_element_t<T_I, std::tuple<T_Args...>>...>
        >
    {
    };

    /** True if the pack is "GPU arguments, then a trailing element name"
     *
     * The all-GPU form is checked first and wins, so that an element whose GPU type
     * genuinely accepts a trailing string keeps today's overload resolution.
     *
     * This is a class template rather than one `&&` expression on purpose. `&&` in a
     * variable-template initializer short-circuits the *value*, not the instantiation of
     * the operands' types, so an empty pack would still form `sizeof...(T_Args) - 1` and
     * ask for a `std::make_index_sequence` of `SIZE_MAX`. Splitting the empty case into
     * the primary template keeps that expression from ever being formed.
     */
    template<typename T_GPU, typename... T_Args>
    struct gpu_with_trailing_name : std::false_type
    {
    };

    template<typename T_GPU, typename T_First, typename... T_Rest>
    struct gpu_with_trailing_name<T_GPU, T_First, T_Rest...>
      : std::bool_constant<
            !gpu_from_all_v<T_GPU, T_First, T_Rest...> &&
            gpu_from_prefix<
                T_GPU,
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

    template<typename T_GPU, typename... T_Args>
    inline constexpr bool gpu_with_trailing_name_v =
        gpu_with_trailing_name<T_GPU, T_Args...>::value;

    /** True if the pack is a single argument of (a type derived from) @c T_Self
     *
     * Used to keep the variadic constructors out of copy and move construction. Without
     * this, @c Element e2(e1) with a non-const lvalue @c e1 prefers the variadic template
     * over the copy constructor, because @c Element& converts to @c T_GPU const&
     * (derived-to-base) and so is an exact match for the template.
     *
     * The failure mode is not a compile error, which is why this guard is easy to lose
     * and important to keep: the variadic would copy-construct only the *GPU* base
     * and default-construct the metadata, so the copy would silently come out unnamed.
     * See tests/python/test_element_metadata.py for the behavioral regression test.
     */
    template<typename T_Self, typename... T_Args>
    inline constexpr bool is_self_like_v =
        sizeof...(T_Args) == 1 &&
        (std::is_base_of_v<T_Self, std::remove_cvref_t<T_Args>> && ...);

} // namespace envelope_detail

    /** A public lattice element: element-specific GPU state plus common host metadata.
     *
     * ImpactX copies an element by value into the closure of a particle-push kernel. Any
     * host-only state reachable from that copy is either dead weight or, for owning state,
     * an outright hazard. This envelope makes the split explicit and generic:
     *
     * - @c T_GPU holds configuration, reference-particle and transfer-map algorithms,
     *   and the particle equations. It is what becomes a pusher.
     * - @see ElementMetadata holds everything that only ever exists on the host.
     *
     * Because both are public bases, the scalar accessors, mixin properties and host
     * entry points -- the particle-container push, the reference-particle push, the
     * covariance-matrix transport and @c transport_map -- remain callable directly on the
     * public element, and no per-element forwarding is needed. Those entry points reach
     * the GPU subobject on their own: each casts @c this to @c T_GPU, so metadata cannot
     * follow them into a kernel. Only pusher preparation takes the subobject
     * explicitly, through @see gpu, and @see make_particle_pusher rejects the envelope.
     *
     * @tparam T_GPU the element-specific GPU implementation
     */
    template<typename T_GPU>
    class Element
      : public T_GPU,
        public ElementMetadata
    {
    public:
        using gpu_type = T_GPU;

        /** Default-construct GPU state and metadata */
        Element () = default;

        /** Construct from GPU arguments only; the element stays unnamed
         *
         * All arguments are forwarded to the GPU constructor.
         */
        template<
            typename... T_Args,
            std::enable_if_t<
                envelope_detail::gpu_from_all_v<T_GPU, T_Args...> &&
                !envelope_detail::is_self_like_v<Element, T_Args...>,
                int
            > = 0
        >
        explicit (sizeof...(T_Args) == 1)
        Element (T_Args&&... args)
          : T_GPU(std::forward<T_Args>(args)...)
        {
        }

        /** Construct from GPU arguments plus a trailing element name
         *
         * The leading arguments are forwarded to the GPU constructor and the last one
         * becomes the element name.
         */
        template<
            typename... T_Args,
            std::enable_if_t<
                envelope_detail::gpu_with_trailing_name_v<T_GPU, T_Args...> &&
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

        /** Access the GPU subobject, i.e. everything that may enter a kernel */
        [[nodiscard]] T_GPU & gpu () noexcept
        {
            return static_cast<T_GPU &>(*this);
        }

        /** Access the GPU subobject, i.e. everything that may enter a kernel */
        [[nodiscard]] T_GPU const & gpu () const noexcept
        {
            return static_cast<T_GPU const &>(*this);
        }

    private:
        /** Peel the trailing name off the argument pack
         *
         * @param args GPU arguments, followed by the element name
         */
        template<typename T_Tuple, std::size_t... T_I>
        Element (T_Tuple&& args, std::index_sequence<T_I...>)
          : T_GPU(std::get<T_I>(std::forward<T_Tuple>(args))...),
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

    template<typename T_GPU>
    inline constexpr bool is_element_envelope_v<Element<T_GPU>> = true;

} // namespace impactx::elements

#endif // IMPACTX_ELEMENTS_MIXIN_GPUELEMENT_H
