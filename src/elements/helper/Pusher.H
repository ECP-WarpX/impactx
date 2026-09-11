/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_HELPER_PUSHER_H
#define IMPACTX_ELEMENTS_HELPER_PUSHER_H

#include "elements/mixin/gpuelement.H"
#include "elements/mixin/hostonly.H"
#include "particles/ReferenceParticle.H"

#include <type_traits>
#include <utility>


namespace impactx::elements
{
// Deliberately not called `detail`: see the note in elements/mixin/gpuelement.H.
namespace pusher_detail
{
    template<typename T_GPU, typename = void>
    struct has_compute_constants : std::false_type {};

    template<typename T_GPU>
    struct has_compute_constants<
        T_GPU,
        std::void_t<decltype(
            std::declval<std::decay_t<T_GPU> &>().compute_constants(std::declval<RefPart const &>())
        )>
    > : std::true_type {};

    template<typename T_GPU, typename = void>
    struct has_own_pusher : std::false_type {};

    template<typename T_GPU>
    struct has_own_pusher<
        T_GPU,
        std::void_t<decltype(
            std::declval<std::decay_t<T_GPU> const &>().make_particle_pusher(std::declval<RefPart const &>())
        )>
    > : std::true_type {};

} // namespace pusher_detail

    /** Can a pusher be prepared for this GPU type?
     *
     * False for elements whose push is an opaque user callback (Programmable): there is
     * nothing to prepare, and the caller keeps using the element directly.
     */
    template<typename T_GPU>
    inline constexpr bool has_particle_pusher_v =
        pusher_detail::has_compute_constants<T_GPU>::value ||
        pusher_detail::has_own_pusher<T_GPU>::value;

    /** Compile-time checks that a type is fit to be copied into a launch closure.
     *
     * Called on every prepared pusher. The @see mixin::HostOnly check is the semantic
     * backstop: no single standard trait proves device-copyability across CUDA, HIP and
     * SYCL, so this combines the strongest portable checks with an explicit statement of
     * intent that a future host-only mixin cannot silently defeat.
     *
     * @tparam T_Pusher the prepared pusher type
     */
    template<typename T_Pusher>
    constexpr void assert_kernel_copyable ()
    {
        using pusher_type = std::decay_t<T_Pusher>;

        static_assert(
            !mixin::is_host_only_v<pusher_type>,
            "A pusher must not carry host-only state (mixin::HostOnly). "
            "Pass the GPU subobject, not the public element: element.gpu()."
        );
        static_assert(
            !is_element_envelope_v<pusher_type>,
            "A pusher must not be a public element envelope. "
            "Pass the GPU subobject instead: element.gpu()."
        );
        static_assert(
            std::is_trivially_destructible_v<pusher_type>,
            "A pusher must be trivially destructible: a launch closure is "
            "torn down without a host destructor running on the device copy."
        );
        static_assert(
            std::is_copy_constructible_v<pusher_type>,
            "A pusher must be copy-constructible: it is copied into the "
            "launch closure of every tile."
        );
    }

    /** Prepare the pusher that particle pushes run against.
     *
     * This is the single preparation boundary. It runs once per element slice, after the
     * reference particle has been advanced for that slice and before any tile is
     * launched, and the result is copied into each tile's launch closure. Particle-push
     * wrappers store the prepared value and execute; they do not prepare.
     *
     * Ordinary, resource-free GPU types are their own pusher: the returned value is
     * a prepared copy. The few GPU types that own dynamic array resources return a
     * separate, resource-free pusher instead, via their own @c make_particle_pusher member.
     *
     * @param gpu the canonical element GPU state; not modified
     * @param ref_part the reference particle, already advanced for this slice
     * @return a prepared, resource-free pusher
     */
    template<typename T_GPU>
    auto make_particle_pusher (T_GPU const & gpu, RefPart const & ref_part)
    {
        static_assert(
            !mixin::is_host_only_v<T_GPU>,
            "make_particle_pusher() takes an element's GPU state, not a type carrying "
            "host-only state. Use element.gpu()."
        );

        if constexpr (pusher_detail::has_own_pusher<T_GPU>::value)
        {
            auto pusher = gpu.make_particle_pusher(ref_part);
            assert_kernel_copyable<decltype(pusher)>();
            return pusher;
        }
        else
        {
            T_GPU pusher = gpu;
            pusher.compute_constants(ref_part);
            assert_kernel_copyable<T_GPU>();
            return pusher;
        }
    }

    /** The public element envelope is deliberately not preparable.
     *
     * It carries host-only metadata, which must never reach a kernel. Unwrap it first:
     * @c make_particle_pusher(element.gpu(), ref_part) .
     */
    template<typename T_GPU>
    auto make_particle_pusher (Element<T_GPU> const &, RefPart const &) = delete;

} // namespace impactx::elements

#endif // IMPACTX_ELEMENTS_HELPER_PUSHER_H
