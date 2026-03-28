/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_MIXIN_SNAPSHOT_H
#define IMPACTX_ELEMENTS_MIXIN_SNAPSHOT_H

#include "elements/Element.H"
#include "hostonly.H"
#include "particles/ReferenceParticle.H"

#include <type_traits>
#include <utility>


namespace impactx::elements
{
// Deliberately not called `detail`: see the note in elements/Element.H.
namespace snapshot_detail
{
    template<typename T_Physics, typename = void>
    struct has_compute_constants : std::false_type {};

    template<typename T_Physics>
    struct has_compute_constants<
        T_Physics,
        std::void_t<decltype(
            std::declval<std::decay_t<T_Physics> &>().compute_constants(std::declval<RefPart const &>())
        )>
    > : std::true_type {};

    template<typename T_Physics, typename = void>
    struct has_own_pusher : std::false_type {};

    template<typename T_Physics>
    struct has_own_pusher<
        T_Physics,
        std::void_t<decltype(
            std::declval<std::decay_t<T_Physics> const &>().make_particle_pusher(std::declval<RefPart const &>())
        )>
    > : std::true_type {};

} // namespace snapshot_detail

    /** Can an execution snapshot be prepared for this physics type?
     *
     * False for elements whose push is an opaque user callback (Programmable): there is
     * nothing to prepare, and the caller keeps using the element directly.
     */
    template<typename T_Physics>
    inline constexpr bool has_particle_pusher_v =
        snapshot_detail::has_compute_constants<T_Physics>::value ||
        snapshot_detail::has_own_pusher<T_Physics>::value;

    /** Compile-time checks that a type is fit to be copied into a launch closure.
     *
     * Called on every prepared pusher. The @see mixin::HostOnly check is the semantic
     * backstop: no single standard trait proves device-copyability across CUDA, HIP and
     * SYCL, so this combines the strongest portable checks with an explicit statement of
     * intent that a future host-only mixin cannot silently defeat.
     *
     * @tparam T_Snapshot the prepared pusher type
     */
    template<typename T_Snapshot>
    constexpr void assert_execution_snapshot ()
    {
        using snapshot_type = std::decay_t<T_Snapshot>;

        static_assert(
            !mixin::is_host_only_v<snapshot_type>,
            "An execution snapshot must not carry host-only state (mixin::HostOnly). "
            "Pass the physics subobject, not the public element: element.physics()."
        );
        static_assert(
            !is_element_envelope_v<snapshot_type>,
            "An execution snapshot must not be a public element envelope. "
            "Pass the physics subobject instead: element.physics()."
        );
        static_assert(
            std::is_trivially_destructible_v<snapshot_type>,
            "An execution snapshot must be trivially destructible: a launch closure is "
            "torn down without a host destructor running on the device copy."
        );
        static_assert(
            std::is_copy_constructible_v<snapshot_type>,
            "An execution snapshot must be copy-constructible: it is copied into the "
            "launch closure of every tile."
        );
    }

    /** Prepare the execution snapshot that particle pushes run against.
     *
     * This is the single preparation boundary. It runs once per element slice, after the
     * reference particle has been advanced for that slice and before any tile is
     * launched, and the result is copied into each tile's launch closure. Particle-push
     * wrappers store the prepared value and execute; they do not prepare.
     *
     * Ordinary, resource-free physics types are their own snapshot: the returned value is
     * a prepared copy. The few physics types that own dynamic array resources return a
     * separate, resource-free core instead, via their own @c make_particle_pusher member.
     *
     * @param physics the canonical element physics; not modified
     * @param ref_part the reference particle, already advanced for this slice
     * @return a prepared, resource-free snapshot
     */
    template<typename T_Physics>
    auto make_particle_pusher (T_Physics const & physics, RefPart const & ref_part)
    {
        static_assert(
            !mixin::is_host_only_v<T_Physics>,
            "make_particle_pusher() takes an element's physics, not a type carrying "
            "host-only state. Use element.physics()."
        );

        if constexpr (snapshot_detail::has_own_pusher<T_Physics>::value)
        {
            auto pusher = physics.make_particle_pusher(ref_part);
            assert_execution_snapshot<decltype(pusher)>();
            return pusher;
        }
        else
        {
            T_Physics pusher = physics;
            pusher.compute_constants(ref_part);
            assert_execution_snapshot<T_Physics>();
            return pusher;
        }
    }

    /** The public element envelope is deliberately not preparable.
     *
     * It carries host-only metadata, which must never reach a kernel. Unwrap it first:
     * @c make_particle_pusher(element.physics(), ref_part) .
     */
    template<typename T_Physics>
    auto make_particle_pusher (Element<T_Physics> const &, RefPart const &) = delete;

    /** The prepared pusher type for a given physics type */
    template<typename T_Physics>
    using particle_pusher_t = decltype(
        make_particle_pusher(std::declval<T_Physics const &>(), std::declval<RefPart const &>())
    );

} // namespace impactx::elements

#endif // IMPACTX_ELEMENTS_MIXIN_SNAPSHOT_H
