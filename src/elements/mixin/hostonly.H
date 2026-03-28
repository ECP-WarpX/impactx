/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_MIXIN_HOSTONLY_H
#define IMPACTX_ELEMENTS_MIXIN_HOSTONLY_H

#include <type_traits>


namespace impactx::elements::mixin
{
    /** Marker for state that must never be copied into a device kernel.
     *
     * Derive from this in any mixin that owns host memory, holds a Python object, or is
     * otherwise meaningless on a device. Execution snapshots are checked against this
     * marker, so a future host-only mixin is rejected at compile time instead of silently
     * inflating every launch closure.
     *
     * This is an empty base: deriving it costs no storage.
     */
    struct HostOnly
    {
    };

    /** True if @c T carries host-only state */
    template<typename T>
    inline constexpr bool is_host_only_v = std::is_base_of_v<HostOnly, std::remove_cvref_t<T>>;

} // namespace impactx::elements::mixin

#endif // IMPACTX_ELEMENTS_MIXIN_HOSTONLY_H
