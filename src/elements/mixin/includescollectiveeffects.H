/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_MIXIN_INCLUDES_COLLECTIVE_EFFECTS_H
#define IMPACTX_ELEMENTS_MIXIN_INCLUDES_COLLECTIVE_EFFECTS_H

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>


namespace impactx::elements::mixin
{
    /** This is a helper class for lattice elements whose push can model collective effects
     *
     * An element with a user-supplied push, e.g., one that calls an external code which
     * computes the beam self-fields or plasma wakefields, declares with this flag that its
     * push already models collective effects over its length. ImpactX then applies no
     * additional collective effect kicks (space charge, wakefields, CSR, ISR) for this element.
     */
    struct IncludesCollectiveEffects
    {
        static constexpr bool DEFAULT_includes_collective_effects = false;

        /** An element that declares if its push already models collective effects
         *
         * @param includes_collective_effects the push already models collective effects
         */
        explicit IncludesCollectiveEffects (
            bool includes_collective_effects = DEFAULT_includes_collective_effects
        )
          : m_includes_collective_effects(includes_collective_effects)
        {
        }

        /** Does the push of this element already model collective effects?
         *
         * @return true if ImpactX applies no additional collective effect kicks for this element
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        bool includes_collective_effects () const
        {
            return m_includes_collective_effects;
        }

        /** The push already models collective effects over the length of this element
         *
         * The default value is false: ImpactX applies the enabled collective effects to this
         * element like to any other element of finite length.
         */
        bool m_includes_collective_effects = DEFAULT_includes_collective_effects;
    };

} // namespace impactx::elements::mixin

#endif // IMPACTX_ELEMENTS_MIXIN_INCLUDES_COLLECTIVE_EFFECTS_H
