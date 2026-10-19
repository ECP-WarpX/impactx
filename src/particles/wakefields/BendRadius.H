/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Alex Bojanich, Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_PARTICLES_WAKEFIELDS_BENDRADIUS_H
#define IMPACTX_PARTICLES_WAKEFIELDS_BENDRADIUS_H

#include "elements/helper/Handle.H"
#include "elements/mixin/traits.H"
#include "particles/ReferenceParticle.H"

#include <AMReX_REAL.H>

#include <cmath>
#include <optional>


namespace impactx::particles::wakefields
{
    /** The bend radius of an element that emits synchrotron radiation
     *
     * CSR and ISR are applied only in these elements, @see elements::mixin::SynchrotronRadiation.
     *
     * @param[in] element_variant the lattice element to check
     * @param[in] refpart reference particle
     * @returns the magnitude of the radius of curvature of the reference orbit in m, or
     *          nothing if the element does not emit synchrotron radiation
     */
    template <typename T_Element>
    std::optional<amrex::Real>
    bend_radius (
        T_Element const & element_variant,
        RefPart const & refpart
    )
    {
        return elements::visit(
            [&refpart](auto const & element) -> std::optional<amrex::Real> {
                if constexpr (elements::mixin::has_synchrotron_radiation_v<decltype(element)>) {
                    return std::abs(element.rc(refpart));
                } else {
                    return std::nullopt;
                }
            },
            element_variant
        );
    }

} // namespace impactx::particles::wakefields

#endif // IMPACTX_PARTICLES_WAKEFIELDS_BENDRADIUS_H
