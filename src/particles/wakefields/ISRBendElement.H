/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Alex Bojanich, Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_PARTICLES_WAKEFIELDS_ISRBENDELEMENT_H
#define IMPACTX_PARTICLES_WAKEFIELDS_ISRBENDELEMENT_H

#include "elements/ExactSbend.H"
#include "elements/Sbend.H"
#include "elements/CFbend.H"
#include "particles/ReferenceParticle.H"

#include <AMReX_REAL.H>

#include <cmath>
#include <tuple>
#include <variant>


namespace impactx::particles::wakefields
{
    /** Function to calculate the radius of curvature (R) and check if an element has ISR
     *
     * @param[in] element_variant the lattice element to check
     * @param[in] refpart reference particle
     * @returns boolean flag indicating if the element has ISR, radius of curvature (or zero)
     */
    template <typename VariantType>
    std::tuple<bool, amrex::Real>
    ISRBendElement (
        VariantType const& element_variant,
        RefPart const & refpart
    )
    {
        amrex::Real R = 0;
        bool element_has_isr = false;

        if (auto const * sbend = elements::get_if<elements::Sbend>(element_variant))
        {
            R = std::abs(sbend->m_rc);
            element_has_isr = true;
        }
        else if (auto const * cfbend = elements::get_if<elements::CFbend>(element_variant))
        {
            R = std::abs(cfbend->m_rc);
            element_has_isr = true;
        }
        else if (auto const * exact_sbend = elements::get_if<elements::ExactSbend>(element_variant))
        {
            R = std::abs(exact_sbend->rc(refpart));
            element_has_isr = true;
        }

        return {element_has_isr, R};
    }

} // namespace impactx::particles::wakefields

#endif // IMPACTX_PARTICLES_WAKEFIELDS_ISRBENDELEMENT_H
