/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Alex Bojanich, Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_PARTICLES_WAKEFIELDS_CSRBENDELEMENT_H_
#define IMPACTX_PARTICLES_WAKEFIELDS_CSRBENDELEMENT_H_

#include "particles/elements/Sbend.H"
#include "particles/elements/CFbend.H"

#include <AMReX_REAL.H>

#include <tuple>
#include <variant>


namespace impactx::particles::wakefields
{
    /** Function to calculate the radius of curvature (R) and check if an element has CSR
     *
     * @param[in] element_variant the lattice element to check
     * @returns boolean flag indicating if the element has CSR, radius of curvature (or zero)
     */
    template <typename VariantType>
    std::tuple<bool, amrex::Real>
    CSRBendElement (
        VariantType const& element_variant
    )
    {
        amrex::Real R = 0;
        bool element_has_csr = false;

        if (std::holds_alternative<Sbend>(element_variant)) // Check if variant holds a value of type Sbend
        {
            auto& element = std::get<Sbend>(element_variant); // Retrieves value stored in variant
            R = std::abs(element.m_rc);
            element_has_csr = true;
        }
        else if (std::holds_alternative<CFbend>(element_variant))
        {
            auto& element = std::get<CFbend>(element_variant); // Check if variant holds a value of type CFbend
            R = std::abs(element.m_rc);
            element_has_csr = true;
        }
        else
        {
            element_has_csr = false;
        }
        return {element_has_csr, R};
    }

} // namespace impactx::particles::wakefields

#endif // IMPACTX_PARTICLES_WAKEFIELDS_CSRBENDELEMENT_H_
