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

#include <variant>
#include <AMReX_REAL.H>
#include "particles/elements/Sbend.H"
 #include "particles/elements/CFbend.H"

namespace impactx::particles::wakefields
{

    /** Function to calculate the radius of curvature (R) and check if an element has CSR
     *
     * @param[in,out] R radius of curvature (updated if element has CSR)
     * @param[in,out] element_has_csr boolean flag indicating if the element has CSR (updated if true)
     * @param[in] element the lattice element to check
     */
    template <typename VariantType>
    void CSRBendElement(amrex::Real& R, bool& element_has_csr, VariantType const& element_variant)
    {
        if (std::holds_alternative<Sbend>(element_variant)) //Check if variant holds a value of type Sbend
        {
            auto& element = std::get<Sbend>(element_variant); //Retrieves value stored in variant
            R = element.m_rc;
            R = std::abs(R);
            element_has_csr = true;
        }
        else if (std::holds_alternative<CFbend>(element_variant))
        {
            auto& element = std::get<CFbend>(element_variant); //Check if variant holds a value of type CFbend
            R = element.m_rc;
            R = std::abs(R);
            element_has_csr = true;
        }
        else
        {
            element_has_csr = false;
        }
    }

} // namespace impactx::particles::wakefields

#endif // IMPACTX_PARTICLES_WAKEFIELDS_CSRBENDELEMENT_H_
