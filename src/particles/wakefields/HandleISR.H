/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Alex Bojanich, Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef HANDLE_ISR_H
#define HANDLE_ISR_H

#include "particles/ImpactXParticleContainer.H"
#include "BendRadius.H"
#include "ISRPush.H"

#include <cmath>
#include <vector>


namespace impactx::particles::wakefields
{
    /** Function to handle ISR in bending dipoles.
     *
     * Returns immediately if slice_ds is zero.
     *
     * @param[in] particle_container the particle species container
     * @param[in] element_variant variant type of the lattice element
     * @param[in] slice_ds slice spacing along s
     */
    template <typename T_Element>
    void HandleISR (
        impactx::ImpactXParticleContainer& particle_container,
        T_Element const& element_variant,
        amrex::Real slice_ds
    )
    {
        if (slice_ds == 0.0) { return; }

        BL_PROFILE("impactx::particles::wakefields::HandleISR")

        amrex::ParmParse pp_algo("algo");
        bool isr = false;
        pp_algo.queryAdd("isr", isr);

        // ISR is applied only in elements that bend the reference orbit
        auto const rc = impactx::particles::wakefields::bend_radius(element_variant, particle_container.GetRefParticle());

        if (isr && rc.has_value())
        {
            amrex::Real const R = *rc;

            int isr_order = 1;
            pp_algo.queryAddWithParser("isr_order", isr_order);

            bool isr_on_ref_part = false;
            pp_algo.queryAdd("isr_on_ref_part", isr_on_ref_part);

            // Call function to kick particles with wake
            impactx::particles::wakefields::ISRPush(particle_container, slice_ds, R, isr_order, isr_on_ref_part);
        }
    }

} // namespace impactx::particles::wakefields

#endif // HANDLE_ISR_H
