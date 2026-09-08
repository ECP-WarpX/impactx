/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Alex Bojanich, Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_PARTICLES_WAKEFIELDS_SHORTRANGEWAKEELEMENT_H_
#define IMPACTX_PARTICLES_WAKEFIELDS_SHORTRANGEWAKEELEMENT_H_

#include "elements/RFCavity.H"
#include "elements/ShortRF.H"
#include "elements/Buncher.H"
#include "elements/ChrUniformAcc.H"
#include "particles/ReferenceParticle.H"

#include <AMReX_REAL.H>

#include <cmath>
#include <variant>


namespace impactx::particles::wakefields
{
    /** Non-owning view on the user-specified single-particle wake function of an element.
     *
     * Trivially copyable, so it can be captured by value into GPU kernels. The
     * pointers address memory owned by the element and are only valid for the
     * duration of the current call: do not cache them, re-query per use instead.
     */
    struct ShortRangeWakeView
    {
        bool has_wake = false;  //!< element contributes a short range structure wake
        int n_data = 0;  //!< number of tabulated values per wake function
        amrex::ParticleReal const * z = nullptr;  //!< longitudinal positions in m
        amrex::ParticleReal const * wake_x = nullptr;  //!< horizontal single-particle wake function
        amrex::ParticleReal const * wake_y = nullptr;  //!< vertical single-particle wake function
        amrex::ParticleReal const * wake_z = nullptr;  //!< longitudinal single-particle wake function
    };

    /** Function to check if an element supports short range structure wakefields and to calculate parameters needed for wake.
     *
     * @param[in] element_variant the lattice element to check
     * @param[in] refpart reference particle
     * @returns a view on the element's wake function data, with has_wake set if the element supplies one
     */
    template <typename VariantType>
    ShortRangeWakeView
    ShortRangeWakeElement (
        VariantType const& element_variant,
        [[maybe_unused]] RefPart const & refpart
    )
    {
        ShortRangeWakeView wake;

        if (std::holds_alternative<elements::RFCavity>(element_variant))
        {
            // Access element parameters as needed.
            auto const & element = std::get<elements::RFCavity>(element_variant);
            auto const & wake_data = *elements::RFCavity::WakeData::get(element.m_wake_id);

            // Upload to device (no-op on CPU builds) and take non-owning views.
            auto const & d_z = wake_data.z_data.device_const();

            // the wake data is optional: an element without a tabulated wake
            // function does not contribute a short-range wake
            if (!d_z.empty())
            {
                wake.has_wake = true;
                wake.n_data = static_cast<int>(d_z.size());
                wake.z = d_z.data();
                wake.wake_x = wake_data.wake_x_data.device_const().data();
                wake.wake_y = wake_data.wake_y_data.device_const().data();
                wake.wake_z = wake_data.wake_z_data.device_const().data();
            }
        }
        else if (std::holds_alternative<elements::ShortRF>(element_variant))
        {
            // Access element parameters as needed.
            // auto& element = std::get<elements::ShortRF>(element_variant);
            // TODO:  not supported yet
        }
        else if (std::holds_alternative<elements::Buncher>(element_variant))
        {
            // Access element parameters as needed.
            // auto& element = std::get<elements::Buncher>(element_variant);
            // TODO:  not supported yet
        }
        else if (std::holds_alternative<elements::ChrAcc>(element_variant))
        {
            // Access element parameters as needed.
            // auto& element = std::get<elements::ChrAcc>(element_variant);
            // TODO:  not supported yet
        }

        return wake;
    }

} // namespace impactx::particles::wakefields

#endif // IMPACTX_PARTICLES_WAKEFIELDS_SHORTRANGEWAKEELEMENT_H_
