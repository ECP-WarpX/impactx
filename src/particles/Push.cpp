/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "Push.H"

#include <AMReX_BLProfiler.H>
#include <AMReX_Extension.H>  // for AMREX_RESTRICT
#include <AMReX_REAL.H>       // for ParticleReal

#include <utility>


namespace impactx
{
    void Push (ImpactXParticleContainer & pc,
               KnownElements const & element_variant)
    {
        // performance profiling per element
        std::string element_name;
        std::visit([&element_name](auto&& element){ element_name = element.name; }, element_variant);
        std::string const profile_name = "impactx::Push::" + element_name;
        BL_PROFILE("impactx::Push");
        BL_PROFILE(profile_name);

        using namespace amrex::literals; // for _rt and _prt

        // preparing to access reference particle data: RefPart
        RefPart & ref_part = pc.GetRefParticle();

        // loop over refinement levels
        int const nLevel = pc.finestLevel();
        for (int lev = 0; lev <= nLevel; ++lev)
        {
            // get simulation geometry information
            //const amrex::Geometry& gm = this->Geom(lev);
            //const auto prob_lo = gm.ProbLo();

            // loop over all particle boxes
            using ParIt = ImpactXParticleContainer::iterator;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (ParIt pti(pc, lev); pti.isValid(); ++pti) {
                // here we just access the element by its respective type
                std::visit(
                    [&pti, &ref_part](auto element) {
                        // push reference particle in global coordinates
                        element(ref_part);

                        // push beam particles relative to reference particle
                        element(pti, ref_part);
                    },
                    element_variant
                );
            } // end loop over all particle boxes
        } // env mesh-refinement level loop
    }

} // namespace impactx
