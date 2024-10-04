/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENT_MARKER_H
#define IMPACTX_ELEMENT_MARKER_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/thin.H"
#include "mixin/named.H"
#include "mixin/nofinalize.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>


namespace impactx
{
    struct Marker
    : public elements::Named,
      public elements::Thin,
      public elements::NoFinalize
    {
        static constexpr auto type = "Marker";
        using PType = ImpactXParticleContainer::ParticleType;

        /** This named element does nothing.
         *
         * @param name a user defined and not necessarily unique name of the element
         */
        Marker (std::string name)
        : Named(name)
        {
        }

        /** Push all particles - nothing to do here */
        void operator() (
            ImpactXParticleContainer & /* pc */,
            int /* step */,
            int /* cycle */
        ) {
            // nothing to do
        }

        /** Push all particles - nothing to do here */
        void operator() (
            ImpactXParticleContainer::iterator & /* pti */,
            RefPart & AMREX_RESTRICT /* ref_part */
        ) {
            // nothing to do
        }

        /** Does nothing to a particle.
         *
         * @param x particle position in x
         * @param y particle position in y
         * @param t particle position in t
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t
         * @param idcpu particle global index (unused)
         * @param refpart reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
            [[maybe_unused]] amrex::ParticleReal & AMREX_RESTRICT x,
            [[maybe_unused]] amrex::ParticleReal & AMREX_RESTRICT y,
            [[maybe_unused]] amrex::ParticleReal & AMREX_RESTRICT t,
            [[maybe_unused]] amrex::ParticleReal & AMREX_RESTRICT px,
            [[maybe_unused]] amrex::ParticleReal & AMREX_RESTRICT py,
            [[maybe_unused]] amrex::ParticleReal & AMREX_RESTRICT pt,
            [[maybe_unused]] uint64_t & AMREX_RESTRICT idcpu,
            [[maybe_unused]] RefPart const & refpart
        ) const
        {
            // nothing to do
        }

        /** This pushes the reference particle. */
        using Thin::operator();
    };

} // namespace impactx

#endif // IMPACTX_ELEMENT_MARKER_H
