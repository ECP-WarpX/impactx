/* Copyright 2022 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_NONE_H
#define IMPACTX_NONE_H

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>


namespace impactx
{
    struct None
    {
        using PType = ImpactXParticleContainer::ParticleType;

        /** This element does nothing.
         */
        None ()
        {
        }

        /** Does nothing to a particle.
         *
         * @param p Particle AoS data for positions and cpu/id
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t
         * @param refpart reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
                [[maybe_unused]] PType& AMREX_RESTRICT p,
                [[maybe_unused]] amrex::ParticleReal & AMREX_RESTRICT px,
                [[maybe_unused]] amrex::ParticleReal & AMREX_RESTRICT py,
                [[maybe_unused]] amrex::ParticleReal & AMREX_RESTRICT pt,
                [[maybe_unused]] RefPart const refpart) const
        {
            // nothing to do
        }


        /** This pushes the reference particle.
         *
         * @param[in,out] refpart reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
                [[maybe_unused]] RefPart & AMREX_RESTRICT refpart) const {

            // nothing to do: this is a zero-length element

        }

    };

} // namespace impactx

#endif // IMPACTX_NONE_H