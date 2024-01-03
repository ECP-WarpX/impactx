/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_DISTRIBUTION_NONE
#define IMPACTX_DISTRIBUTION_NONE

#include "particles/ReferenceParticle.H"

#include <AMReX_Random.H>
#include <AMReX_REAL.H>


namespace impactx::distribution
{
    struct None
    {
        /** This distribution sets all values to zero.
         */
        None()
        {
        }

        /** Initialize the distribution.
         *
         * Nothing to do here.
         *
         * @param bunch_charge charge of the beam in C
         * @param ref the reference particle
         */
        void initialize ([[maybe_unused]] amrex::ParticleReal bunch_charge, [[maybe_unused]] RefPart const & ref)
        {
        }

        /** Return 1 6D particle coordinate
         *
         * Sets all values to zero.
         *
         * @param x particle position in x
         * @param y particle position in y
         * @param t particle position in t
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t
         * @param engine a random number engine (with associated state)
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
                amrex::ParticleReal & x,
                amrex::ParticleReal & y,
                amrex::ParticleReal & t,
                amrex::ParticleReal & px,
                amrex::ParticleReal & py,
                amrex::ParticleReal & pt,
                [[maybe_unused]] amrex::RandomEngine const& engine) const
        {
            using namespace amrex::literals;

            x = 0_prt;
            y = 0_prt;
            t = 0_prt;
            px = 0_prt;
            py = 0_prt;
            pt = 0_prt;
        }
    };

} // namespace impactx::distribution

#endif // IMPACTX_DISTRIBUTION_NONE
