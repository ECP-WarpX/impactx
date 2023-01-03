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

#include <AMReX_Random.H>
#include <AMReX_REAL.H>


namespace impactx
{
namespace distribution
{
    struct None
    {
        /** This distribution does nothing
         */
         None()
         {
         }

        /** Return 1 6D particle coordinate
         *
         * Does nothing to the parameters.
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
                [[maybe_unused]] amrex::ParticleReal & x,
                [[maybe_unused]] amrex::ParticleReal & y,
                [[maybe_unused]] amrex::ParticleReal & t,
                [[maybe_unused]] amrex::ParticleReal & px,
                [[maybe_unused]] amrex::ParticleReal & py,
                [[maybe_unused]] amrex::ParticleReal & pt,
                [[maybe_unused]] amrex::RandomEngine const& engine) const
        {
            /* nothing to do */
        }
    };

} // namespace distribution
} // namespace impactx

#endif // IMPACTX_DISTRIBUTION_NONE
