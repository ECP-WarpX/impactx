/* Copyright 2021 Ji Qiang, Axel Huebl
 *
 * This file is part of ImpactX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_DISTRIBUTION_WATERBAG
#define IMPACTX_DISTRIBUTION_WATERBAG

#include <AMReX_Random.H>
#include <AMReX_REAL.H>

#include <cmath>


namespace impactx
{
namespace distribution
{
    struct Waterbag
    {
        /** A Waterbag distribution
         *
         * return 1 particle coordinates sampling from a 6D Waterbag distribution
         * with mean = 0, std = 1.
         */
        Waterbag() = default;

        /** Return 1 6D particle coordinate
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
                amrex::RandomEngine const& engine) const {

            using namespace amrex::literals;

            amrex::ParticleReal r1=0.0, r2=0.0, r3=0.0;
            amrex::ParticleReal x1=0.0, x2=0.0, x3=0.0;

            for(;;)
            {
                r1 = 2_prt * amrex::Random(engine) - 1.0_prt;
                r2 = 2_prt * amrex::Random(engine) - 1.0_prt;
                r3 = 2_prt * amrex::Random(engine) - 1.0_prt;
                x1 = 2_prt * amrex::Random(engine) - 1.0_prt;
                x2 = 2_prt * amrex::Random(engine) - 1.0_prt;
                x3 = 2_prt * amrex::Random(engine) - 1.0_prt;
                if(r1*r1 + r2*r2 + r3*r3 + x1*x1 + x2*x2 + x3*x3 < 1.0_prt) break;
            }
            x = sqrt(8.0_prt) * r1;
            y = sqrt(8.0_prt) * r2;
            t = sqrt(8.0_prt) * r3;
            px = sqrt(8.0_prt) * x1;
            py = sqrt(8.0_prt) * x2;
            pt = sqrt(8.0_prt) * x3;
        }
    };

} // namespace distribution
} // namespace impactx

#endif // IMPACTX_DISTRIBUTION_WATERBAG
