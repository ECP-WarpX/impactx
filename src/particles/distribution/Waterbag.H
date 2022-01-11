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

//	Waterbag( amrex::ParticleReal const sigx, amrex::ParticleReal const sigpx )
//        : m_sigx(sigx), m_sigpx(sigpx)
//        {
//        }

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
            amrex::ParticleReal pi,ln1,norm,u1,u2;
            amrex::ParticleReal g1,g2,g3,g4,g5,g6;

	    pi = 3.14159265358979_prt
            u1 = amrex::Random(engine); 
            u2 = amrex::Random(engine);
	    ln1 = sqrt(-2_prt*log(u1));
	    g1 = ln1*cos(2_prt*pi*u2);
	    g2 = ln1*sin(2_prt*pi*u2);
            u1 = amrex::Random(engine); 
            u2 = amrex::Random(engine);
	    ln1 = sqrt(-2_prt*log(u1));
	    g3 = ln1*cos(2_prt*pi*u2);
	    g4 = ln1*sin(2_prt*pi*u2);
            u1 = amrex::Random(engine); 
            u2 = amrex::Random(engine);
	    ln1 = sqrt(-2_prt*log(u1));
	    g5 = ln1*cos(2_prt*pi*u2);
	    g6 = ln1*sin(2_prt*pi*u2);
	    norm = sqrt(g1*g1+g2*g2+g3*g3+g4*g4+g5*g5+g6*g6);
	    g1 /= norm;
	    g2 /= norm;
	    g3 /= norm;
	    g4 /= norm;
	    g5 /= norm;
	    g6 /= norm;
            u1 = sqrt(8_prt)*amrex::Random(engine);
            u2 = pow(u1,1_prt/6_prt);
	    x = g1*u2;
	    y = g2*u2;
	    t = g3*u2;
	    px = g4*u2;
	    py = g5*u2;
	    pt = g6*u2;
        }
    };

} // namespace distribution
} // namespace impactx

#endif // IMPACTX_DISTRIBUTION_WATERBAG
