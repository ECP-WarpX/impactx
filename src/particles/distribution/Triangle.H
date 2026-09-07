/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Ji Qiang, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_DISTRIBUTION_TRIANGLE
#define IMPACTX_DISTRIBUTION_TRIANGLE

#include "particles/ReferenceParticle.H"

#include <ablastr/constant.H>

#include <AMReX_Random.H>
#include <AMReX_Algorithm.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <limits>


namespace impactx::distribution
{
    struct Triangle
    {
        /** A Triangle distribution for LPA applications.
         *
         * Return sampling from a ramped, triangular current profile
         * with a Gaussian energy spread (possibly correlated).  The
         * transverse distribution is a 4D waterbag.
         *
         * @param lambdax,lambday,lambdat for zero correlation, these are the related
         *                       RMS sizes (in meters)
         * @param lambdapx,lambdapy,lambdapt RMS momentum
         * @param muxpx,muypy,mutpt correlation length-momentum
         * @param meanx,meany,meant offsets of the mean (centroid) positions from those of the reference particle
         * @param meanpx,meanpy,meanpt offsets of the mean (centroid) momenta from those of the reference particle
         * @param dispx,disppx,dispy,disppy dispersion and its derivative in horizontal and vertical directions
         */
        Triangle (
                amrex::ParticleReal const lambdax, amrex::ParticleReal const lambday,
                amrex::ParticleReal const lambdat, amrex::ParticleReal const lambdapx,
                amrex::ParticleReal const lambdapy, amrex::ParticleReal const lambdapt,
                amrex::ParticleReal const muxpx=0.0, amrex::ParticleReal const muypy=0.0,
                amrex::ParticleReal const mutpt=0.0, amrex::ParticleReal const meanx=0.0,
                amrex::ParticleReal const meany=0.0, amrex::ParticleReal const meant=0.0,
                amrex::ParticleReal const meanpx=0.0, amrex::ParticleReal const meanpy=0.0,
                amrex::ParticleReal const meanpt=0.0, amrex::ParticleReal const dispx=0.0,
                amrex::ParticleReal const disppx=0.0, amrex::ParticleReal const dispy=0.0,
                amrex::ParticleReal const disppy=0.0
        )
          : m_lambdaX(lambdax), m_lambdaY(lambday), m_lambdaT(lambdat), m_lambdaPx(lambdapx), m_lambdaPy(lambdapy),
            m_lambdaPt(lambdapt), m_muxpx(muxpx), m_muypy(muypy), m_mutpt(mutpt), m_meanx(meanx), m_meany(meany),
          m_meant(meant), m_meanpx(meanpx), m_meanpy(meanpy), m_meanpt(meanpt), m_dispx(dispx), m_disppx(disppx),
          m_dispy(dispy), m_disppy(disppy)
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

        /** Close and deallocate all data and handles.
         *
         * Nothing to do here.
         */
        void
        finalize ()
        {
        }

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
            amrex::RandomEngine const& engine
        ) const
        {

            using namespace amrex::literals;
            using ablastr::constant::math::pi;

            amrex::ParticleReal ln1, norm, u0, u1, u2;
            amrex::ParticleReal g1, g2, g3, g4, g5;
            amrex::ParticleReal d, root, a1, a2;

            // Sample the t coordinate for a ramped triangular profile (unit
            // variance):
            u0 = amrex::Random(engine);
            t = std::sqrt(2_prt)*(2_prt-3_prt*std::sqrt(u0));

            // Generate five standard normal random variables using Box-Muller:
            u1 = amrex::Random(engine);
            u1 = amrex::max(u1, std::numeric_limits<amrex::ParticleReal>::min());  // guard log(0): amrex::Random() may be exactly 0 in single precision
            u2 = amrex::Random(engine);
            ln1 = std::sqrt(-2_prt*std::log(u1));
            g1 = ln1 * std::cos(2_prt*pi*u2);
            g2 = ln1 * std::sin(2_prt*pi*u2);
            u1 = amrex::Random(engine);
            u1 = amrex::max(u1, std::numeric_limits<amrex::ParticleReal>::min());  // guard log(0): amrex::Random() may be exactly 0 in single precision
            u2 = amrex::Random(engine);
            ln1 = std::sqrt(-2_prt*std::log(u1));
            g3 = ln1 * std::cos(2_prt*pi*u2);
            g4 = ln1 * std::sin(2_prt*pi*u2);
            u1 = amrex::Random(engine);
            u1 = amrex::max(u1, std::numeric_limits<amrex::ParticleReal>::min());  // guard log(0): amrex::Random() may be exactly 0 in single precision
            u2 = amrex::Random(engine);
            ln1 = std::sqrt(-2_prt*std::log(u1));
            g5 = ln1 * std::cos(2_prt*pi*u2);

            // Use one of these normal random variables for pt:
            pt = g5;

            // Normalize the rest to produce uniform samples on the unit sphere:
            norm = std::sqrt(g1*g1+g2*g2+g3*g3+g4*g4);
            g1 /= norm;
            g2 /= norm;
            g3 /= norm;
            g4 /= norm;

            // Scale to produce uniform samples in a 4D ball (unit variance):
            d = 4_prt;  // unit ball dimension
            u1 = amrex::Random(engine);   // uniform sample
            u2 = std::sqrt(d+2_prt) * std::pow(u1,1_prt/d);
            x = g1*u2;
            y = g2*u2;
            px = g3*u2;
            py = g4*u2;

            // Transform to produce the desired second moments/correlations:
            root = std::sqrt(1.0_prt-m_muxpx*m_muxpx);
            a1 = m_lambdaX * x / root;
            a2 = m_lambdaPx * (-m_muxpx * x / root + px);
            x = a1;
            px = a2;
            root = std::sqrt(1.0_prt-m_muypy*m_muypy);
            a1 = m_lambdaY * y / root;
            a2 = m_lambdaPy * (-m_muypy * y / root + py);
            y = a1;
            py = a2;
            root = std::sqrt(1.0_prt-m_mutpt*m_mutpt);
            a1 = m_lambdaT * t / root;
            a2 = m_lambdaPt * (-m_mutpt * t / root + pt);
            t = a1;
            pt = a2;

            // Apply dispersive correlations
            x = x - m_dispx * pt;
            px = px - m_disppx * pt;
            y = y - m_dispy * pt;
            py = py - m_disppy * pt;

            // Apply centroid translation
            x = x + m_meanx;
            px = px + m_meanpx;
            y = y + m_meany;
            py = py + m_meanpy;
            t = t + m_meant;
            pt = pt + m_meanpt;
        }

        amrex::ParticleReal m_lambdaX, m_lambdaY, m_lambdaT;  //! related position axis intercepts (length) of the phase space ellipse
        amrex::ParticleReal m_lambdaPx, m_lambdaPy, m_lambdaPt;  //! related momentum axis intercepts of the phase space ellipse
        amrex::ParticleReal m_muxpx, m_muypy, m_mutpt; //! correlation length-momentum
        amrex::ParticleReal m_meanx, m_meany, m_meant;  //! spatial coordinates of centroid offset
        amrex::ParticleReal m_meanpx, m_meanpy, m_meanpt;  //! momentum coordinates of centroid offset
        amrex::ParticleReal m_dispx, m_disppx, m_dispy, m_disppy;  //! dispersion and its derivative
    };

} // namespace impactx::distribution

#endif // IMPACTX_DISTRIBUTION_TRIANGLE
