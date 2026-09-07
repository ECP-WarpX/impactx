/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_DISTRIBUTION_GAUSSIAN
#define IMPACTX_DISTRIBUTION_GAUSSIAN

#include "particles/ReferenceParticle.H"

#include <ablastr/constant.H>

#include <AMReX_Random.H>
#include <AMReX_Algorithm.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <limits>


namespace impactx::distribution
{
    struct Gaussian
    {
        /** A Gaussian distribution
         *
         * Return sampling from a 6D Gaussian distribution, optionally with
         * truncation in each phase plane.
         *
         * @param lambdax,lambday,lambdat phase space axis intercepts (position);
         *                                for zero correlation, these are the related RMS sizes (in meters)
         * @param lambdapx,lambdapy,lambdapt phase space axis intercepts (momentum)
         * @param muxpx,muypy,mutpt correlation length-momentum
         * @param meanx,meany,meant offsets of the mean (centroid) positions from those of the reference particle
         * @param meanpx,meanpy,meanpt offsets of the mean (centroid) momenta from those of the reference particle
         * @param dispx,disppx,dispy,disppy dispersion and its derivative in horizontal and vertical directions
         * @param cutx,cuty,cutt dimensionless cut parameter (number of sigma) in each phase plane;
         *                       a value of 0 means no truncation (default)
         */
        Gaussian (
            amrex::ParticleReal lambdax,
            amrex::ParticleReal lambday,
            amrex::ParticleReal lambdat,
            amrex::ParticleReal lambdapx,
            amrex::ParticleReal lambdapy,
            amrex::ParticleReal lambdapt,
            amrex::ParticleReal muxpx=0.0,
            amrex::ParticleReal muypy=0.0,
            amrex::ParticleReal mutpt=0.0,
            amrex::ParticleReal meanx=0.0,
            amrex::ParticleReal meany=0.0,
            amrex::ParticleReal meant=0.0,
            amrex::ParticleReal meanpx=0.0,
            amrex::ParticleReal meanpy=0.0,
            amrex::ParticleReal meanpt=0.0,
            amrex::ParticleReal dispx=0.0,
            amrex::ParticleReal disppx=0.0,
            amrex::ParticleReal dispy=0.0,
            amrex::ParticleReal disppy=0.0,
            amrex::ParticleReal cutx=0.0,
            amrex::ParticleReal cuty=0.0,
            amrex::ParticleReal cutt=0.0
        )
        : m_lambdaX(lambdax), m_lambdaY(lambday), m_lambdaT(lambdat), m_lambdaPx(lambdapx), m_lambdaPy(lambdapy),
          m_lambdaPt(lambdapt), m_muxpx(muxpx), m_muypy(muypy), m_mutpt(mutpt), m_meanx(meanx), m_meany(meany),
          m_meant(meant), m_meanpx(meanpx), m_meanpy(meanpy), m_meanpt(meanpt), m_dispx(dispx), m_disppx(disppx),
          m_dispy(dispy), m_disppy(disppy), m_cutX(cutx), m_cutY(cuty), m_cutT(cutt)
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
            amrex::ParticleReal & AMREX_RESTRICT x,
            amrex::ParticleReal & AMREX_RESTRICT y,
            amrex::ParticleReal & AMREX_RESTRICT t,
            amrex::ParticleReal & AMREX_RESTRICT px,
            amrex::ParticleReal & AMREX_RESTRICT py,
            amrex::ParticleReal & AMREX_RESTRICT pt,
            amrex::RandomEngine const & engine
        ) const
        {
            using namespace amrex::literals;
            using ablastr::constant::math::pi;

            amrex::ParticleReal ln1,u1,u2;
            amrex::ParticleReal root,a1,a2;

            // Generate three pairs of (optionally truncated) normal random variables using Box-Muller:

            u1 = amrex::Random(engine);
            u1 = amrex::max(u1, std::numeric_limits<amrex::ParticleReal>::min());  // guard log(0): amrex::Random() may be exactly 0 in single precision
            u2 = amrex::Random(engine);
            if (m_cutX > 0.0_prt) {
                amrex::ParticleReal f = m_cutX*m_cutX/2.0_prt;
                amrex::ParticleReal ef = std::exp(-f);
                amrex::ParticleReal r2scale = 1_prt - f*ef/(1_prt - ef);
                amrex::ParticleReal r2 = -2_prt * std::log(1_prt - u1*(1_prt - ef));
                ln1 = std::sqrt(r2 / r2scale);
            } else {
                ln1 = std::sqrt(-2_prt*std::log(u1));
            }
            x = ln1 * std::cos(2_prt*pi*u2);
            px = ln1 * std::sin(2_prt*pi*u2);

            u1 = amrex::Random(engine);
            u1 = amrex::max(u1, std::numeric_limits<amrex::ParticleReal>::min());  // guard log(0): amrex::Random() may be exactly 0 in single precision
            u2 = amrex::Random(engine);
            if (m_cutY > 0.0_prt) {
                amrex::ParticleReal f = m_cutY*m_cutY/2.0_prt;
                amrex::ParticleReal ef = std::exp(-f);
                amrex::ParticleReal r2scale = 1_prt - f*ef/(1_prt - ef);
                amrex::ParticleReal r2 = -2_prt * std::log(1_prt - u1*(1_prt - ef));
                ln1 = std::sqrt(r2 / r2scale);
            } else {
                ln1 = std::sqrt(-2_prt*std::log(u1));
            }
            y = ln1 * std::cos(2_prt*pi*u2);
            py = ln1 * std::sin(2_prt*pi*u2);

            u1 = amrex::Random(engine);
            u1 = amrex::max(u1, std::numeric_limits<amrex::ParticleReal>::min());  // guard log(0): amrex::Random() may be exactly 0 in single precision
            u2 = amrex::Random(engine);
            if (m_cutT > 0.0_prt) {
                amrex::ParticleReal f = m_cutT*m_cutT/2.0_prt;
                amrex::ParticleReal ef = std::exp(-f);
                amrex::ParticleReal r2scale = 1_prt - f*ef/(1_prt - ef);
                amrex::ParticleReal r2 = -2_prt * std::log(1_prt - u1*(1_prt - ef));
                ln1 = std::sqrt(r2 / r2scale);
            } else {
                ln1 = std::sqrt(-2_prt*std::log(u1));
            }
            t = ln1 * std::cos(2_prt*pi*u2);
            pt = ln1 * std::sin(2_prt*pi*u2);

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
        amrex::ParticleReal m_muxpx, m_muypy, m_mutpt;  //! correlation length-momentum
        amrex::ParticleReal m_meanx, m_meany, m_meant;  //! spatial coordinates of centroid offset
        amrex::ParticleReal m_meanpx, m_meanpy, m_meanpt;  //! momentum coordinates of centroid offset
        amrex::ParticleReal m_dispx, m_disppx, m_dispy, m_disppy;  //! dispersion and its derivative
        amrex::ParticleReal m_cutX, m_cutY, m_cutT;  //! dimensionless cutoff for each phase plane; 0 means no cut
    };

} // namespace impactx::distribution

#endif // IMPACTX_DISTRIBUTION_GAUSSIAN
