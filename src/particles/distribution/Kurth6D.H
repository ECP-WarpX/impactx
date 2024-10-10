/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_DISTRIBUTION_KURTH6D
#define IMPACTX_DISTRIBUTION_KURTH6D

#include "particles/ReferenceParticle.H"

#include <ablastr/constant.H>

#include <AMReX_Random.H>
#include <AMReX_REAL.H>

#include <cmath>


namespace impactx::distribution
{
    struct Kurth6D
    {
        /** A 6D Kurth distribution
         *
         * Return sampling from a 6D Kurth distribution.
         *
         * R. Kurth, Quarterly of Applied Mathematics vol. 32, pp. 325-329 (1978)
         * C. Mitchell, K. Hwang and R. D. Ryne, IPAC2021, WEPAB248 (2021)
         *
         * @param lambdax,lambday,lambdat for zero correlation, these are the related
         *                                   RMS sizes (in meters)
         * @param lambdapx,lambdapy,lambdapt RMS momentum
         * @param muxpx,muypy,mutpt correlation length-momentum
         */
        Kurth6D (
            amrex::ParticleReal lambdax,
            amrex::ParticleReal lambday,
            amrex::ParticleReal lambdat,
            amrex::ParticleReal lambdapx,
            amrex::ParticleReal lambdapy,
            amrex::ParticleReal lambdapt,
            amrex::ParticleReal muxpx=0.0,
            amrex::ParticleReal muypy=0.0,
            amrex::ParticleReal mutpt=0.0
        )
        : m_lambdaX(lambdax), m_lambdaY(lambday), m_lambdaT(lambdat), m_lambdaPx(lambdapx), m_lambdaPy(lambdapy),
          m_lambdaPt(lambdapt), m_muxpx(muxpx), m_muypy(muypy), m_mutpt(mutpt)
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

            amrex::ParticleReal v,costheta,sintheta,phi,r;
            amrex::ParticleReal L,alpha,pmax,pr,beta,p1,p2;
            amrex::ParticleReal root,a1,a2;

            // Random samples used to define (x,y,z):
            v = amrex::Random(engine);
            costheta = amrex::Random(engine);
            costheta = 2_prt*(costheta-0.5_prt);
            sintheta = std::sqrt(1_prt-pow(costheta,2));
            phi = amrex::Random(engine);
            phi = 2_prt*pi*phi;

            // Transformations for (x,y,t):
            r = std::pow(v,1_prt/3_prt);
            x = r*sintheta * std::cos(phi);
            y = r*sintheta * std::sin(phi);
            t = r*costheta;

            // Random samples used to define L:
            L = amrex::Random(engine);
            L = r*std::sqrt(L);

            // Random samples used to define pr:
            alpha = amrex::Random(engine);
            alpha = pi*alpha;
            pmax = 1_prt - std::pow(L/r,2) - std::pow(r,2) + std::pow(L,2);
            pmax = std::sqrt(pmax);
            pr = pmax * std::cos(alpha);

            // Random samples used to define ptangent:
            beta = amrex::Random(engine);
            beta = 2_prt*pi*beta;
            p1 = L/r * std::cos(beta);  // This is phi component
            p2 = L/r * std::sin(beta);  // This is theta component

            // Transformation from spherical to Cartesian coord.:
            px = pr*sintheta * std::cos(phi) + p2*costheta * std::cos(phi) - p1 * std::sin(phi);
            py = pr*sintheta * std::sin(phi) + p2*costheta * std::sin(phi) + p1 * std::cos(phi);
            pt = pr*costheta - p2*sintheta;

            // Scale to produce the identity covariance matrix:
            amrex::ParticleReal const c = std::sqrt(5.0_prt);
            x = c*x;
            y = c*y;
            t = c*t;
            px = c*px;
            py = c*py;
            pt = c*pt;

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
        }

        amrex::ParticleReal m_lambdaX, m_lambdaY, m_lambdaT;  //! related position axis intercepts (length) of the phase space ellipse
        amrex::ParticleReal m_lambdaPx, m_lambdaPy, m_lambdaPt;  //! related momentum axis intercepts of the phase space ellipse
        amrex::ParticleReal m_muxpx, m_muypy, m_mutpt;  //! correlation length-momentum
    };

} // namespace impactx::distribution

#endif // IMPACTX_DISTRIBUTION_KURTH6D
