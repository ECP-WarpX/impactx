/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_DISTRIBUTION_KURTH4D
#define IMPACTX_DISTRIBUTION_KURTH4D

#include "particles/ReferenceParticle.H"

#include <ablastr/constant.H>

#include <AMReX_Random.H>
#include <AMReX_REAL.H>

#include <cmath>


namespace impactx::distribution
{
    struct Kurth4D
    {
        /** A 4D Kurth distribution transversely + a uniform distribution
         *  it t + a Gaussian distribution in pt
         *
         * Return sampling from a 4D Kurth + uniform distribution.
         *
         * @param sigx,sigy,sigt for zero correlation, these are the related
         *                                   RMS sizes (in meters)
         * @param sigpx,sigpy,sigpt RMS momentum
         * @param muxpx,muypy,mutpt correlation length-momentum
         */
        Kurth4D (
            amrex::ParticleReal sigx,
            amrex::ParticleReal sigy,
            amrex::ParticleReal sigt,
            amrex::ParticleReal sigpx,
            amrex::ParticleReal sigpy,
            amrex::ParticleReal sigpt,
            amrex::ParticleReal muxpx=0.0,
            amrex::ParticleReal muypy=0.0,
            amrex::ParticleReal mutpt=0.0
        )
        : m_sigmaX(sigx), m_sigmaY(sigy), m_sigmaT(sigt), m_sigmaPx(sigpx), m_sigmaPy(sigpy),
          m_sigmaPt(sigpt), m_muxpx(muxpx), m_muypy(muypy), m_mutpt(mutpt)
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

            amrex::ParticleReal v,phi,r,u1,u2,ln1;
            amrex::ParticleReal alpha,u,Lz,pmax,pr,pphi;
            amrex::ParticleReal root,a1,a2;

            // Sample and transform to define (x,y):
            v = amrex::Random(engine);
            phi = amrex::Random(engine);
            phi = 2_prt*pi*phi;
            r = sqrt(v);
            x = r*cos(phi);
            y = r*sin(phi);

            // Random samples used to define Lz:
            u = amrex::Random(engine);
            Lz = r*(2.0_prt*u-1.0_prt);

            // Random samples used to define pr:
            alpha = amrex::Random(engine);
            alpha = pi*alpha;
            pmax = 1.0_prt - pow((Lz/r),2) - pow(r,2) + pow(Lz,2);
            pmax = sqrt(pmax);
            pr = pmax*cos(alpha);
            pphi = Lz/r;

            // Transformations used to obtain (px,py):
            px = pr*cos(phi)-pphi*sin(phi);
            py = pr*sin(phi)+pphi*cos(phi);

            // Sample and transform to define (t,pt):
            t = amrex::Random(engine);
            t = 2.0_prt*(t-0.5_prt);
            u1 = amrex::Random(engine);
            u2 = amrex::Random(engine);
            ln1 = sqrt(-2_prt*log(u1));
            pt = ln1*cos(2_prt*pi*u2);

            // Scale to produce the identity covariance matrix:
            amrex::ParticleReal const c = sqrt(3.0_prt);
            x = 2_prt*x;
            y = 2_prt*y;
            t = c*t;
            px = 2_prt*px;
            py = 2_prt*py;
            // pt = pt;

            // Transform to produce the desired second moments/correlations:
            root = sqrt(1.0_prt-m_muxpx*m_muxpx);
            a1 = m_sigmaX*x/root;
            a2 = m_sigmaPx*(-m_muxpx*x/root+px);
            x = a1;
            px = a2;
            root = sqrt(1.0_prt-m_muypy*m_muypy);
            a1 = m_sigmaY*y/root;
            a2 = m_sigmaPy*(-m_muypy*y/root+py);
            y = a1;
            py = a2;
            root = sqrt(1.0_prt-m_mutpt*m_mutpt);
            a1 = m_sigmaT*t/root;
            a2 = m_sigmaPt*(-m_mutpt*t/root+pt);
            t = a1;
            pt = a2;
        }

    private:
        amrex::ParticleReal m_sigmaX, m_sigmaY, m_sigmaT;  //! related RMS sizes (length)
        amrex::ParticleReal m_sigmaPx, m_sigmaPy, m_sigmaPt;  //! RMS momentum
        amrex::ParticleReal m_muxpx, m_muypy, m_mutpt;  //! correlation length-momentum
    };

} // namespace impactx::distribution

#endif // IMPACTX_DISTRIBUTION_KURTH4D
