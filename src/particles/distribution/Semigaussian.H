/* Copyright 2022 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_DISTRIBUTION_SEMIGAUSSIAN
#define IMPACTX_DISTRIBUTION_SEMIGAUSSIAN

#include <AMReX_Random.H>
#include <AMReX_REAL.H>

#include <cmath>


namespace impactx
{
namespace distribution
{
    struct Semigaussian
    {
        /** A Semi-Gaussian distribution
         *
         * Return sampling from a 6D Semi-Gaussian distribution (uniform in
         * position, Gaussian in momentum).
         *
         * @param sigx,sigy,sigt for zero correlation, these are the related
         *                                   RMS sizes (in meters)
         * @param sigpx,sigpy,sigpt RMS momentum
         * @param muxpx,muypy,mutpt correlation length-momentum
         */
         Semigaussian(amrex::ParticleReal const sigx, amrex::ParticleReal const sigy,
         amrex::ParticleReal const sigt,amrex::ParticleReal const sigpx,
         amrex::ParticleReal const sigpy,amrex::ParticleReal const sigpt,
         amrex::ParticleReal const muxpx=0.0, amrex::ParticleReal const muypy=0.0,
         amrex::ParticleReal const mutpt=0.0
         )
          : m_sigmaX(sigx),m_sigmaY(sigy),m_sigmaT(sigt),m_sigmaPx(sigpx),m_sigmaPy(sigpy),
            m_sigmaPt(sigpt),m_muxpx(muxpx),m_muypy(muypy),m_mutpt(mutpt)
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
                amrex::RandomEngine const& engine) const {

            using namespace amrex::literals;

            amrex::ParticleReal ln1,u1,u2;
            amrex::ParticleReal root,a1,a2;
            amrex::ParticleReal phi,v,r;

            constexpr amrex::ParticleReal pi = 3.14159265358979_prt;

            // Generate a 3D uniform distribution (x,y,t) within a cylinder:

            phi = amrex::Random(engine);
            phi = 2_prt*pi*phi;
            v = amrex::Random(engine);
            r = sqrt(v);
            x = r*cos(phi);
            y = r*sin(phi);
            t = amrex::Random(engine);
            t = 2_prt*(t-0.5_prt);

            // Scale to produce the identity covariance matrix:
            amrex::ParticleReal c = sqrt(3.0_prt);
            x = 2_prt*x;
            y = 2_prt*y;
            t = c*t;

            // Generate three normal random variables (px,py,pt) using Box-Muller:

            u1 = amrex::Random(engine);
            u2 = amrex::Random(engine);
            ln1 = sqrt(-2_prt*log(u1));
            px = ln1*cos(2_prt*pi*u2);
            py = ln1*sin(2_prt*pi*u2);

            u1 = amrex::Random(engine);
            u2 = amrex::Random(engine);
            ln1 = sqrt(-2_prt*log(u1));
            pt = ln1*cos(2_prt*pi*u2);

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
        amrex::ParticleReal m_sigmaX,m_sigmaY,m_sigmaT; //! related RMS sizes (length)
        amrex::ParticleReal m_sigmaPx,m_sigmaPy,m_sigmaPt; //! RMS momentum
        amrex::ParticleReal m_muxpx,m_muypy,m_mutpt; //! correlation length-momentum
    };

} // namespace distribution
} // namespace impactx

#endif // IMPACTX_DISTRIBUTION_SEMIGAUSSIAN
