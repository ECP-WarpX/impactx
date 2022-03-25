/* Copyright 2022 Chad Mitchell, Axel Huebl
 *
 * This file is part of ImpactX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_DISTRIBUTION_KURTH6D
#define IMPACTX_DISTRIBUTION_KURTH6D

#include <AMReX_Random.H>
#include <AMReX_REAL.H>

#include <cmath>


namespace impactx
{
namespace distribution
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
         * @param sigx,sigy,sigt for zero correlation, these are the related
         *                                   RMS sizes (in meters)
         * @param sigpx,sigpy,sigpt RMS momentum
         * @param muxpx,muypy,mutpt correlation length-momentum
         */
         Kurth6D(amrex::ParticleReal const sigx, amrex::ParticleReal const sigy,
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

            amrex::ParticleReal v,costheta,sintheta,phi,r;
            amrex::ParticleReal L,alpha,pmax,pr,beta,p1,p2;
            amrex::ParticleReal root,a1,a2;

            constexpr amrex::ParticleReal pi = 3.14159265358979_prt;

            // Random samples used to define (x,y,z):
            v = amrex::Random(engine);
            costheta = amrex::Random(engine);
            costheta = 2_prt*(costheta-0.5_prt);
            sintheta = sqrt(1_prt-pow(costheta,2));
            phi = amrex::Random(engine);
            phi = 2_prt*pi*phi;

            // Transformations for (x,y,t):
            r = pow(v,1_prt/3_prt);
            x = r*sintheta*cos(phi);
            y = r*sintheta*sin(phi);
            t = r*costheta;

            // Random samples used to define L:
            L = amrex::Random(engine);
            L = r*sqrt(L);

            // Random samples used to define pr:
            alpha = amrex::Random(engine);
            alpha = pi*alpha;
            pmax = 1_prt - pow(L/r,2) - pow(r,2) + pow(L,2);
            pmax = sqrt(pmax);
            pr = pmax*cos(alpha);

            // Random samples used to define ptangent:
            beta = amrex::Random(engine);
            beta = 2_prt*pi*beta;
            p1 = L/r*cos(beta);  // This is phi component
            p2 = L/r*sin(beta);  // This is theta component

            // Transformation from spherical to Cartesian coord.:
            px = pr*sintheta*cos(phi) + p2*costheta*cos(phi) - p1*sin(phi);
            py = pr*sintheta*sin(phi) + p2*costheta*sin(phi) + p1*cos(phi);
            pt = pr*costheta - p2*sintheta;

            // Scale to produce the identity covariance matrix:
            amrex::ParticleReal c = sqrt(5.0_prt);
            x = c*x;
            y = c*y;
            t = c*t;
            px = c*px;
            py = c*py;
            pt = c*pt;

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

#endif // IMPACTX_DISTRIBUTION_KURTH6D
