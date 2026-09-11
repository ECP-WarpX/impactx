/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_THINDIPOLE_H
#define IMPACTX_THINDIPOLE_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/gpuelement.H"
#include "mixin/alignment.H"
#include "mixin/beamoptic.H"
#include "mixin/thin.H"
#include "mixin/lineartransport.H"
#include "mixin/spintransport.H"
#include "mixin/nofinalize.H"
#include "mixin/pipeaperture.H"

#include <AMReX_Extension.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>
#include <AMReX_SIMD.H>

#include <cmath>
#include <stdexcept>

namespace impactx::elements
{
    struct ThinDipoleGPU
    : public mixin::BeamOptic<ThinDipoleGPU>,
      public mixin::LinearTransport<ThinDipoleGPU>,
      public mixin::Thin,
      public mixin::Alignment,
      public mixin::PipeAperture,
      public mixin::SpinTransport,
      public mixin::NoFinalize,
      public amrex::simd::Vectorized<amrex::simd::native_simd_size_particlereal>
    {
        static constexpr auto type = "ThinDipole";
        using PType = ImpactXParticleContainer::ParticleType;

        /** A general thin-kick dipole element with chromatic effects
         *
         * This model is equivalent to the model described in Section 3.1 of
         * G. Ripken, F. Schmidt, "A Symplectic Six-Dimensional Thin-Lens Formalism
         * for Tracking," CERN/SL/95-12 (AP), 1995.  It is intended to replicate
         * the thin-lens dipole model in MAD-X.
         *
         * @param theta - the total bending angle (degrees)
         * @param rc - the curvature radius (m)
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param aperture_x horizontal half-aperture in m
         * @param aperture_y vertical half-aperture in m
         */
        ThinDipoleGPU (
            amrex::ParticleReal theta,
            amrex::ParticleReal rc,
            amrex::ParticleReal dx = DEFAULT_dx,
            amrex::ParticleReal dy = DEFAULT_dy,
            amrex::ParticleReal rotation_degree = DEFAULT_rotation_degree,
            amrex::ParticleReal aperture_x = DEFAULT_aperture_x,
            amrex::ParticleReal aperture_y = DEFAULT_aperture_y
        )
        : Alignment(dx, dy, rotation_degree),
          PipeAperture(aperture_x, aperture_y),
          m_theta(theta * degree2rad),
          m_rc(rc)
        {
        }

        /** Reverse the element (negate bending angle) */
        void reverse () { m_theta = -m_theta; }

        /** Push all particles */
        using BeamOptic::operator();

        /** Compute and cache the constants for the push.
         *
         * In particular, used to pre-compute and cache variables that are
         * independent of the individually tracked particle.
         *
         * @param refpart reference particle
         */
        void compute_constants (RefPart const & refpart)
        {
            using namespace amrex::literals; // for _rt and _prt

            Alignment::compute_constants(refpart);

            // access reference particle to find relativistic beta
            m_beta = refpart.beta();

            // compute the effective (equivalent) arc length and curvature
            m_ds = m_theta * m_rc;
            m_kx = 1.0_prt / m_rc;
        }

        /** This is a multipole functor, so that a variable of this type can be used like a
         *  multipole function.
         *
         * The @see compute_constants method must be called before pushing particles through this operator.
         *
         * @param x particle position in x
         * @param y particle position in y (unused)
         * @param t particle position in t
         * @param px particle momentum in x
         * @param py particle momentum in y (unused)
         * @param pt particle momentum in t
         * @param idcpu particle global index (unused)
         * @param refpart reference particle (unused)
         */
        template<typename T_Real=amrex::ParticleReal, typename T_IdCpu=uint64_t>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
            T_Real const & AMREX_RESTRICT x,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT y,
            T_Real & AMREX_RESTRICT t,
            T_Real & AMREX_RESTRICT px,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT py,
            T_Real const & AMREX_RESTRICT pt,
            [[maybe_unused]] T_IdCpu const & AMREX_RESTRICT idcpu,
            [[maybe_unused]] RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using amrex::Math::powi;
            using namespace std; // for cmath(float)

            // access position data
            T_Real const xout = x;
            // T_Real const yout = y;
            T_Real const tout = t;

            // initialize output values of momenta
            T_Real pxout = px;
            // T_Real pyout = py;
            // T_Real ptout = pt;

            // compute the function expressing dp/p in terms of pt (labeled f in Ripken etc.)
            T_Real const f = -1.0_prt + sqrt(1.0_prt - 2.0_prt * pt / m_beta + powi<2>(pt));
            T_Real const fprime = (1.0_prt - m_beta*pt) / (m_beta * (1.0_prt + f));

            // advance position and momentum
            // x = xout;
            pxout = px - powi<2>(m_kx) * m_ds * xout + m_kx * m_ds * f; //eq (3.2b)

            // y = yout;
            // pyout = py;

            t = tout + m_kx * xout * m_ds * fprime; //eq (3.2e)
            // ptout = pt;

            // assign updated momenta
            px = pxout;
            // py = pyout;
            // pt = ptout;
        }

        /** This pushes the reference particle. */
        using Thin::operator();

        /** This pushes the reference particle.
         *
         * @param[in,out] refpart reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (RefPart & AMREX_RESTRICT refpart) const
        {
            // assign input reference particle values
            amrex::ParticleReal const px = refpart.px;
            amrex::ParticleReal const pz = refpart.pz;

            // calculate expensive terms once
            auto const [sin_theta, cos_theta] = amrex::Math::sincos(m_theta);

            // advance position and momentum (thin dipole)
            refpart.px = px*cos_theta - pz*sin_theta;
            refpart.pz = pz*cos_theta + px*sin_theta;
        }


        /** This operator pushes the particle spin.
         *
         * @param x particle position in x
         * @param y particle position in y
         * @param t particle position in t
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t
         * @param sx particle spin x-component
         * @param sy particle spin y-component
         * @param sz particle spin z-component
         * @param idcpu particle global index
         * @param refpart reference particle
         */
        template<typename T_Real=amrex::ParticleReal, typename T_IdCpu=uint64_t>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void spin_and_phasespace_push (
            T_Real const & AMREX_RESTRICT x,
            T_Real const & AMREX_RESTRICT y,
            T_Real & AMREX_RESTRICT t,
            T_Real & AMREX_RESTRICT px,
            T_Real const & AMREX_RESTRICT py,
            T_Real const & AMREX_RESTRICT pt,
            T_Real & AMREX_RESTRICT sx,
            T_Real & AMREX_RESTRICT sy,
            T_Real & AMREX_RESTRICT sz,
            T_IdCpu const & AMREX_RESTRICT idcpu,
            RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using namespace std; // for cmath(float)

            // reference particle relativistic factors
            amrex::ParticleReal const gamma_ref = refpart.gamma();
            amrex::ParticleReal const gyro_anomaly = refpart.gyromagnetic_anomaly;

            // initialize the three components of the axis-angle vector
            T_Real lambdax = 0_prt;
            T_Real lambday = 0_prt;
            T_Real lambdaz = 0_prt;

            // Magnetic field normalized by q/mc [1/m]:
            T_Real const Bx = 0_prt;
            T_Real const By = m_kx * m_beta * gamma_ref;
            T_Real const Bz = 0_prt;

            // Electric field normalized by q/mc^2 [1/m]:
            T_Real const Ex = 0.0_prt;
            T_Real const Ey = 0.0_prt;
            T_Real const Ez = 0.0_prt;

            // The Thomas-BMT generator is evaluated at the horizontal momentum
            // halfway through the kick, so that reversing the element (negating
            // the bend angle) produces the exact inverse rotation.
            T_Real const px_in = px;

            // phase space push
            (*this)(x, y, t, px, py, pt, idcpu, refpart);

            // horizontal momentum at the kick midpoint
            T_Real const px_mid = 0.5_prt * (px_in + px);

            // Quantities required to evaluate the full Thomas-BMT precession
            // vector, at the midpoint horizontal momentum.
            T_Real const Pnorm = sqrt(1_prt - 2_prt * pt / m_beta + pt*pt);
            T_Real const iPnorm = 1_prt/Pnorm;
            T_Real const gamma = gamma_ref * (1_prt - pt*m_beta);
            T_Real const Ps = sqrt(Pnorm*Pnorm - px_mid*px_mid - py*py);
            T_Real const ux = px_mid * iPnorm;
            T_Real const uy = py * iPnorm;
            T_Real const uz = Ps * iPnorm;

            // Design-orbit curvature [1/m]
            amrex::ParticleReal const h = m_kx;

            // Evaluation of the full Thomas-BMT precession vector.
            tbmt_precession_vector(x,ux,uy,uz,gamma,h,gyro_anomaly,Bx,By,Bz,Ex,Ey,Ez,lambdax,lambday,lambdaz);

            // generator over the effective (equivalent) arc length
            lambdax *= m_ds;
            lambday *= m_ds;
            lambdaz *= m_ds;

            // push the spin vector using the generator just determined
            rotate_spin(lambdax, lambday, lambdaz, sx, sy, sz);
        }


        /** This pushes the covariance matrix. */
        using LinearTransport::operator();

        /** This function returns the linear transport map.
         *
         * Linearization of the Ripken-Schmidt thin-dipole kick
         * about the reference trajectory. With the
         * effective arc length ds = theta * rc and curvature kx = 1/rc,
         * the thin kick is
         *   dpx = - kx^2 * ds * x + kx * ds * f(pt)
         *   dt  = + kx * ds * x * f'(pt)
         * where f(pt) = -1 + sqrt(1 - 2*pt/beta + pt^2) satisfies
         * f(0) = 0, f'(0) = 1/beta. The Jacobian at the origin yields:
         *   R(2,1) = - kx^2 * ds     = - theta / rc
         *   R(2,6) = - kx * ds / beta = - theta / beta
         *   R(5,1) = + kx * ds / beta = + theta / beta
         * All other entries are identity.
         *
         * Ideal behavior only: alignment and rotation errors and
         * feed-down from a nonzero reference orbit are ignored; the
         * linearization is taken at (x, y, t) = (0, 0, 0) and pt = 0.
         *
         * @param[in] refpart reference particle
         * @returns 6x6 transport matrix
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        Map6x6
        transport_map (RefPart const & AMREX_RESTRICT refpart) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using amrex::Math::powi;

            amrex::ParticleReal const beta = refpart.beta();
            amrex::ParticleReal const ds = m_theta * m_rc;     // effective arc length
            amrex::ParticleReal const kx = 1.0_prt / m_rc;     // curvature

            Map6x6 R = Map6x6::Identity();
            R(2,1) = -powi<2>(kx) * ds;      // geometric focusing
            R(2,6) = -kx * ds / beta;        // chromatic dispersion kick
            R(5,1) = +kx * ds / beta;        // path-length coupling to x

            // apply the transverse rotation (roll) alignment error
            return rotate_aligned_map(R);
        }

        amrex::ParticleReal m_theta;  //! dipole bending angle (rad)
        amrex::ParticleReal m_rc;     //! curvature radius (m)

      private:
        // constants that are independent of the individually tracked particle,
        // see: compute_constants() to refresh
        amrex::ParticleReal m_beta;  //! beta of the reference particle
        amrex::ParticleReal m_ds;    //! effective (equivalent) arc length
        amrex::ParticleReal m_kx;    //! effective (equivalent) curvature
    };

    /** ThinDipole: GPU state plus common host-only metadata */
    using ThinDipole = Element<ThinDipoleGPU>;

} // namespace impactx

IMPACTX_PUSH_EXTERN_TEMPLATE(impactx::elements::ThinDipoleGPU)

#endif // IMPACTX_THINDIPOLE_H
