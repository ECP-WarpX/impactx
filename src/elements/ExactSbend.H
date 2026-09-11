/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_EXACTSBEND_H
#define IMPACTX_EXACTSBEND_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/gpuelement.H"
#include "mixin/alignment.H"
#include "mixin/pipeaperture.H"
#include "mixin/beamoptic.H"
#include "mixin/thick.H"
#include "mixin/lineartransport.H"
#include "mixin/spintransport.H"
#include "mixin/nofinalize.H"

#include <AMReX_Extension.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>
#include <AMReX_SIMD.H>

#include <cmath>
#include <stdexcept>

namespace impactx::elements
{
    struct ExactSbendGPU
    : public mixin::BeamOptic<ExactSbendGPU>,
      public mixin::LinearTransport<ExactSbendGPU>,
      public mixin::Thick,
      public mixin::Alignment,
      public mixin::PipeAperture,
      public mixin::SpinTransport,
      public mixin::NoFinalize,
      public amrex::simd::Vectorized<amrex::simd::native_simd_size_particlereal>
    {
        static constexpr auto type = "ExactSbend";
        using PType = ImpactXParticleContainer::ParticleType;

        static constexpr amrex::ParticleReal DEFAULT_B = 0;

        /** The body of an ideal sector bend, using the exact nonlinear transfer map.
         *  The map corresponds to the map described in:
         *  D. L. Bruhwiler et al, in Proc. of EPAC 98, pp. 1171-1173 (1998).
         *
         *  In the ultrarelativistic limit, this is equivalent to the map described in:
         *  E. Forest et al, Part. Accel. 45, pp. 65-94 (1994).
         *
         *  The model consists of a uniform bending field B_y with a hard edge.
         *  Pole faces are normal to the entry and exit velocity of the reference particle.
         *
         * @param ds Segment length in m.
         * @param phi Segment bend angle in degrees.
         * @param B Magnetic field strength in T.
         *  When B = 0 (default), the reference bending radius is defined by r0 = length / (angle in rad), corresponding to
         *  a magnetic field of B = rigidity / r0; otherwise the reference bending radius is defined by r0 = rigidity / B.
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param aperture_x horizontal half-aperture in m
         * @param aperture_y vertical half-aperture in m
         * @param nslice number of slices used for the application of space charge
         */
        ExactSbendGPU (
            amrex::ParticleReal ds,
            amrex::ParticleReal phi,
            amrex::ParticleReal B = DEFAULT_B,
            amrex::ParticleReal dx = DEFAULT_dx,
            amrex::ParticleReal dy = DEFAULT_dy,
            amrex::ParticleReal rotation_degree = DEFAULT_rotation_degree,
            amrex::ParticleReal aperture_x = DEFAULT_aperture_x,
            amrex::ParticleReal aperture_y = DEFAULT_aperture_y,
            int nslice = DEFAULT_nslice
        )
        : Thick(ds, nslice),
          Alignment(dx, dy, rotation_degree),
          PipeAperture(aperture_x, aperture_y),
          m_phi(phi * degree2rad), m_B(B)
        {
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::ParticleReal
        rc (RefPart const & refpart) const
        {
            using namespace amrex::literals; // for _rt and _prt

            return m_B != 0_prt ? refpart.rigidity_Tm() / m_B : m_ds / m_phi;
        }

        /** Reverse the element (negate ds and phi) */
        void reverse () { Thick::reverse(); m_phi = -m_phi; }

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
            using amrex::Math::powi;

            Alignment::compute_constants(refpart);

            // reference particle's orbital radius
            m_rc = this->rc(refpart);

            // arc length and angle of the current slice
            m_slice_ds = m_ds / nslice();
            amrex::ParticleReal const bend_sign = m_ds * m_rc / std::abs(m_ds * m_rc);
            m_slice_phi = bend_sign * std::abs(m_phi) / nslice();  // the sign is determined by rc/ds

            // trigonometric evaluations
            auto const [sin_phi, cos_phi] = amrex::Math::sincos(m_slice_phi);
            m_sin_phi = sin_phi;
            m_cos_phi = cos_phi;

            // access reference particle values to find relativistic factors
            m_ibeta = 1_prt / refpart.beta();
            m_ibetagam2 = 1_prt / powi<2>(refpart.beta_gamma());

            // gyromagnetic constant (used during spin push)
            m_gyro = refpart.gyromagnetic_anomaly;
            m_beta_gamma = refpart.beta_gamma();
            m_beta = refpart.beta();

        }

        /** This is an ExactSbend functor, so that a variable of this type can be used like an ExactSbend function.
         *
         * @param x particle position in x
         * @param y particle position in y
         * @param t particle position in t
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t
         * @param idcpu particle global index
         * @param refpart reference particle (unused)
         */
        template<typename T_Real=amrex::ParticleReal, typename T_IdCpu=uint64_t>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
            T_Real & AMREX_RESTRICT x,
            T_Real & AMREX_RESTRICT y,
            T_Real & AMREX_RESTRICT t,
            T_Real & AMREX_RESTRICT px,
            T_Real const & AMREX_RESTRICT py,
            T_Real const & AMREX_RESTRICT pt,
            T_IdCpu & AMREX_RESTRICT idcpu,
            [[maybe_unused]] RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using amrex::Math::powi;
            using namespace std; // for cmath(float)
            namespace stdx = amrex::simd::stdx;

            // access position data
            T_Real const xout = x;
            T_Real const yout = y;
            T_Real const tout = t;

            // initialize output values of momenta
            T_Real pxout = px;
            // T_Real pyout = py;
            // T_Real ptout = pt;

            // treat the special case of zero field
            if (m_phi==0_prt && m_B==0_prt) {
               // compute the radical in the denominator (= pz):
               T_Real const inv_pzden = 1_prt / sqrt(
                   powi<2>(pt - m_ibeta) -
                   m_ibetagam2 -
                   powi<2>(px) -
                   powi<2>(py)
               );

               // advance position and momentum (exact drift)
               x = xout + m_slice_ds * px * inv_pzden;
               // pxout = px;
               y = yout + m_slice_ds * py * inv_pzden;
               // pyout = py;
               t = tout - m_slice_ds * (m_ibeta + (pt - m_ibeta) * inv_pzden);
               // ptout = pt;

               // assign updated momenta
               // px = pxout;
               // py = pyout;
               // pt = ptout;

            } else {

               // assign intermediate quantities
               T_Real const pperp2 = powi<2>(pt)-2.0_prt * m_ibeta * pt - powi<2>(py)+1.0_prt;
               T_Real const px2 = powi<2>(px);
               // determine if particle lies within the domain of map definition
               auto const mask = pperp2 <= px2;
               amrex::ParticleIDWrapper<T_IdCpu>{idcpu}.make_invalid(mask);
               {
                   T_Real const pperp = sqrt(pperp2);
                   T_Real const pzi = sqrt(powi<2>(pperp) - powi<2>(px));
                   T_Real const rho = m_rc + xout;

                   // update momenta
                   pxout = px * m_cos_phi + (pzi - rho / m_rc) * m_sin_phi;
                   // pyout = py;
                   // ptout = pt;

                   // angle of momentum rotation
                   T_Real const px2f = powi<2>(pxout);
                   // determine if particle lies within domain of map definition
                   auto const mask2 = pperp2 <= px2f;
                   amrex::ParticleIDWrapper<T_IdCpu>{idcpu}.make_invalid(mask2);
                   {
                       T_Real const pzf = sqrt(powi<2>(pperp)-powi<2>(pxout));
                       T_Real const theta = m_slice_phi + asin(px/pperp) - asin(pxout/pperp);

                       // update position coordinates
                       x = -m_rc + rho*m_cos_phi + m_rc*(pzf + px*m_sin_phi - pzi*m_cos_phi);
                       y = yout + theta * m_rc * py;
                       t = tout - theta * m_rc * (pt - 1.0_prt * m_ibeta) - m_slice_phi * m_rc * m_ibeta;

                       // assign updated momenta
                       px = pxout;
                       // py = pyout;
                       // pt = ptout;
                   }
               }
            }
        }

        /** This pushes the reference particle.
         *
         * @param[in,out] refpart reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (RefPart & AMREX_RESTRICT refpart) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using amrex::Math::powi;

            // assign input reference particle values
            amrex::ParticleReal const x = refpart.x;
            amrex::ParticleReal const px = refpart.px;
            amrex::ParticleReal const y = refpart.y;
            amrex::ParticleReal const py = refpart.py;
            amrex::ParticleReal const z = refpart.z;
            amrex::ParticleReal const pz = refpart.pz;
            amrex::ParticleReal const t = refpart.t;
            amrex::ParticleReal const pt = refpart.pt;
            amrex::ParticleReal const s = refpart.s;

            // length of the current slice
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // special case of zero field = an exact drift
            if (m_phi==0_prt && m_B==0_prt) {
               // advance position and momentum (drift)
               amrex::ParticleReal const step = slice_ds /std::sqrt(powi<2>(pt)-1.0_prt);
               refpart.x = x + step*px;
               refpart.y = y + step*py;
               refpart.z = z + step*pz;
               refpart.t = t - step*pt;

            } else {

               // assign intermediate parameters
               amrex::ParticleReal const rc = (m_B != 0_prt) ? refpart.rigidity_Tm() / m_B : m_ds / m_phi;
               amrex::ParticleReal const B = refpart.beta_gamma() /rc;
               amrex::ParticleReal const bend_sign = m_ds * rc / std::abs(m_ds * rc);
               amrex::ParticleReal const theta = bend_sign * std::abs(m_phi) / nslice(); // sign is determined by rc*ds

               // calculate expensive terms once
               auto const [sin_theta, cos_theta] = amrex::Math::sincos(theta);

               // advance position and momentum (bend)
               refpart.px = px*cos_theta - pz*sin_theta;
               refpart.py = py;
               refpart.pz = pz*cos_theta + px*sin_theta;
               refpart.pt = pt;

               refpart.x = x + (refpart.pz - pz)/B;
               refpart.y = y + (theta/B)*py;
               refpart.z = z - (refpart.px - px)/B;
               refpart.t = t - (theta/B)*pt;

            }

            // advance integrated path length
            refpart.s = s + slice_ds;
        }

        /** This pushes the covariance matrix. */
        using LinearTransport::operator();

        /** This function returns the linear transport map.
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

            // length of the current slice
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // access reference particle values to find beta*gamma^2
            amrex::ParticleReal const pt_ref = refpart.pt;
            amrex::ParticleReal const betgam2 = powi<2>(pt_ref) - 1_prt;
            amrex::ParticleReal const bet = refpart.beta();
            amrex::ParticleReal const inv_bet = 1_prt / bet;

            // initialize linear map matrix elements
            Map6x6 R = Map6x6::Identity();

            // treat the special case of zero field
            if (m_phi == 0.0_prt && m_B == 0.0_prt)
            {
               R(1,2) = slice_ds;
               R(3,4) = slice_ds;
               R(5,6) = slice_ds / betgam2;
            } else
            {
               // reference particle's orbital radius
               amrex::ParticleReal const rc = this->rc(refpart);

               // calculate expensive terms once
               amrex::ParticleReal const theta = slice_ds / rc;
               auto const [sin_theta, cos_theta] = amrex::Math::sincos(theta);

               // assign linear map matrix elements
               R(1,1) = cos_theta;
               R(1,2) = rc * sin_theta;
               R(1,6) = -rc * inv_bet * (1_prt - cos_theta);
               R(2,1) = -sin_theta / rc;
               R(2,2) = cos_theta;
               R(2,6) = -sin_theta * inv_bet;
               R(3,4) = rc * theta;
               R(5,1) = sin_theta * inv_bet;
               R(5,2) = rc * inv_bet * (1_prt - cos_theta);
               R(5,6) = rc * (sin_theta * inv_bet * inv_bet - theta);
            }

            // apply the transverse rotation (roll) alignment error
            return rotate_aligned_map(R);
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
            T_Real & AMREX_RESTRICT x,
            T_Real & AMREX_RESTRICT y,
            T_Real & AMREX_RESTRICT t,
            T_Real & AMREX_RESTRICT px,
            T_Real const & AMREX_RESTRICT py,
            T_Real const & AMREX_RESTRICT pt,
            T_Real & AMREX_RESTRICT sx,
            T_Real & AMREX_RESTRICT sy,
            T_Real & AMREX_RESTRICT sz,
            T_IdCpu & AMREX_RESTRICT idcpu,
            RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using namespace std; // for cmath(float)

            // special case of zero field: exact drift, spin is unaffected
            if (m_phi == 0_prt && m_B == 0_prt) {
                (*this)(x, y, t, px, py, pt, idcpu, refpart);
                return;
            }

            // initialize the three components of the axis-angle vector
            T_Real lambdax = 0_prt;
            T_Real lambday = 0_prt;
            T_Real lambdaz = 0_prt;

            // store initial variables needed for spin push
            T_Real const ti = t;
            T_Real const Pnorm = sqrt(1_prt + pt*pt - 2_prt*pt*m_ibeta);
            T_Real const gamma = (1_prt - pt*m_beta)*refpart.gamma();
            T_Real const ux = px/Pnorm;
            T_Real const uy = py/Pnorm;
            T_Real const uz = sqrt(1_prt - ux*ux - uy*uy);

            // For this element, note that the time-of-flight from entry to
            // exit is needed for evaluation of the spin map, so the phase
            // space push is applied first, in order to obtain the change in t.

            // phase space push
            (*this)(x, y, t, px, py, pt, idcpu, refpart);

            // particle transit time through the slice
            T_Real const dt = t - ti + m_slice_phi * m_rc * m_ibeta;

            // advance in cyclotron phase
            T_Real const theta = dt / m_rc * m_beta_gamma / gamma;

            lambdax = m_gyro * theta * (gamma - 1_prt) * uy * ux;
            lambday = m_gyro * theta * ( -gamma + (gamma - 1_prt) * uy * uy );
            lambdaz = m_gyro * theta * (gamma - 1_prt) * uy * uz;

            // push the spin vector using the generator just determined (first factor)
            rotate_spin(lambdax,lambday,lambdaz,sx,sy,sz);

            lambdax = 0_prt;
            lambday = m_slice_phi - theta;
            lambdaz = 0_prt;

            // push the spin vector using the generator just determined (second factor)
            rotate_spin(lambdax,lambday,lambdaz,sx,sy,sz);

        }


        amrex::ParticleReal m_phi; //! bend angle in radians
        amrex::ParticleReal m_B;  //! magnetic field in T

      private:
        // constants that are independent of the individually tracked particle,
        // see: compute_constants() to refresh
        amrex::ParticleReal m_slice_ds, m_slice_phi, m_ibeta, m_ibetagam2, m_rc, m_sin_phi, m_cos_phi;
        amrex::ParticleReal m_gyro, m_beta, m_beta_gamma;

    };

    /** ExactSbend: GPU state plus common host-only metadata */
    using ExactSbend = Element<ExactSbendGPU>;

} // namespace impactx

IMPACTX_PUSH_EXTERN_TEMPLATE(impactx::elements::ExactSbendGPU)

#endif // IMPACTX_EXACTSBEND_H
