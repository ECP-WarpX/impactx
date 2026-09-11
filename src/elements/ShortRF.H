/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_SHORTRF_H
#define IMPACTX_SHORTRF_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/gpuelement.H"
#include "mixin/alignment.H"
#include "mixin/beamoptic.H"
#include "mixin/thin.H"
#include "mixin/nofinalize.H"
#include "mixin/pipeaperture.H"
#include "mixin/lineartransport.H"
#include "mixin/spintransport.H"

#include <AMReX_Extension.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>
#include <AMReX_SIMD.H>

#include <cmath>

namespace impactx::elements
{
    struct ShortRFGPU
    : public mixin::BeamOptic<ShortRFGPU>,
      public mixin::LinearTransport<ShortRFGPU>,
      public mixin::Thin,
      public mixin::Alignment,
      public mixin::PipeAperture,
      public mixin::SpinTransport,
      public mixin::NoFinalize,
      public amrex::simd::Vectorized<amrex::simd::native_simd_size_particlereal>
    {
        static constexpr auto type = "ShortRF";
        using PType = ImpactXParticleContainer::ParticleType;

        static constexpr amrex::ParticleReal DEFAULT_phase = -90;

        /** A short RF cavity element
         *
         * @param V Normalized RF voltage V = maximum energy gain/(m*c^2)
         * @param freq RF frequency in Hz
         * @param phase RF synchronous phase of the reference particle in deg
         *     phi = 0:  maximum energy gain (on-crest)
         *     phi = -90 deg:  zero-crossing for bunching
         *     phi = 90 deg:  zero-crossing for debunching
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param aperture_x horizontal half-aperture in m
         * @param aperture_y vertical half-aperture in m
         */
        ShortRFGPU (
            amrex::ParticleReal V,
            amrex::ParticleReal freq,
            amrex::ParticleReal phase = DEFAULT_phase,
            amrex::ParticleReal dx = DEFAULT_dx,
            amrex::ParticleReal dy = DEFAULT_dy,
            amrex::ParticleReal rotation_degree = DEFAULT_rotation_degree,
            amrex::ParticleReal aperture_x = DEFAULT_aperture_x,
            amrex::ParticleReal aperture_y = DEFAULT_aperture_y
        )
        : Alignment(dx, dy, rotation_degree),
          PipeAperture(aperture_x, aperture_y),
          m_V(V), m_freq(freq), m_phase(phase)
        {
        }

        /** Reverse the element (negate V) */
        void reverse () { m_V = -m_V; }

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

            // Define parameters and intermediate constants
            using ablastr::constant::math::pi;
            using ablastr::constant::SI::c;
            m_k = 2.0_prt * pi / c * m_freq;
            m_phi = m_phase * pi / 180.0_prt;

            // access reference particle values (final, initial):
            amrex::ParticleReal const ptf_ref = refpart.pt;
            m_V_cos_phi = m_V * std::cos(m_phi);
            amrex::ParticleReal const pti_ref = ptf_ref + m_V_cos_phi;
            m_bgf = std::sqrt(powi<2>(ptf_ref) - 1.0_prt);
            m_bgi = std::sqrt(powi<2>(pti_ref) - 1.0_prt);
        }

        /** This is a shortrf functor, so that a variable of this type can be used like a
         *  shortrf function.
         *
         * The @see compute_constants method must be called before pushing particles through this operator.
         *
         * @param x particle position in x (unused)
         * @param y particle position in y (unused)
         * @param t particle position in t
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t
         * @param idcpu particle global index (unused)
         * @param refpart reference particle (unused)
         */
        template<typename T_Real=amrex::ParticleReal, typename T_IdCpu=uint64_t>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
            [[maybe_unused]] T_Real const & AMREX_RESTRICT x,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT y,
            T_Real const & AMREX_RESTRICT t,
            T_Real & AMREX_RESTRICT px,
            T_Real & AMREX_RESTRICT py,
            T_Real & AMREX_RESTRICT pt,
            [[maybe_unused]] T_IdCpu const & AMREX_RESTRICT idcpu,
            [[maybe_unused]] RefPart const & AMREX_RESTRICT refpart
        ) const
        {

            using namespace amrex::literals; // for _rt and _prt
            using namespace std; // for cmath(float)

            // initial conversion from static to dynamic units:
            px = px * m_bgi;
            py = py * m_bgi;
            pt = pt * m_bgi;

            // initialize output values
            // T_Real xout = x;
            // T_Real yout = y;
            // T_Real tout = t;
            // T_Real pxout = px;
            // T_Real pyout = py;
            T_Real ptout = pt;

            // advance position and momentum in dynamic units
            // xout = x;
            // pxout = px;

            // yout = y;
            // pyout = py;

            // tout = t;
            ptout = pt - m_V * cos(m_k * t + m_phi) + m_V_cos_phi;

            // assign updated values
            // x = xout;
            // y = yout;
            // t = tout;
            // px = pxout;
            // py = pyout;
            pt = ptout;

            // final conversion from dynamic to static units:
            px = px / m_bgf;
            py = py / m_bgf;
            pt = pt / m_bgf;
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

            // Define parameters and intermediate constants
            using ablastr::constant::math::pi;
            amrex::ParticleReal const phi = m_phase*(pi/180.0_prt);

            // compute initial value of beta*gamma
            amrex::ParticleReal const bgi = std::sqrt(powi<2>(pt) - 1.0_prt);

            // advance pt
            refpart.pt = pt - m_V * std::cos(phi);

            // compute final value of beta*gamma
            amrex::ParticleReal const ptf = refpart.pt;
            amrex::ParticleReal const bgf = std::sqrt(powi<2>(ptf) - 1.0_prt);

            // advance position (x,y,z,t)
            refpart.x = x;
            refpart.y = y;
            refpart.z = z;
            refpart.t = t;

            // advance momentum (px,py,pz)
            refpart.px = px*bgf/bgi;
            refpart.py = py*bgf/bgi;
            refpart.pz = pz*bgf/bgi;

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
            T_Real const & AMREX_RESTRICT t,
            T_Real & AMREX_RESTRICT px,
            T_Real & AMREX_RESTRICT py,
            T_Real & AMREX_RESTRICT pt,
            T_Real & AMREX_RESTRICT sx,
            T_Real & AMREX_RESTRICT sy,
            T_Real & AMREX_RESTRICT sz,
            T_IdCpu const & AMREX_RESTRICT idcpu,
            RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using namespace std; // for cmath(float)

            // reference particle relativistic factors, before and after the
            // cavity: the reference particle carries pt = -gamma of the final
            // state, and gained m_V_cos_phi passing through
            amrex::ParticleReal const gyro_anomaly = refpart.gyromagnetic_anomaly;
            amrex::ParticleReal const gamf = -refpart.pt;
            amrex::ParticleReal const gami = gamf - m_V_cos_phi;

            // initialize the three components of the axis-angle vector
            T_Real lambdax = 0_prt;
            T_Real lambday = 0_prt;
            T_Real lambdaz = 0_prt;

            // The ShortRF model is based on the impulse associated with the
            // on-axis longitudinal electric field Ez.  As a result, only this
            // contribution is considered below.

            // Integrated magnetic field normalized by q/mc (half-step):
            T_Real const Bx = 0_prt;
            T_Real const By = 0_prt;
            T_Real const Bz = 0_prt;

            // Integrated electric field normalized by q/mc^2 (half-step):
            T_Real const Ex = 0.0_prt;
            T_Real const Ey = 0.0_prt;
            T_Real const Ez = 0.5_prt * m_V * cos(m_k * t + m_phi);

            // Quantities required to evaluate the full Thomas-BMT precession vector (before phase space kick).
            amrex::ParticleReal beta = m_bgi / gami;
            amrex::ParticleReal ibeta = gami / m_bgi;
            T_Real Pnorm = sqrt(1_prt - 2_prt * pt * ibeta + pt*pt);
            T_Real iPnorm = 1_prt/Pnorm;
            T_Real gamma = gami * (1_prt - pt*beta);
            T_Real Ps = sqrt(Pnorm*Pnorm - px*px - py*py);
            T_Real ux = px * iPnorm;
            T_Real uy = py * iPnorm;
            T_Real uz = Ps * iPnorm;

            // Zero curvature
            amrex::ParticleReal const h = 0.0_prt;

            // Evaluation of the full Thomas-BMT precession vector.
            tbmt_precession_vector(x,ux,uy,uz,gamma,h,gyro_anomaly,Bx,By,Bz,Ex,Ey,Ez,lambdax,lambday,lambdaz);

            // push the spin vector using the generator just determined (half spin-kick)
            rotate_spin(lambdax, lambday, lambdaz, sx, sy, sz);

            // phase space push
            (*this)(x, y, t, px, py, pt, idcpu, refpart);

            // Quantities required to evaluate the full Thomas-BMT precession vector (after phase space kick).
            beta = m_bgf / gamf;
            ibeta = gamf / m_bgf;
            Pnorm = sqrt(1_prt - 2_prt * pt * ibeta + pt*pt);
            iPnorm = 1_prt/Pnorm;
            gamma = gamf * (1_prt - pt*beta);
            Ps = sqrt(Pnorm*Pnorm - px*px - py*py);
            ux = px * iPnorm;
            uy = py * iPnorm;
            uz = Ps * iPnorm;

            // Evaluation of the full Thomas-BMT precession vector.
            tbmt_precession_vector(x,ux,uy,uz,gamma,h,gyro_anomaly,Bx,By,Bz,Ex,Ey,Ez,lambdax,lambday,lambdaz);

            // push the spin vector using the generator just determined (half spin-kick)
            rotate_spin(lambdax, lambday, lambdaz, sx, sy, sz);
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
        transport_map ([[maybe_unused]] RefPart const & AMREX_RESTRICT refpart) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using amrex::Math::powi;

            // Define parameters and intermediate constants
            using ablastr::constant::math::pi;
            using ablastr::constant::SI::c;
            amrex::ParticleReal const k = (2.0_prt*pi/c)*m_freq;
            amrex::ParticleReal const phi = m_phase*(pi/180.0_prt);

            // access reference particle values (final, initial):
            amrex::ParticleReal const ptf_ref = refpart.pt;
            amrex::ParticleReal const pti_ref = ptf_ref + m_V * std::cos(phi);
            amrex::ParticleReal const bgf = std::sqrt(powi<2>(ptf_ref) - 1.0_prt);
            amrex::ParticleReal const bgi = std::sqrt(powi<2>(pti_ref) - 1.0_prt);

            // initialize transport matrix
            Map6x6 R = Map6x6::Identity();

            // This is a nonlinear element: the linearized map is returned here.
            // assign linear map matrix elements
            R(2,2) = bgi/bgf;
            R(4,4) = bgi/bgf;
            R(6,5) = k*m_V*std::sin(phi)/bgf;
            R(6,6) = bgi/bgf;

            // apply the transverse rotation (roll) alignment error
            return rotate_aligned_map(R);
        }

        amrex::ParticleReal m_V; //! normalized (max) RF voltage drop.
        amrex::ParticleReal m_freq; //! RF frequency in Hz.
        amrex::ParticleReal m_phase; //! reference RF phase in degrees.

      private:
        // constants that are independent of the individually tracked particle,
        // see: compute_constants() to refresh
        amrex::ParticleReal m_k, m_phi, m_V_cos_phi, m_bgi, m_bgf;
    };

    /** ShortRF: GPU state plus common host-only metadata */
    using ShortRF = Element<ShortRFGPU>;

} // namespace impactx

IMPACTX_PUSH_EXTERN_TEMPLATE(impactx::elements::ShortRFGPU)

#endif // IMPACTX_SHORTRF_H
