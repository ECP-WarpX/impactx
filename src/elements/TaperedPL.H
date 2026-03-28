/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_TAPEREDPL_H
#define IMPACTX_TAPEREDPL_H

#include "particles/ImpactXParticleContainer.H"
#include "Element.H"
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
    struct TaperedPLPhysics
    : public mixin::BeamOptic<TaperedPLPhysics>,
      public mixin::LinearTransport<TaperedPLPhysics>,
      public mixin::Thin,
      public mixin::Alignment,
      public mixin::PipeAperture,
      public mixin::SpinTransport,
      public mixin::NoFinalize,
      // At least on Intel AVX512, there is a small overhead to vectorize this element
      // for pure phase-space pushes. But a win for spin-tracking. See:
      //   https://github.com/BLAST-ImpactX/impactx/pull/1002
      //   https://github.com/BLAST-ImpactX/impactx/pull/1499
      public amrex::simd::Vectorized<amrex::simd::native_simd_size_particlereal>
    {
        static constexpr auto type = "TaperedPL";
        using PType = ImpactXParticleContainer::ParticleType;

        static constexpr int DEFAULT_unit = 0;

        /** A short segment of a nonlinear plasma lens with a transverse taper.
         *
         * @param k  integrated plasma lens (linear) focusing strength in m^(-1)
         *           = (length in m) * (azimuthal field gradient in T/m) / (rigidity in T-m)
         *        OR integrated azimuthal magnetic field gradient in T (k > 0)
         * @param unit  Unit specification
         *           unit = 0 integrated linear focusing strength in m^(-1)
         *           unit = 1 integrated linear focusing strength in T
         * @param taper  horizontal taper parameter in m^(-1)
                    = 1 / (horizontal dispersion in m)
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param aperture_x horizontal half-aperture in m
         * @param aperture_y vertical half-aperture in m
         */
        TaperedPLPhysics (
            amrex::ParticleReal k,
            amrex::ParticleReal taper,
            int unit = DEFAULT_unit,
            amrex::ParticleReal dx = DEFAULT_dx,
            amrex::ParticleReal dy = DEFAULT_dy,
            amrex::ParticleReal rotation_degree = DEFAULT_rotation_degree,
            amrex::ParticleReal aperture_x = DEFAULT_aperture_x,
            amrex::ParticleReal aperture_y = DEFAULT_aperture_y
        )
        : Alignment(dx, dy, rotation_degree),
          PipeAperture(aperture_x, aperture_y),
          m_k(k), m_taper(taper), m_unit(unit)
        {
        }

        /** Reverse the element (negate strength) */
        void reverse () { m_k = -m_k; }

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

            // normalize focusing strength units to MAD-X convention if needed
            m_g = m_unit == 1 ? m_k / refpart.rigidity_Tm() : m_k;

            // additional constants needed for spin push
            m_beta = refpart.beta();
            m_ibeta = 1_prt / m_beta;

        }

        /** Push all particles */
        using BeamOptic::operator();

        /** This is a tapered PL functor, so that a variable of this type can be used like a
         *  shortrf function.
         *
         * @param x particle position in x
         * @param y particle position in y
         * @param t particle position in t (unused)
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t (unused)
         * @param idcpu particle global index (unused)
         * @param refpart reference particle (unused)
         */
        template<typename T_Real=amrex::ParticleReal, typename T_IdCpu=uint64_t>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
            T_Real const & AMREX_RESTRICT x,
            T_Real const & AMREX_RESTRICT y,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT t,
            T_Real & AMREX_RESTRICT px,
            T_Real & AMREX_RESTRICT py,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT pt,
            [[maybe_unused]] T_IdCpu const & AMREX_RESTRICT idcpu,
            [[maybe_unused]] RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using amrex::Math::powi;
            using namespace std; // for cmath(float)

            // initialize output values
            // T_Real xout = x;
            // T_Real yout = y;
            // T_Real tout = t;
            T_Real pxout = px;
            T_Real pyout = py;
            // T_Real ptout = pt;

            // advance position and momentum
            // xout = x;
            pxout = px - m_g * ( x + m_taper*0.5_prt * (powi<2>(x) + powi<2>(y)) );

            // yout = y;
            pyout = py - m_g * ( y + m_taper * x * y );

            // tout = t;
            // ptout = pt;

            // assign updated values
            // x = xout;
            // y = yout;
            // t = tout;
            px = pxout;
            py = pyout;
            // pt = ptout;
        }

        /** This pushes the reference particle. */
        using Thin::operator();


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
            using amrex::Math::powi;

            // The Thomas-BMT generator is evaluated at the transverse momentum
            // halfway through the kick, so that reversing the element (negating
            // the strength) produces the exact inverse rotation.
            T_Real const px_in = px;
            T_Real const py_in = py;

            // pure phase-space push (shared with the non-spin operator())
            (*this)(x, y, t, px, py, pt, idcpu, refpart);

            // transverse momenta at the kick midpoint
            T_Real const px_mid = 0.5_prt * (px_in + px);
            T_Real const py_mid = 0.5_prt * (py_in + py);

            // initialize the three components of the axis-angle vector
            T_Real lambdax = 0_prt;
            T_Real lambday = 0_prt;
            T_Real lambdaz = 0_prt;

            // compute multipole magnetic field normalized by q/mc:
            amrex::ParticleReal const gamma_ref = refpart.gamma();
            T_Real const Bx = -m_g * ( y + m_taper * x * y ) * gamma_ref * m_beta;
            T_Real const By = m_g * ( x + m_taper*0.5_prt * (powi<2>(x) + powi<2>(y)) ) * gamma_ref * m_beta;
            T_Real const Bz = 0_prt;

            // Electric field normalized by q/mc^2:
            T_Real const Ex = 0.0_prt;
            T_Real const Ey = 0.0_prt;
            T_Real const Ez = 0.0_prt;

            // Quantities required to evaluate the full Thomas-BMT precession
            // vector, at the midpoint transverse momentum.
            T_Real const Pnorm = sqrt(1_prt - 2_prt * pt * m_ibeta + pt*pt);
            T_Real const iPnorm = 1_prt/Pnorm;
            T_Real const Ps = sqrt(Pnorm*Pnorm - px_mid*px_mid - py_mid*py_mid);
            T_Real const gamma = gamma_ref * (1_prt - pt*m_beta);
            T_Real const ux = px_mid * iPnorm;
            T_Real const uy = py_mid * iPnorm;
            T_Real const uz = Ps * iPnorm;

            // Zero curvature in a quadrupole
            amrex::ParticleReal const h = 0.0_prt;

            // Evaluation of the full Thomas-BMT precession vector.
            amrex::ParticleReal const gyro_anomaly = refpart.gyromagnetic_anomaly;
            tbmt_precession_vector(x, ux, uy, uz, gamma, h, gyro_anomaly, Bx, By, Bz, Ex, Ey, Ez, lambdax, lambday, lambdaz);

            // push the spin vector using the generator just determined
            rotate_spin(lambdax, lambday, lambdaz, sx, sy, sz);

        }


        /** This pushes the covariance matrix. */
        using LinearTransport::operator();

        /** This function returns the linear transport map.
         *
         * The tapered plasma-lens kick implemented by @see operator() is
         *   dpx = -g * ( x + taper/2 * (x^2 + y^2) )
         *   dpy = -g * ( y + taper * x * y )
         * where g is the integrated focusing strength (in units that
         * match @c m_unit: either m^-1 directly, or rigidity-normalized
         * from T via refpart.rigidity_Tm()).
         *
         * Linearized at the reference trajectory (x = y = 0), the taper
         * contributes only quadratic terms and drops out. The element
         * therefore reduces to a thin axisymmetric focusing kick:
         *   R(2,1) = -g,  R(4,3) = -g.
         *
         * Ideal behavior only: alignment and rotation errors, and any
         * feed-down from an off-axis reference orbit, are ignored; the
         * map is taken at (x, y) = (0, 0).
         *
         * MAD-X convention: no direct analogue. Its linear limit
         * coincides with a thin axisymmetric focusing quadrupole-like
         * kick (equivalent to a MAD-X MULTIPOLE with normal KNL of
         * equal strength in both transverse planes), but MAD-X has no
         * element that carries the quadratic taper term.
         *
         * Physics references:
         * - C. Lindstrøm et al., "Achromatic optics using nonlinear plasma lenses
for beam-quality preservation between plasma-accelerator stages",
         *   submitted (2026), https://arxiv.org/abs/2604.17605
         *
         * @param[in] refpart reference particle (used only for the
         *                    unit conversion when m_unit == 1)
         * @returns 6x6 transport matrix
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        Map6x6
        transport_map (RefPart const & AMREX_RESTRICT refpart) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // Normalize focusing strength to MAD-X convention if needed.
            amrex::ParticleReal const g = (m_unit == 1)
                ? m_k / refpart.rigidity_Tm()
                : m_k;

            // Linearized axisymmetric focusing kick; taper is quadratic
            // and therefore absent from the linear map.
            Map6x6 R = Map6x6::Identity();
            R(2,1) = -g;
            R(4,3) = -g;

            // apply the transverse rotation (roll) alignment error
            return rotate_aligned_map(R);
        }

        amrex::ParticleReal m_k; //! linear focusing strength (field gradient)
        amrex::ParticleReal m_taper; //! horizontal taper parameter
        int m_unit;  //! units for linear focusing strength (field gradient)

      private:
        // constants that are independent of the individually tracked particle
        amrex::ParticleReal m_beta;           //! relativistic beta
        amrex::ParticleReal m_ibeta;          //! inverse of beta
        amrex::ParticleReal m_g;              //! quadrupole strength
    };

    /** TaperedPL: physics plus common host-only metadata */
    using TaperedPL = Element<TaperedPLPhysics>;

} // namespace impactx

IMPACTX_PUSH_EXTERN_TEMPLATE(impactx::elements::TaperedPLPhysics)

#endif // IMPACTX_TAPEREDPL_H
