/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_DIPEDGE_H
#define IMPACTX_DIPEDGE_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/gpuelement.H"
#include "mixin/alignment.H"
#include "mixin/beamoptic.H"
#include "mixin/thin.H"
#include "mixin/lineartransport.H"
#include "mixin/nofinalize.H"
#include "mixin/pipeaperture.H"
#include "mixin/spintransport.H"

#include <ablastr/constant.H>

#include <AMReX_Enum.H>
#include <AMReX_Extension.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>
#include <AMReX_SIMD.H>

#include <cmath>
#include <stdexcept>

namespace impactx::elements
{
    struct DipEdgeGPU
    : public mixin::BeamOptic<DipEdgeGPU>,
      public mixin::LinearTransport<DipEdgeGPU>,
      public mixin::Thin,
      public mixin::Alignment,
      public mixin::PipeAperture,
      public mixin::SpinTransport,
      public mixin::NoFinalize,
      public amrex::simd::Vectorized<amrex::simd::native_simd_size_particlereal>
    {
        static constexpr auto type = "DipEdge";
        using PType = ImpactXParticleContainer::ParticleType;

        AMREX_ENUM_IN_CLASS(Model,
            linear,   //! for the linear fringe field model
            nonlinear, //! for the 6th order fringe field model
        );

        AMREX_ENUM_IN_CLASS(Location,
            entry,  //! for leading (entry) fringe field
            exit,   //! for trailing (exit) fringe field
        );

        static constexpr amrex::ParticleReal DEFAULT_R = 1;
        static constexpr amrex::ParticleReal DEFAULT_K0 = ablastr::constant::math::pi *
            ablastr::constant::math::pi / 6.0;
        static constexpr amrex::ParticleReal DEFAULT_K1 = 0;
        static constexpr amrex::ParticleReal DEFAULT_K2 = 1;
        static constexpr amrex::ParticleReal DEFAULT_K3 = amrex::ParticleReal(1.0 / 6.0);
        static constexpr amrex::ParticleReal DEFAULT_K4 = 0;
        static constexpr amrex::ParticleReal DEFAULT_K5 = 0;
        static constexpr amrex::ParticleReal DEFAULT_K6 = 0;
        static constexpr Model DEFAULT_model = Model::linear;
        static constexpr Location DEFAULT_location = Location::entry;
        static constexpr bool DEFAULT_modify_ref_part = false;

        /** Edge focusing associated with bend entry or exit
         *
         * This model assumes a first-order effect of nonzero gap.
         * Here we use the linear fringe field map, given to first order in g/rc (gap / radius of curvature).
         *
         * References:
         *   K. L. Brown, SLAC Report No. 75 (1982).
         *   K. Hwang and S. Y. Lee, PRAB 18, 122401 (2015).
         *
         * @param psi Pole face entry/exit angle in rad.
         * @param rc Radius of curvature in m.
         * @param g Gap parameter in m.
         * @param R Length scale parameter in m.
         * @param K0 Fringe field integeral (unitless).
         * @param K1 Fringe field integeral (unitless).
         * @param K2 Fringe field integral (unitless).
         * @param K3 Fringe field integral (unitless).
         * @param K4 Fringe field integral (unitless).
         * @param K5 Fringe field integral (unitless).
         * @param K6 Fringe field integral (unitless).
         * @param model Fringe field model (linear or nonlinear).
         * @param location Fringe field location (entry or exit).
         * @param modify_ref_part Apply fringe field to the reference particle (boolean).
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param aperture_x horizontal half-aperture in m
         * @param aperture_y vertical half-aperture in m
         */
        DipEdgeGPU (
            amrex::ParticleReal psi,
            amrex::ParticleReal rc,
            amrex::ParticleReal g,
            amrex::ParticleReal R = DEFAULT_R,
            amrex::ParticleReal K0 = DEFAULT_K0,
            amrex::ParticleReal K1 = DEFAULT_K1,
            amrex::ParticleReal K2 = DEFAULT_K2,
            amrex::ParticleReal K3 = DEFAULT_K3,
            amrex::ParticleReal K4 = DEFAULT_K4,
            amrex::ParticleReal K5 = DEFAULT_K5,
            amrex::ParticleReal K6 = DEFAULT_K6,
            Model model = DEFAULT_model,
            Location location = DEFAULT_location,
            bool modify_ref_part = DEFAULT_modify_ref_part,
            amrex::ParticleReal dx = DEFAULT_dx,
            amrex::ParticleReal dy = DEFAULT_dy,
            amrex::ParticleReal rotation_degree = DEFAULT_rotation_degree,
            amrex::ParticleReal aperture_x = DEFAULT_aperture_x,
            amrex::ParticleReal aperture_y = DEFAULT_aperture_y
        )
        : Alignment(dx, dy, rotation_degree),
          PipeAperture(aperture_x, aperture_y),
          m_psi(psi), m_rc(rc), m_g(g), m_R(R), m_K0(K0), m_K1(K1), m_K2(K2), m_K3(K3), m_K4(K4), m_K5(K5), m_K6(K6), m_model(model), m_location(location), m_modify_ref_part(modify_ref_part)
        {
        }

        /** Reverse the element (negate psi and K2, swap entry/exit)
         *
         * @todo the m_modify_ref_part=true parameter is not yet handled
         *       see https://github.com/BLAST-ImpactX/impactx/pull/1431
         */
        void reverse ()
        {
            m_psi = -m_psi;
            m_K2 = -m_K2;
            m_location = (m_location == Location::entry)
                ? Location::exit
                : Location::entry;
        }

        /** Push all particles */
        using BeamOptic::operator();

        /** Compute and cache the constants for the push.
         *
         * In particular, used to pre-compute and cache variables that are
         * independent of the individually tracked particle.
         *
         * @param refpart reference particle (unused)
         */
        void compute_constants ([[maybe_unused]] RefPart const & refpart)
        {
            using namespace amrex::literals; // for _rt and _prt
            using amrex::Math::powi;

            Alignment::compute_constants(refpart);

            // constants for linear map:

            // edge focusing matrix elements (zero gap)
            m_R21 = std::tan(m_psi) / m_rc;

            // first-order effect of nonzero gap
            auto const [sin_psi, cos_psi] = amrex::Math::sincos(m_psi);
            amrex::ParticleReal const vf = (1.0_prt + powi<2>(sin_psi))
                / powi<3>(cos_psi)
                * m_g * m_K2 / powi<2>(m_rc);

            m_R43 = vf - m_R21;

            // constants for nonlinear map:

            // access reference particle values to find beta
            m_beta = refpart.beta();

            // expressions from eq (14) of Mitchell and Hwang, NAPAC2025
            amrex::ParticleReal const tan_psi = sin_psi/cos_psi;
            amrex::ParticleReal const sec_psi = 1_prt/cos_psi;

            // the scale constant specifying entry or exit
            amrex::ParticleReal const loc = m_location == Location::entry ? 1_prt : -1_prt;

            m_c1 = (m_g / m_rc) * m_K1 / cos_psi;
            m_c2_times_1plusdelta = powi<2>(m_g/m_rc) * m_K0 * sin_psi * powi<3>(sec_psi)/2_prt;
            m_c3_times_1plusdelta = loc * powi<2>(m_g)/m_rc * m_K0 * powi<2>(sec_psi);
            m_c4_times_1plusdelta = loc * (m_g / m_rc) * m_K1 * sin_psi / (powi<2>(cos_psi));
            m_c5 = tan_psi / (2_prt * m_rc);
            m_c6_times_1plusdelta = 0.5_prt * powi<2>(sin_psi)/(2_prt * m_rc * powi<3>(cos_psi)) * m_g/m_rc * m_K1;
            m_c7_times_1plusdelta = 0.5_prt * powi<3>(sec_psi) * (m_g * m_K1 / (2_prt*powi<2>(m_rc)) + (1_prt+powi<2>(sin_psi))*m_g*m_K2/powi<2>(m_rc));
            m_c8_times_1plusdelta = (1_prt/6_prt) * powi<3>(tan_psi) / (2_prt * powi<2>(m_rc));
            m_c9_times_1plusdelta = 0.5_prt * tan_psi * powi<2>(sec_psi) / (2_prt * powi<2>(m_rc));
            m_c10_times_1plusdelta = loc * powi<2>(tan_psi) / (2_prt * m_rc);
            m_c11_times_1plusdelta = loc * 1_prt / (2_prt * m_rc);
            m_c12_times_1plusdelta = (m_g > 0_prt) ? 1_prt / (24_prt) * (4_prt/cos_psi - 8_prt/(powi<3>(cos_psi))) * m_K3/(powi<2>(m_rc) * m_g) : 0_prt;
            m_c13 = powi<2>(sin_psi) / (2_prt*powi<3>(cos_psi)) * powi<2>(m_g)/(m_rc*m_R) * m_K4;
            m_c14 = 0.5_prt * sin_psi/(powi<3>(cos_psi)) * m_g/(m_rc*m_R) * m_K5;
            m_c15 = (m_K6 / (m_rc*m_R)) / (powi<3>(cos_psi));
        }

        /** This is a dipedge functor, so that a variable of this type can be used like a
         *  dipedge function.
         *
         * The @see compute_constants method must be called before pushing particles through this operator.
         *
         * @param x particle position in x
         * @param y particle position in y
         * @param t particle position in t (unchanged by the linear model)
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t (unchanged)
         * @param idcpu particle global index (unused)
         * @param refpart reference particle (unused)
         */
        template<typename T_Real=amrex::ParticleReal, typename T_IdCpu=uint64_t>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
            T_Real & AMREX_RESTRICT x,
            T_Real & AMREX_RESTRICT y,
            T_Real & AMREX_RESTRICT t,
            T_Real & AMREX_RESTRICT px,
            T_Real & AMREX_RESTRICT py,
            T_Real const & AMREX_RESTRICT pt,
            [[maybe_unused]] T_IdCpu const & AMREX_RESTRICT idcpu,
            [[maybe_unused]] RefPart const & AMREX_RESTRICT refpart
        ) const
        {

            using namespace amrex::literals; // for _rt and _prt
            using namespace std; // for cmath(float)
            using amrex::Math::powi;

            if (m_model == Model::linear) {
               // apply linear model
               px = px + m_R21 * x;
               py = py + m_R43 * y;
            } else {
               // apply nonlinear model
               T_Real xout = x;
               T_Real pxout = px;
               T_Real yout = y;
               T_Real pyout = py;
               T_Real tout = t;
               // T_Real ptout = pt;

               // compute particle momentum deviation delta + 1 (reciprocal)
               T_Real const inv_delta1 = 1_prt / sqrt(1_prt - 2_prt*pt/m_beta + powi<2>(pt));

               T_Real const c2 = m_c2_times_1plusdelta * inv_delta1;
               T_Real const c3 = m_c3_times_1plusdelta * inv_delta1;
               //T_Real const c4 = m_c4_times_1plusdelta * inv_delta1; //not currently used
               //T_Real const c6 = m_c6_times_1plusdelta * inv_delta1; //not currently used
               T_Real const c7 = m_c7_times_1plusdelta * inv_delta1;
               T_Real const c8 = m_c8_times_1plusdelta * inv_delta1;
               T_Real const c9 = m_c9_times_1plusdelta * inv_delta1;
               T_Real const c10 = m_c10_times_1plusdelta * inv_delta1;
               T_Real const c11 = m_c11_times_1plusdelta * inv_delta1;
               T_Real const c12 = m_c12_times_1plusdelta * inv_delta1;

               // update transverse coordinates
               xout = x - c3 - c10*powi<2>(x) + (c10 + c11)*powi<2>(y);
               yout = y + 2_prt * c10 * x * y;

               // update transverse momenta
               T_Real const D = 1_prt - 4_prt*c10*c11*powi<2>(y) - 4_prt*powi<2>(c10)*(powi<2>(x)+powi<2>(y));
               T_Real const dRdx = -m_c13 + c2 + c3*m_c5 + 2_prt*(m_c14 - m_c5)*x + 3_prt*(m_c15/6_prt + c10*m_c5 + c8)*powi<2>(x)
                                   + (-m_c15/2_prt + c10*m_c5 - c11*m_c5 - c9)*powi<2>(y);
               T_Real const dRdy = 2_prt*(-m_c14 + m_c5 - c7)*y + 2_prt*(-m_c15/2_prt + c10*m_c5 - c11*m_c5 - c9)*x*y - 4_prt*c12*powi<3>(y);
               T_Real const dXdx = 1_prt - 2_prt*c10*x;
               T_Real const dXdy = 2_prt*(c10 + c11)*y;
               T_Real const dYdx = 2_prt*c10*y;
               T_Real const dYdy = 1_prt + 2_prt*c10*x;

               pxout = (dYdy * (px - dRdx) - dYdx * (py - dRdy)) / D;
               pyout = (-dXdy * (px - dRdx) + dXdx * (py - dRdy)) / D;

               // update temporal coordinate
               T_Real const dDelta_dpt = (pt - 1_prt/m_beta) * inv_delta1;
               tout = t - dDelta_dpt * inv_delta1 * ((c2 + c3*m_c5)*x + (c10*m_c5 + c8)*powi<3>(x) - c7*powi<2>(y) - c12*powi<4>(y)
                      + (c10*m_c5 - c11*m_c5 - c9)*x*powi<2>(y) + pxout*(-c3 + (c11 + c10)*powi<2>(y) - c10*powi<2>(x)) + 2_prt*pyout*c10*x*y);

               // optional removal of reference particle shift
               amrex::ParticleReal const c2_ref = m_c2_times_1plusdelta;
               amrex::ParticleReal const c3_ref = m_c3_times_1plusdelta;
               T_Real const dx_ref = (m_modify_ref_part) ? -c3_ref : 0_prt;
               T_Real const dpx_ref = (m_modify_ref_part) ? m_c13 - c2_ref - c3_ref*m_c5 : 0_prt;
               T_Real const dt_ref = -c3_ref * dpx_ref / m_beta;
               xout -= dx_ref;
               pxout -= dpx_ref;
               tout -= dt_ref;

               x = xout;
               px = pxout;
               y = yout;
               py = pyout;
               t = tout;
               // pt = ptout;
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
            amrex::ParticleReal const z = refpart.z;
            amrex::ParticleReal const t = refpart.t;
            amrex::ParticleReal const px = refpart.px;
            amrex::ParticleReal const pz = refpart.pz;
            amrex::ParticleReal const pt = refpart.pt;

            if (m_modify_ref_part && m_model == Model::nonlinear) {
                // calculate expensive terms
                amrex::ParticleReal const loc = m_location == Location::entry ? 1_prt : -1_prt;
                auto const [sin_psi, cos_psi] = amrex::Math::sincos(m_psi);
                amrex::ParticleReal const tan_psi = sin_psi/cos_psi;
                amrex::ParticleReal const sec_psi = 1_prt/cos_psi;
                amrex::ParticleReal const c2_ref = powi<2>(m_g/m_rc) * m_K0 * sin_psi * powi<3>(sec_psi)/2_prt;
                amrex::ParticleReal const c3_ref = loc * powi<2>(m_g)/m_rc * m_K0 * powi<2>(sec_psi);
                amrex::ParticleReal const c5_ref = tan_psi / (2_prt * m_rc);
                amrex::ParticleReal const c13_ref = powi<2>(sin_psi) / (2_prt*powi<3>(cos_psi)) * powi<2>(m_g)/(m_rc*m_R) * m_K4;

                amrex::ParticleReal const pnorm = std::sqrt(pt*pt-1_prt);
                amrex::ParticleReal const dx_ref = -c3_ref;
                amrex::ParticleReal const dpx_ref = c13_ref - c2_ref - c3_ref*c5_ref;
                amrex::ParticleReal const dpz_ref = pnorm*(-1_prt + std::sqrt(1_prt - powi<2>(dpx_ref/pnorm)));
                amrex::ParticleReal const dt_ref = -c3_ref * dpx_ref / refpart.beta();

                // advance position and momentum
                refpart.x = x + (pz/pnorm)*dx_ref;
                refpart.z = z - (px/pnorm)*dx_ref;
                refpart.t = t + dt_ref;
                refpart.px = px + (pz/pnorm)*dpx_ref + (px/pnorm)*dpz_ref;
                refpart.pz = pz - (px/pnorm)*dpx_ref + (pz/pnorm)*dpz_ref;
                // refpart.pt = pt;  // Unchanged
            }

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

            // additional constants required
            amrex::ParticleReal const loc = m_location == Location::entry ? 1_prt : -1_prt;
            auto const [sin_psi, cos_psi] = amrex::Math::sincos(m_psi);

            // reference particle relativistic factor and gyromagnetic anomaly
            amrex::ParticleReal const gamma_minus_1 = refpart.gamma() - 1_prt;
            amrex::ParticleReal const gyro_anomaly = refpart.gyromagnetic_anomaly;

            // The required momentum components in the coordinate frame defined by the pole-face:
            T_Real pz = sqrt(1_prt - 2_prt*pt/m_beta + pt*pt - px*px - py*py);
            T_Real px_face = px * cos_psi - pz * sin_psi;

            // The axis-angle vector in the coordinate frame defined by the pole-face.
            // This expression is taken from D. Abell, PRAB 18, 024001 (2015), eqs (41-42):
            T_Real lambdaz_face = -loc / m_rc * (1_prt + gyro_anomaly) * y;
            T_Real lambdax_face =  loc / m_rc * gyro_anomaly * gamma_minus_1 * px_face;

            // the axis-angle vector in the local coordinate frame defined by the reference particle
            T_Real const lambdax = lambdax_face * cos_psi + lambdaz_face * sin_psi;
            T_Real const lambday = 0_prt;
            T_Real const lambdaz = lambdaz_face * cos_psi - lambdax_face * sin_psi;

            // push the spin vector using the generator just determined
            rotate_spin(lambdax, lambday, lambdaz, sx, sy, sz);

            // phase space push
            (*this)(x, y, t, px, py, pt, idcpu, refpart);
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

            // initialize linear map matrix elements
            Map6x6 R = Map6x6::Identity();

            // edge focusing matrix elements (zero gap)
            amrex::ParticleReal const R21 = std::tan(m_psi) / m_rc;

            // first-order effect of nonzero gap
            auto const [sin_psi, cos_psi] = amrex::Math::sincos(m_psi);
            amrex::ParticleReal const vf = (1.0_prt + powi<2>(sin_psi))
                / powi<3>(cos_psi)
                * m_g * m_K2 / powi<2>(m_rc);

            amrex::ParticleReal const R43 = vf - R21;

            // set non-identity matrix elements
            R(2,1) = R21;
            R(4,3) = R43;

            // apply the transverse rotation (roll) alignment error
            return rotate_aligned_map(R);
        }

        amrex::ParticleReal m_psi; //! pole face angle in rad
        amrex::ParticleReal m_rc; //! bend radius in m
        amrex::ParticleReal m_g; //! gap parameter in m
        amrex::ParticleReal m_R; //! scale length in m
        amrex::ParticleReal m_K0; //! fringe field integral
        amrex::ParticleReal m_K1; //! fringe field integral
        amrex::ParticleReal m_K2; //! fringe field integral
        amrex::ParticleReal m_K3; //! fringe field integral
        amrex::ParticleReal m_K4; //! fringe field integral
        amrex::ParticleReal m_K5; //! fringe field integral
        amrex::ParticleReal m_K6; //! fringe field integral
        Model m_model; //! fringe field model: linear or nonlinear
        Location m_location; //! fringe field location: entry, or exit
        bool m_modify_ref_part; //! apply fringe field to reference particle (true/false)

      private:
        // constants that are independent of the individually tracked particle,
        // see: compute_constants() to refresh
        amrex::ParticleReal m_R21, m_R43, m_beta;
        amrex::ParticleReal m_c1, m_c5, m_c13, m_c14, m_c15;
        amrex::ParticleReal m_c2_times_1plusdelta, m_c3_times_1plusdelta, m_c4_times_1plusdelta, m_c6_times_1plusdelta;
        amrex::ParticleReal m_c7_times_1plusdelta, m_c8_times_1plusdelta, m_c9_times_1plusdelta, m_c10_times_1plusdelta;
        amrex::ParticleReal m_c11_times_1plusdelta, m_c12_times_1plusdelta;
    };

    /** DipEdge: GPU state plus common host-only metadata */
    using DipEdge = Element<DipEdgeGPU>;

} // namespace impactx

IMPACTX_PUSH_EXTERN_TEMPLATE(impactx::elements::DipEdgeGPU)

#endif // IMPACTX_DIPEDGE_H
