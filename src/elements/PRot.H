/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_PROT_H
#define IMPACTX_PROT_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/gpuelement.H"
#include "mixin/alignment.H"
#include "mixin/beamoptic.H"
#include "mixin/thin.H"
#include "mixin/lineartransport.H"
#include "mixin/spintransport.H"
#include "mixin/nofinalize.H"

#include <ablastr/constant.H>

#include <AMReX_Extension.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>
#include <AMReX_SIMD.H>

#include <cmath>
#include <stdexcept>

namespace impactx::elements
{
    struct PRotGPU
    : public mixin::BeamOptic<PRotGPU>,
      public mixin::LinearTransport<PRotGPU>,
      public mixin::Thin,
      public mixin::SpinTransport,
      public mixin::NoFinalize,
      public amrex::simd::Vectorized<amrex::simd::native_simd_size_particlereal>
    {
        static constexpr auto type = "PRot";
        using PType = ImpactXParticleContainer::ParticleType;

        /** An exact pole face rotation in the x-z plane, from a frame
         *  in which the reference orbit has angle phi_in with the z-axis,
         *  to a frame in which the reference orbit has angle phi_out
         *  with the z-axis.  The net rotation angle is phi_out - phi_in.
         *
         * @param phi_in Initial angle of reference trajectory w/r/t z (degrees)
         * @param phi_out Final angle of reference trajectory w/r/t/ z (degrees)
         */
        PRotGPU (
            amrex::ParticleReal phi_in,
            amrex::ParticleReal phi_out
        )
        : m_phi_in(phi_in * mixin::Alignment::degree2rad),
          m_phi_out(phi_out * mixin::Alignment::degree2rad)
        {
        }

        /** Reverse the element (swap phi_in and phi_out) */
        void reverse () { std::swap(m_phi_in, m_phi_out); }

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

            // access reference particle values to find 1/beta
            m_inv_beta = 1.0_prt / refpart.beta();

            // store rotation angle and initial, final values of pz
            amrex::ParticleReal const theta = m_phi_out - m_phi_in;
            auto const [sin_theta, cos_theta] = amrex::Math::sincos(theta);
            m_sin_theta = sin_theta;
            m_cos_theta = cos_theta;
            auto const [sin_phi_in, cos_phi_in] = amrex::Math::sincos(m_phi_in);
            m_sin_phi_in = sin_phi_in;
            m_cos_phi_in = cos_phi_in;
        }

        /** This is a prot functor, so that a variable of this type can be used like a
         *  prot function.
         *
         * The @see compute_constants method must be called before pushing particles through this operator.
         *
         * @param x particle position in x
         * @param y particle position in y
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
            T_Real & AMREX_RESTRICT x,
            T_Real & AMREX_RESTRICT y,
            T_Real & AMREX_RESTRICT t,
            T_Real & AMREX_RESTRICT px,
            T_Real const & AMREX_RESTRICT py,
            T_Real const & AMREX_RESTRICT pt,
            [[maybe_unused]] T_IdCpu const & AMREX_RESTRICT idcpu,
            [[maybe_unused]] RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using namespace std; // for cmath(float)
            using amrex::Math::powi;

            // initialize output values
            T_Real xout = x;
            T_Real yout = y;
            T_Real tout = t;
            T_Real pxout = px;
            // T_Real pyout = py;
            // T_Real ptout = pt;

            T_Real const pz = sqrt(
                1.0_prt -
                2.0_prt * pt * m_inv_beta +
                powi<2>(pt) -
                powi<2>(py) -
                powi<2>(px + m_sin_phi_in)
            );
            T_Real const pzf = pz * m_cos_theta - (px + m_sin_phi_in) * m_sin_theta;

            // advance position and momentum
            xout = x * pz / pzf;
            pxout = px * m_cos_theta + (pz - m_cos_phi_in) * m_sin_theta;

            yout = y + py * x * m_sin_theta / pzf;
            // pyout = py;

            tout = t - (pt - m_inv_beta) * x * m_sin_theta / pzf;
            // ptout = pt;

            // assign updated values
            x = xout;
            y = yout;
            t = tout;
            px = pxout;
            // py = pyout;
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
         * @param idcpu particle global index (unused)
         * @param refpart reference particle (unused)
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
            [[maybe_unused]] T_IdCpu const & AMREX_RESTRICT idcpu,
            [[maybe_unused]] RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // A simple rotation about the y-axis (in the x-z plane).
            T_Real const lambdax = 0_prt;
            T_Real const lambday = m_phi_out - m_phi_in;
            T_Real const lambdaz = 0_prt;

            // push the spin vector using the generator just determined
            rotate_spin(lambdax, lambday, lambdaz, sx, sy, sz);

            // phase space push
            (*this)(x, y, t, px, py, pt, idcpu, refpart);
        }

        /** This pushes the covariance matrix. */
        using LinearTransport::operator();

        /** This function returns the linear transport map.
         *
         * Jacobian at the reference trajectory of the exact pole-face
         * rotation implemented by @see operator(). With net rotation
         * angle theta = phi_out - phi_in, the nonzero entries are
         *   R(1,1) = cos(phi_in) / cos(phi_out),
         *   R(2,2) = cos(phi_out) / cos(phi_in),
         *   R(2,6) = -sin(theta) / (beta * cos(phi_in)),
         *   R(5,1) = +sin(theta) / (beta * cos(phi_out)),
         * and all other diagonal entries are 1. The derivation uses the
         * trig identity
         *   cos(phi_in) * cos(phi_in - phi_out)
         *     + sin(phi_in) * sin(phi_in - phi_out) = cos(phi_out),
         * which simplifies the x-denominator pz|_0 -> cos(phi_out).
         * Symplecticity of the reduced (x, px, t, pt) block was
         * verified with SymPy (R^T J R = J).
         *
         * The (y, py) sub-block is the identity: pole-face rotation in
         * the x-z plane does not mix the vertical plane at linear order.
         *
         * Reference: MaryLie 3.0 Users' Manual, Sect. 6.5 ("Plane Rotation");
         * A. J. Dragt, Lie Methods for Nonlinear Dynamics with Applications
         * to Accelerator Physics, v. Nov. 2020,
         * Appendix V ("PROT without and in the Presence of a Magnetic Field").
         *
         * Ideal behavior only: alignment and rotation errors are
         * ignored; the linearization is at (x, y, t) = (0, 0, 0) and
         * (px, py, pt) = (0, 0, 0).
         *
         * @param[in] refpart reference particle (used for beta)
         * @returns 6x6 transport matrix
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        Map6x6
        transport_map (RefPart const & AMREX_RESTRICT refpart) const
        {
            using namespace amrex::literals; // for _rt and _prt

            amrex::ParticleReal const beta = refpart.beta();

            // Precompute trigonometric quantities.
            amrex::ParticleReal const theta = m_phi_out - m_phi_in;
            auto const [sin_theta, cos_theta] = amrex::Math::sincos(theta);
            (void) cos_theta;  // not needed for the linear map
            auto const [sin_phi_in, cos_phi_in] = amrex::Math::sincos(m_phi_in);
            (void) sin_phi_in;
            auto const [sin_phi_out, cos_phi_out] = amrex::Math::sincos(m_phi_out);
            (void) sin_phi_out;

            Map6x6 R = Map6x6::Identity();
            R(1,1) = cos_phi_in / cos_phi_out;
            R(2,2) = cos_phi_out / cos_phi_in;
            R(2,6) = -sin_theta / (beta * cos_phi_in);
            R(5,1) = +sin_theta / (beta * cos_phi_out);

            return R;
        }

        amrex::ParticleReal m_phi_in; //! normalized (max) RF voltage drop.
        amrex::ParticleReal m_phi_out; //! RF wavenumber in 1/m.

      private:
        // constants that are independent of the individually tracked particle,
        // see: compute_constants() to refresh
        amrex::ParticleReal m_inv_beta, m_sin_theta, m_cos_theta, m_sin_phi_in, m_cos_phi_in;
    };

    /** PRot: GPU state plus common host-only metadata */
    using PRot = Element<PRotGPU>;

} // namespace impactx

IMPACTX_PUSH_EXTERN_TEMPLATE(impactx::elements::PRotGPU)

#endif // IMPACTX_PROT_H
