/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl, Kyrre Ness Sjobak
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_QUAD_H
#define IMPACTX_QUAD_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/gpuelement.H"
#include "mixin/alignment.H"
#include "mixin/pipeaperture.H"
#include "mixin/beamoptic.H"
#include "mixin/thick.H"
#include "mixin/nofinalize.H"
#include "mixin/lineartransport.H"
#include "mixin/spintransport.H"

#include <AMReX_Extension.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>
#include <AMReX_SIMD.H>

#include <cmath>

namespace impactx::elements
{
    struct QuadGPU
    : public mixin::BeamOptic<QuadGPU>,
      public mixin::LinearTransport<QuadGPU>,
      public mixin::Thick,
      public mixin::Alignment,
      public mixin::PipeAperture,
      public mixin::SpinTransport,
      public mixin::NoFinalize,
      public amrex::simd::Vectorized<amrex::simd::native_simd_size_particlereal>
    {
        static constexpr auto type = "Quad";
        using PType = ImpactXParticleContainer::ParticleType;

        // The light no-spin push in double precision does not benefit from SIMD,
        // see https://github.com/BLAST-ImpactX/impactx/pull/1002
        // This overwrites in mixin::BeamOptic the amrex::simd::Vectorized flag.
#ifndef AMREX_SINGLE_PRECISION_PARTICLES
        static constexpr bool simd_ps = false;
#endif

        /** A Quadrupole magnet
         *
         * @param ds Segment length in m.
         * @param k  Quadrupole strength in m^(-2) (MADX convention)
         *           = (gradient in T/m) / (rigidity in T-m)
         *           k > 0 horizontal focusing
         *           k < 0 horizontal defocusing
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param aperture_x horizontal half-aperture in m
         * @param aperture_y vertical half-aperture in m
         * @param nslice number of slices used for the application of space charge
         */
        QuadGPU (
            amrex::ParticleReal ds,
            amrex::ParticleReal k,
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
          m_k(k)
        {
        }

        /** Reverse the element (negate ds) */
        void reverse () { Thick::reverse(); }

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

            // length of the current slice
            m_slice_ds = m_ds / nslice();

            amrex::ParticleReal const pt_ref = refpart.pt;
            // find beta*gamma^2
            m_betgam2 = powi<2>(pt_ref) - 1.0_prt;
            m_slice_bg = m_slice_ds / m_betgam2;

            // gyromagnetic constant (used during spin push)
            amrex::ParticleReal const G = refpart.gyromagnetic_anomaly;
            m_gyro_const = 1_prt + G * refpart.gamma();

            // Resolve the focusing/defocusing/drift branch here, once per
            // element, so the per-particle operator() and spin push are
            // branch-free. omega = sqrt(|k|) is independent of the tracked
            // particle, so all trigonometric quantities can be cached.
            // The transverse 6x6 transfer-map blocks are:
            //   xout = m_R11*x + m_R12*px ;  pxout = m_R21*x + m_R22*px
            //   yout = m_R33*y + m_R34*py ;  pyout = m_R43*y + m_R44*py
            // and the spin generator is:
            //   lambdax = m_spin_py*py + m_spin_y*y
            //   lambday = m_spin_px*px + m_spin_x*x
            amrex::ParticleReal const omega = std::sqrt(std::abs(m_k));
            amrex::ParticleReal const sin_omega_ds = std::sin(omega * m_slice_ds);
            amrex::ParticleReal const cos_omega_ds = std::cos(omega * m_slice_ds);
            amrex::ParticleReal const sinh_omega_ds = std::sinh(omega * m_slice_ds);
            amrex::ParticleReal const cosh_omega_ds = std::cosh(omega * m_slice_ds);

            if (m_k > 0.0_prt)
            {
                // horizontally focusing quad (x focusing, y defocusing)
                m_R11 = cos_omega_ds;
                m_R12 = sin_omega_ds / omega;
                m_R21 = -omega * sin_omega_ds;
                m_R22 = cos_omega_ds;
                m_R33 = cosh_omega_ds;
                m_R34 = sinh_omega_ds / omega;
                m_R43 = omega * sinh_omega_ds;
                m_R44 = cosh_omega_ds;

                m_spin_py = -m_gyro_const * (cosh_omega_ds - 1_prt);
                m_spin_y  = -m_gyro_const * omega * sinh_omega_ds;
                m_spin_px = -m_gyro_const * (1_prt - cos_omega_ds);
                m_spin_x  = -m_gyro_const * omega * sin_omega_ds;
            }
            else if (m_k < 0.0_prt)
            {
                // horizontally defocusing quad (x defocusing, y focusing)
                m_R11 = cosh_omega_ds;
                m_R12 = sinh_omega_ds / omega;
                m_R21 = omega * sinh_omega_ds;
                m_R22 = cosh_omega_ds;
                m_R33 = cos_omega_ds;
                m_R34 = sin_omega_ds / omega;
                m_R43 = -omega * sin_omega_ds;
                m_R44 = cos_omega_ds;

                m_spin_py = m_gyro_const * (1_prt - cos_omega_ds);
                m_spin_y  = m_gyro_const * omega * sin_omega_ds;
                m_spin_px = m_gyro_const * (cosh_omega_ds - 1_prt);
                m_spin_x  = m_gyro_const * omega * sinh_omega_ds;
            }
            else
            {
                // zero strength = drift
                m_R11 = 1.0_prt;
                m_R12 = m_slice_ds;
                m_R21 = 0.0_prt;
                m_R22 = 1.0_prt;
                m_R33 = 1.0_prt;
                m_R34 = m_slice_ds;
                m_R43 = 0.0_prt;
                m_R44 = 1.0_prt;

                m_spin_py = 0.0_prt;
                m_spin_y  = 0.0_prt;
                m_spin_px = 0.0_prt;
                m_spin_x  = 0.0_prt;
            }
        }

        /** This is a quad functor, so that a variable of this type can be used like a quad function.
         *
         * The @see compute_constants method must be called before pushing particles through this operator.
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
            T_Real & AMREX_RESTRICT py,
            T_Real const & AMREX_RESTRICT pt,
            [[maybe_unused]] T_IdCpu const & AMREX_RESTRICT idcpu,
            [[maybe_unused]] RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // copies of original input values
            T_Real const xin = x;
            T_Real const yin = y;

            // advance position and momentum
            x  = m_R11 * xin + m_R12 * px;
            px = m_R21 * xin + m_R22 * px;

            y  = m_R33 * yin + m_R34 * py;
            py = m_R43 * yin + m_R44 * py;

            // R56 = m_slice_bg ;  R55 = R66 = 1
            t = t + m_slice_bg * pt;
            // pt unchanged
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

            // assign intermediate parameter
            amrex::ParticleReal const step = slice_ds / std::sqrt(powi<2>(pt) - 1.0_prt);

            // advance position and momentum (straight element)
            refpart.x = x + step*px;
            refpart.y = y + step*py;
            refpart.z = z + step*pz;
            refpart.t = t - step*pt;

            // advance integrated path length
            refpart.s = s + slice_ds;
        }

        /** This function returns the linear transport map.
         *
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
            amrex::ParticleReal const betgam2 = powi<2>(pt_ref) - 1.0_prt;

            // compute phase advance per unit length in s (in rad/m)
            amrex::ParticleReal const omega = std::sqrt(std::abs(m_k));

            // initialize linear map matrix elements
            Map6x6 R = Map6x6::Identity();

            if (m_k > 0.0) {
                R(1,1) = std::cos(omega*slice_ds);
                R(1,2) = std::sin(omega*slice_ds)/omega;
                R(2,1) = -omega*std::sin(omega*slice_ds);
                R(2,2) = std::cos(omega*slice_ds);
                R(3,3) = std::cosh(omega*slice_ds);
                R(3,4) = std::sinh(omega*slice_ds)/omega;
                R(4,3) = omega*std::sinh(omega*slice_ds);
                R(4,4) = std::cosh(omega*slice_ds);
                R(5,6) = slice_ds/betgam2;
            } else if (m_k < 0.0) {
                R(1,1) = std::cosh(omega*slice_ds);
                R(1,2) = std::sinh(omega*slice_ds)/omega;
                R(2,1) = omega*std::sinh(omega*slice_ds);
                R(2,2) = std::cosh(omega*slice_ds);
                R(3,3) = std::cos(omega*slice_ds);
                R(3,4) = std::sin(omega*slice_ds)/omega;
                R(4,3) = -omega*std::sin(omega*slice_ds);
                R(4,4) = std::cos(omega*slice_ds);
                R(5,6) = slice_ds/betgam2;
            } else {
                R(1,2) = slice_ds;
                R(3,4) = slice_ds;
                R(5,6) = slice_ds / betgam2;
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
         * @param refpart reference particle (unused)
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

            // the axis-angle vector
            // The focusing/defocusing/drift selection is fully resolved in
            // compute_constants(), see the m_spin_* coefficients there.
            T_Real const lambdax = m_spin_py * py + m_spin_y * y;
            T_Real const lambday = m_spin_px * px + m_spin_x * x;
            T_Real const lambdaz = 0_prt;

            // push the spin vector using the generator just determined
            rotate_spin(lambdax, lambday, lambdaz, sx, sy, sz);

            // phase space push
            (*this)(x, y, t, px, py, pt, idcpu, refpart);
        }

        /** This pushes the covariance matrix. */
        using LinearTransport::operator();

        amrex::ParticleReal m_k; //! quadrupole strength in 1/m

    private:
        // constants that are independent of the individually tracked particle,
        // see: compute_constants() to refresh
        amrex::ParticleReal m_slice_ds;       //! m_ds / nslice();
        amrex::ParticleReal m_betgam2;        //! beta*gamma^2
        amrex::ParticleReal m_slice_bg;       //! R56 = m_slice_ds / m_betgam2
        amrex::ParticleReal m_gyro_const;     //! (1 + G*gamma) for spin calculation
        // pre-resolved 6x6 transfer-map elements (transverse blocks), see compute_constants()
        amrex::ParticleReal m_R11, m_R12, m_R21, m_R22; //! x--px block
        amrex::ParticleReal m_R33, m_R34, m_R43, m_R44; //! y--py block
        // branch-free spin generator coefficients, see compute_constants()
        amrex::ParticleReal m_spin_px;        //! lambday += m_spin_px*px
        amrex::ParticleReal m_spin_x;         //! lambday += m_spin_x*x
        amrex::ParticleReal m_spin_py;        //! lambdax += m_spin_py*py
        amrex::ParticleReal m_spin_y;         //! lambdax += m_spin_y*y

    };

    /** Quad: GPU state plus common host-only metadata */
    using Quad = Element<QuadGPU>;

} // namespace impactx

IMPACTX_PUSH_EXTERN_TEMPLATE(impactx::elements::QuadGPU)

#endif // IMPACTX_QUAD_H
