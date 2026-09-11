/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENT_LINEAR_MAP_H
#define IMPACTX_ELEMENT_LINEAR_MAP_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/gpuelement.H"
#include "mixin/alignment.H"
#include "mixin/beamoptic.H"
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
#include <type_traits>

namespace impactx::elements
{
    struct LinearMapGPU
    : public mixin::BeamOptic<LinearMapGPU>,
      public mixin::LinearTransport<LinearMapGPU>,
      public mixin::Alignment,
      public mixin::PipeAperture,
      public mixin::SpinTransport,
      public mixin::NoFinalize,
      public amrex::simd::Vectorized<amrex::simd::native_simd_size_particlereal>
    {
        static constexpr auto type = "LinearMap";
        using PType = ImpactXParticleContainer::ParticleType;

        static constexpr amrex::ParticleReal DEFAULT_ds = 0;

        /** A thin element that applies a user-provided linear transport map R
         *  to the 6-vector of phase space coordinates (x,px,y,py,t,pt).
         *  Thus  x_final = R(1,1)*x + R(1,2)*px + R(1,3)*y + ...,
         *       px_final = R(2,1)*x + R(2,2)*px + R(2,3)*y + ..., etc.
         *
         * @param R user-provided transport map
         * @param ds Segment length in m
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param aperture_x horizontal half-aperture in m
         * @param aperture_y vertical half-aperture in m
         */
        LinearMapGPU (
            Map6x6 const & R,
            amrex::ParticleReal ds = DEFAULT_ds,
            amrex::ParticleReal dx = DEFAULT_dx,
            amrex::ParticleReal dy = DEFAULT_dy,
            amrex::ParticleReal rotation_degree = DEFAULT_rotation_degree,
            amrex::ParticleReal aperture_x = DEFAULT_aperture_x,
            amrex::ParticleReal aperture_y = DEFAULT_aperture_y
        )
        : Alignment(dx, dy, rotation_degree),
          PipeAperture(aperture_x, aperture_y)
        {
            m_transport_map = R;
            m_ds = ds;
        }

        /** Check if the transport map is symplectic
         *
         * A matrix R is symplectic if R^T Omega R = Omega, where Omega is the
         * standard 6x6 symplectic form.
         *
         * @see https://en.wikipedia.org/wiki/Symplectic_matrix
         * @return true if the map satisfies the symplectic condition
         */
        bool symplectic () const
        {
            using namespace amrex::literals;

            Map6x6 const J = impactx::symplectic_form();
            Map6x6 const RT = m_transport_map.transpose();
            Map6x6 const check = RT * J * m_transport_map;
            amrex::ParticleReal max_err = 0_prt;
            for (int i = 1; i <= 6; ++i) {
                for (int j = 1; j <= 6; ++j) {
                    amrex::ParticleReal const err = std::abs(check(i,j) - J(i,j));
                    max_err = std::max(max_err, err);
                }
            }

            amrex::ParticleReal constexpr tol = std::is_same_v<amrex::ParticleReal, float>
                ? 1.0e-4_prt : 1.0e-10_prt;
            return max_err <= tol;
        }

        /** Reverse the element (compute symplectic inverse of transport map)
         *
         * For a symplectic matrix R, the inverse is
         *   R^{-1} = J^T R^T J = -J R^T J
         * where J is the standard 6x6 symplectic form.
         *
         * @see https://en.wikipedia.org/wiki/Symplectic_matrix
         * Throws if the transport map is not symplectic.
         */
        void reverse ()
        {
            if (!symplectic()) {
                throw std::runtime_error(
                    "LinearMap::reverse: transport map is not symplectic.");
            }

            Map6x6 const Omega = impactx::symplectic_form();

            // R_inv = Omega^T * R^T * Omega = -Omega * R^T * Omega
            Map6x6 const RT = m_transport_map.transpose();
            Map6x6 const R_inv = (-Omega) * RT * Omega;
            m_transport_map = R_inv;
            m_ds = -m_ds;
        }

        /** Push all particles */
        using BeamOptic::operator();

        /** This is a LinearMap functor, so that a variable of this type can be used like a
         *  LinearMap function.
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
            T_Real & AMREX_RESTRICT py,
            T_Real & AMREX_RESTRICT pt,
            [[maybe_unused]] T_IdCpu const & AMREX_RESTRICT idcpu,
            [[maybe_unused]] RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // input and output phase space vectors
            amrex::SmallVector<T_Real, 6, 1> const vectorin{
                x, px, y, py, t, pt
            };

            // TODO amrex::SmallMatrix::operator* overloads on T for SIMD*scalar / scalar*SIMD
            amrex::SmallVector<T_Real, 6, 1> const vectorout = m_transport_map * vectorin;

            // assign updated values
            x = vectorout(1);
            px = vectorout(2);
            y = vectorout(3);
            py = vectorout(4);
            t = vectorout(5);
            pt = vectorout(6);
        }

        /** This pushes the reference particle.
         *
         * @param[in,out] refpart reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (RefPart & AMREX_RESTRICT refpart) const
        {
            if (m_ds != 0)  // Drift forward or backward
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
                amrex::ParticleReal const step = slice_ds / std::sqrt(powi<2>(pt)-1.0_prt);

                // advance position and momentum (drift)
                refpart.x = x + step*px;
                refpart.y = y + step*py;
                refpart.z = z + step*pz;
                refpart.t = t - step*pt;

                // advance integrated path length
                refpart.s = s + slice_ds;
            }
            // else nothing to do for a zero-length element
        }

        /** Number of slices used for the application of space charge
         *
         * @return one, because we do not support slicing of this element
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        int nslice () const
        {
            return 1;
        }

        /** Return the segment length
         *
         * @return by default zero, but users can set a corresponding ds for bookkeeping
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::ParticleReal ds () const
        {
            return m_ds;
        }

        /** This operator pushes the particle spin.
         *
         * @param x particle position in x
         * @param y particle position in y
         * @param t particle position in t
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t
         * @param idcpu particle global index
         * @param refpart reference particle (unused)
         * @param sx particle spin x-component (unused)
         * @param sy particle spin y-component (unused)
         * @param sz particle spin z-component (unused)
         */
        template<typename T_Real=amrex::ParticleReal, typename T_IdCpu=uint64_t>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void spin_and_phasespace_push (
            T_Real & AMREX_RESTRICT x,
            T_Real & AMREX_RESTRICT y,
            T_Real & AMREX_RESTRICT t,
            T_Real & AMREX_RESTRICT px,
            T_Real & AMREX_RESTRICT py,
            T_Real & AMREX_RESTRICT pt,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT sx,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT sy,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT sz,
            T_IdCpu const & AMREX_RESTRICT idcpu,
            RefPart const & AMREX_RESTRICT refpart
        ) const
        {
             // spin is unaffected

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
            Map6x6 R = Map6x6::Identity();
            R = m_transport_map;

            // apply the transverse rotation (roll) alignment error
            return rotate_aligned_map(R);
        }

        Map6x6 m_transport_map;  // 6x6 transport map
        amrex::ParticleReal m_ds;  // finite ds allowed for bookkeeping, but we do not allow slicing
    };

    /** LinearMap: GPU state plus common host-only metadata */
    using LinearMap = Element<LinearMapGPU>;

} // namespace impactx

IMPACTX_PUSH_EXTERN_TEMPLATE(impactx::elements::LinearMapGPU)

#endif // IMPACTX_ELEMENT_LINEAR_MAP_H
