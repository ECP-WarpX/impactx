/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENT_SPIN_MAP_H
#define IMPACTX_ELEMENT_SPIN_MAP_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/gpuelement.H"
#include "mixin/alignment.H"
#include "mixin/beamoptic.H"
#include "mixin/lineartransport.H"
#include "mixin/spintransport.H"
#include "mixin/nofinalize.H"

#include <AMReX_Extension.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>
#include <AMReX_SIMD.H>
#include <AMReX_SmallMatrix.H>

#include <cmath>
#include <stdexcept>

/** This defines the 3x6 and 3x1 matrix types used here. */

using Map3x6 = amrex::SmallMatrix<amrex::ParticleReal, 3, 6, amrex::Order::F, 1>;
using Vector3 = amrex::SmallMatrix<amrex::ParticleReal, 3, 1, amrex::Order::F, 1>;

namespace impactx::elements
{
    struct SpinMapGPU
    : public mixin::BeamOptic<SpinMapGPU>,
      public mixin::LinearTransport<SpinMapGPU>,
      public mixin::Alignment,
      public mixin::SpinTransport,
      public mixin::NoFinalize,
      public amrex::simd::Vectorized<amrex::simd::native_simd_size_particlereal>
    {
        static constexpr auto type = "SpinMap";
        using PType = ImpactXParticleContainer::ParticleType;

        static constexpr amrex::ParticleReal DEFAULT_ds = 0;

        /** A thin element that applies a user-provided spin map M
         *  to the spin 3-vector (sx, sy, sz).  The map M is an
         *  orthogonal transformation expressed in Lie-algebraic
         *  (axis-angle) form, whose Lie generators may depend
         *  linearly on the phase space (orbit) variables (x,px,y,py,t,pt).
         *
         * @param v user-provided 3-vector describing spin rotation at the design point
         * @param A user-provided 3x6 spin-orbit coupling matrix
         * @param ds Segment length in m
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         */
        SpinMapGPU (
            Vector3 const & v,
            Map3x6 const & A,
            amrex::ParticleReal ds = DEFAULT_ds,
            amrex::ParticleReal dx = DEFAULT_dx,
            amrex::ParticleReal dy = DEFAULT_dy,
            amrex::ParticleReal rotation_degree = DEFAULT_rotation_degree
        )
        : Alignment(dx, dy, rotation_degree)
        {
            m_spin_rotation_vector = v;
            m_spin_orbit_coupling = A;
            m_ds = ds;
        }

        /** Reverse the element (reverse orbit-dependent and design spin rotations) */
        void reverse ()
        {
            throw std::runtime_error("SpinMap: reverse() is not yet supported.");

            /* Reversibility: This works but is hacky
            for (int i = 1; i <= 3; ++i) {
                m_spin_rotation_vector(i) = -m_spin_rotation_vector(i);
                for (int j = 1; j <= 6; ++j) {
                    m_spin_orbit_coupling(i,j) = -m_spin_orbit_coupling(i,j);
                }
            }
            m_reversed = !m_reversed;
            */
        }

        /** Push all particles */
        using BeamOptic::operator();

        /** This is a SpinMap functor, so that a variable of this type can be used like a
         *  SpinMap function.
         *
         * @param x particle position in x (unused)
         * @param y particle position in y (unused)
         * @param t particle position in t (unused)
         * @param px particle momentum in x (unused)
         * @param py particle momentum in y (unused)
         * @param pt particle momentum in t (unused)
         * @param idcpu particle global index (unused)
         * @param refpart reference particle (unused)
         */
        template<typename T_Real=amrex::ParticleReal, typename T_IdCpu=uint64_t>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
            [[maybe_unused]] T_Real const & AMREX_RESTRICT x,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT y,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT t,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT px,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT py,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT pt,
            [[maybe_unused]] T_IdCpu const & AMREX_RESTRICT idcpu,
            [[maybe_unused]] RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            // This element does nothing to the phase space variables.
        }

        /** This pushes the reference particle.
         *
         * @param[in,out] refpart reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() ([[maybe_unused]] RefPart & AMREX_RESTRICT refpart) const
        {
            // This element does nothing to the reference particle.
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
         * @param sx particle spin x-component
         * @param sy particle spin y-component
         * @param sz particle spin z-component
         * @param idcpu particle global index
         * @param refpart reference particle (unused)
         */
        template<typename T_Real=amrex::ParticleReal, typename T_IdCpu=uint64_t>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void spin_and_phasespace_push (
            T_Real const & AMREX_RESTRICT x,
            T_Real const & AMREX_RESTRICT y,
            T_Real const & AMREX_RESTRICT t,
            T_Real const & AMREX_RESTRICT px,
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

            // initialize the three components of the axis-angle vector
            T_Real lambdax = 0_prt;
            T_Real lambday = 0_prt;
            T_Real lambdaz = 0_prt;

            // input phase space 6-vector
            amrex::SmallVector<T_Real, 6, 1> const vectorin{
                x, px, y, py, t, pt
            };

            // determine the rotation vector for off-design particles
            amrex::SmallVector<T_Real, 3, 1> const vectorout = m_spin_orbit_coupling * vectorin;

            // assign updated values
            lambdax = vectorout(1);
            lambday = vectorout(2);
            lambdaz = vectorout(3);

            /* Reversibility: This works but is hacky
            T_Real const design_lambdax = m_spin_rotation_vector(1);
            T_Real const design_lambday = m_spin_rotation_vector(2);
            T_Real const design_lambdaz = m_spin_rotation_vector(3);

            if (m_reversed) {
                // The inverse map applies the inverse design rotation first,
                // then the inverse orbit-dependent rotation.
                rotate_spin(
                    design_lambdax,
                    design_lambday,
                    design_lambdaz,
                    sx, sy, sz
                );
                rotate_spin(lambdax,lambday,lambdaz,sx,sy,sz);
            } else {
                rotate_spin(lambdax,lambday,lambdaz,sx,sy,sz);
                rotate_spin(
                    design_lambdax,
                    design_lambday,
                    design_lambdaz,
                    sx, sy, sz
                );
            }
            */

            { /* non-reversible block, replace with the above block to test */
            // push the spin vector using the generator just determined
            rotate_spin(lambdax,lambday,lambdaz,sx,sy,sz);

            // specify the design rotation vector
            lambdax = m_spin_rotation_vector(1);
            lambday = m_spin_rotation_vector(2);
            lambdaz = m_spin_rotation_vector(3);

            // push the spin vector using the generator just determined
            rotate_spin(lambdax,lambday,lambdaz,sx,sy,sz);
            } /* end of non-reversible block, replace with the above block to test */

            // phase space push (the identity for this element)
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
            // apply the transverse rotation (roll) alignment error
            return rotate_aligned_map(R);
        }

        Vector3 m_spin_rotation_vector;  // 3x1 axis-angle vector for design spin rotation
        Map3x6 m_spin_orbit_coupling;  // 3x6 spin-orbit coupling matrix
        amrex::ParticleReal m_ds;  // finite ds allowed for bookkeeping, but we do not allow slicing
        /* Reversibility: This works but is hacky */
        // bool m_reversed = false;  // apply inverse spin rotations in reverse order after reverse()
    };

    /** SpinMap: GPU state plus common host-only metadata */
    using SpinMap = Element<SpinMapGPU>;

} // namespace impactx

IMPACTX_PUSH_EXTERN_TEMPLATE(impactx::elements::SpinMapGPU)

#endif // IMPACTX_ELEMENT_SPIN_MAP_H
