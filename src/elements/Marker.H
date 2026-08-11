/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENT_MARKER_H
#define IMPACTX_ELEMENT_MARKER_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/gpuelement.H"
#include "mixin/thin.H"
#include "mixin/lineartransport.H"
#include "mixin/spintransport.H"
#include "mixin/nofinalize.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>
#include <AMReX_SIMD.H>

namespace impactx::elements
{
    struct MarkerGPU
    : public mixin::LinearTransport<MarkerGPU>,
      public mixin::Thin,
      public mixin::SpinTransport,
      public mixin::NoFinalize,
      public amrex::simd::Vectorized<amrex::simd::native_simd_size_particlereal>
    {
        static constexpr auto type = "Marker";
        using PType = ImpactXParticleContainer::ParticleType;

        /** This named element does nothing.
         *
         * The name is carried by the element envelope, so the GPU part is stateless.
         * It must stay default-constructible: that is what lets the envelope build the
         * GPU part from an empty argument prefix and route the name to the metadata.
         */
        MarkerGPU () = default;

        /** Reverse the element (no-op) */
        void reverse () { }

        /** Push all particles - nothing to do here */
        void operator() (
            ImpactXParticleContainer & /* pc */,
            int /* step */,
            int /* period */
        ) {
            // nothing to do
        }

        /** Compute and cache the constants for the push - nothing to do here */
        void compute_constants (RefPart const & /* refpart */)
        {
            // nothing to do
        }

        /** Push all particles - nothing to do here */
        void operator() (
            ImpactXParticleContainer::iterator & /* pti */,
            RefPart & AMREX_RESTRICT /* ref_part */
        ) {
            // nothing to do
        }

        /** Does nothing to a particle.
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
            // nothing to do
        }

        /** This pushes the reference particle. */
        using Thin::operator();

        /** This operator pushes the particle spin.
         *
         * A marker does nothing, so the spin is unaffected.
         *
         * @param x particle position in x (unused)
         * @param y particle position in y (unused)
         * @param t particle position in t (unused)
         * @param px particle momentum in x (unused)
         * @param py particle momentum in y (unused)
         * @param pt particle momentum in t (unused)
         * @param sx particle spin x-component (unused)
         * @param sy particle spin y-component (unused)
         * @param sz particle spin z-component (unused)
         * @param idcpu particle global index (unused)
         * @param refpart reference particle (unused)
         */
        template<typename T_Real=amrex::ParticleReal, typename T_IdCpu=uint64_t>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void spin_and_phasespace_push (
            [[maybe_unused]] T_Real const & AMREX_RESTRICT x,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT y,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT t,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT px,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT py,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT pt,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT sx,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT sy,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT sz,
            [[maybe_unused]] T_IdCpu const & AMREX_RESTRICT idcpu,
            [[maybe_unused]] RefPart const & AMREX_RESTRICT refpart
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
            // nothing to do
            return Map6x6::Identity();
        }
    };

    /** Marker: GPU state plus common host-only metadata */
    using Marker = Element<MarkerGPU>;

} // namespace impactx

#endif // IMPACTX_ELEMENT_MARKER_H
