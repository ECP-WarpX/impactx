/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_PROGRAMMABLE_H
#define IMPACTX_ELEMENTS_PROGRAMMABLE_H

#include "particles/ImpactXParticleContainer.H"
#include "Element.H"
#include "particles/CovarianceMatrix.H"
#include "mixin/hostonly.H"
#include "mixin/lineartransport.H"
#include "mixin/thick.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>
// TODO: #include <AMReX_SIMD.H>
//   SIMD will require better complex math support in vir-simd / C++26
//   https://github.com/mattkretz/vir-simd/issues/42

#include <functional>
#include <stdexcept>

namespace impactx::elements
{
    /** An element whose push is supplied from Python
     *
     * Host-only: the callbacks below hold Python objects, which no kernel may touch. This
     * element is never turned into an execution snapshot -- its push is the callback
     * itself -- and the marker says so rather than leaving it to that fact.
     */
    struct ProgrammablePhysics
    : public mixin::LinearTransport<ProgrammablePhysics, false>,
      public mixin::HostOnly
      // TODO: public amrex::simd::Vectorized<amrex::simd::native_simd_size_particlereal>
    {
        static constexpr auto type = "Programmable";
        using PType = ImpactXParticleContainer::ParticleType;

        static constexpr amrex::ParticleReal DEFAULT_ds = 0;
        static constexpr int DEFAULT_nslice = mixin::Thick::DEFAULT_nslice;

        /** This element can be programmed
         *
         * @param ds Segment length in m
         * @param nslice number of slices used for the application of space charge
         */
        ProgrammablePhysics (
            amrex::ParticleReal ds = DEFAULT_ds,
            int nslice = DEFAULT_nslice
        )
          : m_ds(ds), m_nslice(nslice)
        {}

        /** Reverse the element (not supported for Programmable) */
        void reverse ()
        {
            throw std::runtime_error("Programmable: reverse() is not supported.");
        }

        /** Push all particles relative to the reference particle
         *
         * @param[in,out] pc particle container to push
         * @param[in] step global step for diagnostics
         * @param[in] period for periodic lattices, this is the current period (turn or cycle)
         */
        void operator() (
            ImpactXParticleContainer & pc,
            int step,
            int period
        ) const;

        /** Push all particles relative to the reference particle */
        void operator() (
            ImpactXParticleContainer::iterator & pti,
            RefPart & ref_part,
            [[maybe_unused]] bool spin
        ) const;

        /** This pushes the reference particle.
         *
         * @param[in,out] refpart reference particle
         */
        void operator() (
            RefPart & refpart
        ) const;

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
            // TODO: can/should we fuse the two elements?
            throw std::runtime_error("The Programmable element cannot be used for linear transport maps. Use the LinearMap element instead.");

            return Map6x6::Identity();
        }

        /** Number of slices used for the application of space charge
         *
         * @return positive integer
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        int nslice () const
        {
            return m_nslice;
        }

        /** Return the segment length
         *
         * @return value in meters
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::ParticleReal ds () const
        {
            return m_ds;
        }

        /** Close and deallocate all data and handles.
         */
        void
        finalize ();

        amrex::ParticleReal m_ds = 0.0;  //! segment length in m
        int m_nslice = DEFAULT_nslice;  //! number of slices used for the application of space charge

        /** Allow threading via OpenMP for the particle iterator loop
         *
         * This will only affect threading if the OMP backend is active.
         *
         * The default value is false to do the safe thing by default. Users
         * must opt-in their guarantee that their code is thread-safe.
         */
        bool m_threadsafe = false;

        std::function<void(ImpactXParticleContainer *, int, int)> m_push;  //! hook for push of whole container (pc, step, period)
        std::function<void(ImpactXParticleContainer::iterator *, RefPart &)> m_beam_particles;  //! hook for beam particles (pti, ref)
        std::function<void(RefPart &)> m_ref_particle;  //! hook for reference particle
        std::function<void()> m_finalize;  //! hook for finalize cleanup
    };

    /** Programmable: physics plus common host-only metadata */
    using Programmable = Element<ProgrammablePhysics>;

} // namespace impactx

#endif // IMPACTX_ELEMENTS_PROGRAMMABLE_H
