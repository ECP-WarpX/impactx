/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
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
#include "mixin/named.H"
#include "mixin/thick.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>

#include <functional>


namespace impactx
{
    struct Programmable
    : public elements::Named
    {
        static constexpr auto type = "Programmable";
        using PType = ImpactXParticleContainer::ParticleType;

        /** This element can be programmed
         *
         * @param name a user defined and not necessarily unique name of the element
         * @param ds Segment length in m
         * @param nslice number of slices used for the application of space charge
         */
        Programmable (
            std::optional<std::string> name,
            amrex::ParticleReal ds=0.0,
            int nslice=1
        )
          : Named(name), m_ds(ds), m_nslice(nslice)
        {}

        /** Push all particles relative to the reference particle
         *
         * @param[in,out] pc particle container to push
         * @param[in] step global step for diagnostics
         */
        void operator() (
            ImpactXParticleContainer & pc,
            int step
        ) const;

        /** Push all particles relative to the reference particle */
        void operator() (
            ImpactXParticleContainer::iterator & pti,
            RefPart & ref_part
        ) const;

        /** This pushes the reference particle.
         *
         * @param[in,out] refpart reference particle
         */
        void operator() (
            RefPart & refpart
        ) const;

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

        amrex::ParticleReal m_ds = 0.0; //! segment length in m
        int m_nslice = 1; //! number of slices used for the application of space charge

        /** Allow threading via OpenMP for the particle iterator loop
         *
         * This will only affect threading if the OMP backend is active.
         *
         * The default value is false to do the safe thing by default. Users
         * must opt-in their guarantee that their code is thread-safe.
         */
        bool m_threadsafe = false;

        std::function<void(ImpactXParticleContainer *, int)> m_push; //! hook for push of whole container
        std::function<void(ImpactXParticleContainer::iterator *, RefPart &)> m_beam_particles; //! hook for beam particles
        std::function<void(RefPart &)> m_ref_particle; //! hook for reference particle
        std::function<void()> m_finalize; //! hook for finalize cleanup
    };

} // namespace impactx

#endif // IMPACTX_ELEMENTS_PROGRAMMABLE_H
