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

#include "mixin/thick.H"
#include "particles/ImpactXParticleContainer.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>

#include <functional>


namespace impactx
{
    struct Programmable
    {
        static constexpr auto name = "Programmable";
        using PType = ImpactXParticleContainer::ParticleType;

        /** This element writes the particle beam out to openPMD data.
         */
        Programmable () {}

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

    public:
        amrex::ParticleReal m_ds = 0.0; //! segment length in m
        int m_nslice = 1; //! number of slices used for the application of space charge

        std::function<void(ImpactXParticleContainer::iterator *, RefPart &)> m_beam_particles; //! hook for beam particles
        std::function<void(RefPart &)> m_ref_particle; //! hook for reference particle
    };

} // namespace impactx

#endif // IMPACTX_ELEMENTS_PROGRAMMABLE_H