/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_MIXIN_THICK_H
#define IMPACTX_ELEMENTS_MIXIN_THICK_H

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>


namespace impactx::elements
{
    /** This is a helper class for lattice elements with finite length
     */
    struct Thick
    {
        /** A finite-length element
         *
         * @param ds Segment length in m
         * @param nslice number of slices used for the application of space charge
         */
        Thick (
            amrex::ParticleReal ds,
            int nslice
        )
        : m_ds(ds), m_nslice(nslice)
        {
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

        amrex::ParticleReal m_ds; //! segment length in m
        int m_nslice; //! number of slices used for the application of space charge
    };

} // namespace impactx::elements

#endif // IMPACTX_ELEMENTS_MIXIN_THICK_H
