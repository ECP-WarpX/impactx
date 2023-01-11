/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_MIXIN_THIN_H
#define IMPACTX_ELEMENTS_MIXIN_THIN_H

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>


namespace impactx::elements
{
    /** This is a helper class for lattice elements with zero length
     */
    struct Thin
    {
        /** This pushes the reference particle.
         *
         * @param[in,out] refpart reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() ([[maybe_unused]] RefPart & AMREX_RESTRICT refpart) const
        {
            // nothing to do: this is a zero-length element
        }

        /** Number of slices used for the application of space charge
         *
         * @return one, because this is a zero-length element
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        int nslice () const
        {
            return 1;
        }

        /** Return the segment length
         *
         * @return zero, because this is a zero-length element
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::ParticleReal ds () const
        {
            using namespace amrex::literals;
            return 0.0_prt;
        }
    };

} // namespace impactx::elements

#endif // IMPACTX_ELEMENTS_MIXIN_THIN_H
