/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Andrew Myers
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_INITIALIZATION_INITDISTRIBUTION_H
#define IMPACTX_INITIALIZATION_INITDISTRIBUTION_H

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_Extension.H> // for AMREX_RESTRICT
#include <AMReX_REAL.H>

#include <type_traits>

namespace impactx::initialization
{
    /** Initialize a single particle's data using the given distribution
     *
     * Note: we usually would just write a C++ lambda below in ParallelFor. But, due to restrictions
     * in NVCC as of 11.5, we cannot write a lambda in a lambda as we also std::visit the element
     * types of our lattice elements list.
     *    error #3206-D: An extended __device__ lambda cannot be defined inside a generic lambda expression("operator()").
     * Thus, we fall back to writing a C++ functor here, instead of nesting two lambdas.
     *
     * Nvidia bug report: 3458976
     * Minimal demonstrator: https://cuda.godbolt.org/z/39e4q53Ye
     *
     * @tparam T_Distribution This can be a \see Gaussian, \see Waterbag, \see Kurth6D, \see Thermal etc.
     */
    template <typename T_Distribution>
    struct InitSingleParticleData
    {
        /** Constructor taking in pointers to particle data
         *
         * @param distribution the type of distribution function to call
         * @param part_x the array to the particle position (x)
         * @param part_y the array to the particle position (y)
         * @param part_t the array to the particle position (t)
         * @param part_px the array to the particle momentum (x)
         * @param part_py the array to the particle momentum (y)
         * @param part_pt the array to the particle momentum (t)
         */
        InitSingleParticleData (
            T_Distribution distribution,
            amrex::ParticleReal* AMREX_RESTRICT part_x,
            amrex::ParticleReal* AMREX_RESTRICT part_y,
            amrex::ParticleReal* AMREX_RESTRICT part_t,
            amrex::ParticleReal* AMREX_RESTRICT part_px,
            amrex::ParticleReal* AMREX_RESTRICT part_py,
            amrex::ParticleReal* AMREX_RESTRICT part_pt
        )
        : m_distribution(std::move(distribution)),
          m_part_x(part_x), m_part_y(part_y), m_part_t(part_t),
          m_part_px(part_px), m_part_py(part_py), m_part_pt(part_pt)
        {
        }

        InitSingleParticleData () = delete;
        InitSingleParticleData (InitSingleParticleData const &) = default;
        InitSingleParticleData (InitSingleParticleData &&) = default;
        ~InitSingleParticleData () = default;

        /** Initialize the data for a single particle
         *
         * @param i particle index in the current box
         * @param engine a random number engine (with associated state)
         */
        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        void
        operator() (
            amrex::Long i,
            amrex::RandomEngine const & engine
        ) const
        {
            m_distribution(m_part_x[i], m_part_y[i], m_part_t[i],
                           m_part_px[i], m_part_py[i], m_part_pt[i],
                           engine);
        }

    private:
        T_Distribution const m_distribution;
        amrex::ParticleReal* const AMREX_RESTRICT m_part_x;
        amrex::ParticleReal* const AMREX_RESTRICT m_part_y;
        amrex::ParticleReal* const AMREX_RESTRICT m_part_t;
        amrex::ParticleReal* const AMREX_RESTRICT m_part_px;
        amrex::ParticleReal* const AMREX_RESTRICT m_part_py;
        amrex::ParticleReal* const AMREX_RESTRICT m_part_pt;
    };
}  // namespace impactx::initialization

#endif  // IMPACTX_INITIALIZATION_INITDISTRIBUTION_H
