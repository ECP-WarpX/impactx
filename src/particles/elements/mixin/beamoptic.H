/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_MIXIN_BEAMOPTIC_H
#define IMPACTX_ELEMENTS_MIXIN_BEAMOPTIC_H

#include "particles/ImpactXParticleContainer.H"
#include "particles/PushAll.H"

#include <AMReX_Extension.H> // for AMREX_RESTRICT
#include <AMReX_REAL.H>

#include <type_traits>


namespace impactx::elements
{
namespace detail
{
    /** Push a single particle through an element
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
     * @tparam T_Element This can be a \see Drift, \see Quad, \see Sbend, etc.
     */
    template <typename T_Element>
    struct PushSingleParticle
    {
        using PType = ImpactXParticleContainer::ParticleType;
        using ParticleTileType = amrex::ParticleTile<amrex::SoAParticle<RealSoA::nattribs, IntSoA::nattribs>,RealSoA::nattribs, IntSoA::nattribs>;

        /** Constructor taking in pointers to particle data
         *
         * @param element the beamline element to push through
         * @param part_x the array to the particle position (x)
         * @param part_y the array to the particle position (y)
         * @param part_t the array to the particle position (t)
         * @param part_px the array to the particle momentum (x)
         * @param part_py the array to the particle momentum (y)
         * @param part_pt the array to the particle momentum (t)
         * @param part_idcpu the array to the particle global index
         * @param ref_part the struct containing the reference particle
         */
        PushSingleParticle (T_Element element,
                            amrex::ParticleReal* AMREX_RESTRICT part_x,
                            amrex::ParticleReal* AMREX_RESTRICT part_y,
                            amrex::ParticleReal* AMREX_RESTRICT part_t,
                            amrex::ParticleReal* AMREX_RESTRICT part_px,
                            amrex::ParticleReal* AMREX_RESTRICT part_py,
                            amrex::ParticleReal* AMREX_RESTRICT part_pt,
                            uint64_t* AMREX_RESTRICT part_idcpu,
                            RefPart ref_part)
            : m_element(std::move(element)),
              m_part_x(part_x), m_part_y(part_y), m_part_t(part_t),
              m_part_px(part_px), m_part_py(part_py), m_part_pt(part_pt),
              m_part_idcpu(part_idcpu),
              m_ref_part(std::move(ref_part))
        {
        }

        PushSingleParticle () = delete;
        PushSingleParticle (PushSingleParticle const &) = default;
        PushSingleParticle (PushSingleParticle &&) = default;
        ~PushSingleParticle () = default;

        /** Push a single particle through an element
         *
         * @param i particle index in the current box
         */
        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        void
        operator() (long i) const
        {
            // access SoA Real data
            // note: an optimizing compiler will eliminate loads of unused parameters
            amrex::ParticleReal & AMREX_RESTRICT x = m_part_x[i];
            amrex::ParticleReal & AMREX_RESTRICT y = m_part_y[i];
            amrex::ParticleReal & AMREX_RESTRICT t = m_part_t[i];
            amrex::ParticleReal & AMREX_RESTRICT px = m_part_px[i];
            amrex::ParticleReal & AMREX_RESTRICT py = m_part_py[i];
            amrex::ParticleReal & AMREX_RESTRICT pt = m_part_pt[i];
            uint64_t & AMREX_RESTRICT idcpu = m_part_idcpu[i];

            // push through element
            m_element(x, y, t, px, py, pt, idcpu, m_ref_part);

        }

    private:
        T_Element const m_element;
        amrex::ParticleReal* const AMREX_RESTRICT m_part_x;
        amrex::ParticleReal* const AMREX_RESTRICT m_part_y;
        amrex::ParticleReal* const AMREX_RESTRICT m_part_t;
        amrex::ParticleReal* const AMREX_RESTRICT m_part_px;
        amrex::ParticleReal* const AMREX_RESTRICT m_part_py;
        amrex::ParticleReal* const AMREX_RESTRICT m_part_pt;
        uint64_t* const AMREX_RESTRICT m_part_idcpu;
        RefPart const m_ref_part;
    };

    /** This pushes all particles on a particle iterator tile/box
     */
    template< typename T_Element >
    void push_all_particles (
            ImpactXParticleContainer::iterator & pti,
            RefPart & AMREX_RESTRICT ref_part,
            T_Element & element
    ) {
        const int np = pti.numParticles();

        // preparing access to particle data: SoA of Reals
        auto& soa_real = pti.GetStructOfArrays().GetRealData();
        amrex::ParticleReal* const AMREX_RESTRICT part_x = soa_real[RealSoA::x].dataPtr();
        amrex::ParticleReal* const AMREX_RESTRICT part_y = soa_real[RealSoA::y].dataPtr();
        amrex::ParticleReal* const AMREX_RESTRICT part_t = soa_real[RealSoA::t].dataPtr();
        amrex::ParticleReal* const AMREX_RESTRICT part_px = soa_real[RealSoA::px].dataPtr();
        amrex::ParticleReal* const AMREX_RESTRICT part_py = soa_real[RealSoA::py].dataPtr();
        amrex::ParticleReal* const AMREX_RESTRICT part_pt = soa_real[RealSoA::pt].dataPtr();

        uint64_t* const AMREX_RESTRICT part_idcpu = pti.GetStructOfArrays().GetIdCPUData().dataPtr();

        detail::PushSingleParticle<T_Element> const pushSingleParticle(
                element, part_x, part_y, part_t, part_px, part_py, part_pt, part_idcpu, ref_part);
        //   loop over beam particles in the box
        amrex::ParallelFor(np, pushSingleParticle);
    }
} // namespace detail

    /** Mixin class for a regular beam optics lattice element.
     *
     * This element pushes particles relative to a reference particle.
     * All particles are pushed independently with the same logic.
     */
    template<typename T_Element>
    struct BeamOptic
    {
        /** Push first the reference particle, then all other particles
         *
         * @param[inout] pc container of the particles to push
         * @param[in] step global step for diagnostics
         * @param[in] period for periodic lattices, this is the current period (turn or cycle)
         */
        void operator() (
            ImpactXParticleContainer & pc,
            int step,
            int period
        )
        {
            static_assert(
                std::is_base_of_v<BeamOptic, T_Element>,
                "BeamOptic can only be used as a mixin class!"
            );

            T_Element& element = *static_cast<T_Element*>(this);
            push_all(pc, element, step, period);
        }

        /** This pushes the particles on a particle iterator tile or box.
         *
         * Particles are relative to the reference particle.
         *
         * @param[in] pti particle iterator for a current tile or box.
         * @param[in] ref_part reference particle
         */
         void operator() (
            ImpactXParticleContainer::iterator & pti,
            RefPart & AMREX_RESTRICT ref_part
         )
         {
            static_assert(
                std::is_base_of_v<BeamOptic, T_Element>,
                "BeamOptic can only be used as a mixin class!"
            );

            T_Element& element = *static_cast<T_Element*>(this);
            detail::push_all_particles<T_Element>(pti, ref_part, element);
         }
    };

} // namespace impactx::elements

#endif // IMPACTX_ELEMENTS_MIXIN_BEAMOPTIC_H
