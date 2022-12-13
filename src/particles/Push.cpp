/* Copyright 2022 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "Push.H"

#include <AMReX_BLProfiler.H>
#include <AMReX_Extension.H>  // for AMREX_RESTRICT
#include <AMReX_REAL.H>       // for ParticleReal

#include <utility>


namespace impactx
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

        /** Constructor taking in pointers to particle data
         *
         * @param element the beamline element to push through
         * @param aos_ptr the array-of-struct with position and ids
         * @param part_px the array to the particle momentum (x)
         * @param part_py the array to the particle momentum (y)
         * @param part_pt the array to the particle momentum (t)
         * @param ref_part the struct containing the reference particle
         */
        PushSingleParticle (T_Element element,
                            PType* AMREX_RESTRICT aos_ptr,
                            amrex::ParticleReal* AMREX_RESTRICT part_px,
                            amrex::ParticleReal* AMREX_RESTRICT part_py,
                            amrex::ParticleReal* AMREX_RESTRICT part_pt,
                            RefPart ref_part)
            : m_element(std::move(element)), m_aos_ptr(aos_ptr),
              m_part_px(part_px), m_part_py(part_py), m_part_pt(part_pt),
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
            // access AoS data such as positions and cpu/id
            PType& AMREX_RESTRICT p = m_aos_ptr[i];

            // access SoA Real data
            amrex::ParticleReal & AMREX_RESTRICT px = m_part_px[i];
            amrex::ParticleReal & AMREX_RESTRICT py = m_part_py[i];
            amrex::ParticleReal & AMREX_RESTRICT pt = m_part_pt[i];

            // push through element
            m_element(p, px, py, pt, m_ref_part);

        }

    private:
        T_Element const m_element;
        PType* const AMREX_RESTRICT m_aos_ptr;
        amrex::ParticleReal* const AMREX_RESTRICT m_part_px;
        amrex::ParticleReal* const AMREX_RESTRICT m_part_py;
        amrex::ParticleReal* const AMREX_RESTRICT m_part_pt;
        RefPart const m_ref_part;
    };
} // namespace detail

    void Push (ImpactXParticleContainer & pc,
               KnownElements const & element_variant)
    {
        // performance profiling per element
        std::string element_name;
        std::visit([&element_name](auto&& element){ element_name = element.name; }, element_variant);
        std::string const profile_name = "impactx::Push::" + element_name;
        BL_PROFILE("impactx::Push");
        BL_PROFILE(profile_name);

        using namespace amrex::literals; // for _rt and _prt

        // preparing to access reference particle data: RefPart
        RefPart & ref_part = pc.GetRefParticle();

        // loop over refinement levels
        int const nLevel = pc.finestLevel();
        for (int lev = 0; lev <= nLevel; ++lev)
        {
            // get simulation geometry information
            //const amrex::Geometry& gm = this->Geom(lev);
            //const auto prob_lo = gm.ProbLo();

            // loop over all particle boxes
            using ParIt = ImpactXParticleContainer::iterator;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (ParIt pti(pc, lev); pti.isValid(); ++pti) {
                const int np = pti.numParticles();
                //const auto t_lev = pti.GetLevel();
                //const auto index = pti.GetPairIndex();
                // ...

                // preparing access to particle data: AoS
                using PType = ImpactXParticleContainer::ParticleType;
                auto& aos = pti.GetArrayOfStructs();
                PType* AMREX_RESTRICT aos_ptr = aos().dataPtr();

                // preparing access to particle data: SoA of Reals
                auto& soa_real = pti.GetStructOfArrays().GetRealData();
                amrex::ParticleReal* const AMREX_RESTRICT part_px = soa_real[RealSoA::ux].dataPtr();
                amrex::ParticleReal* const AMREX_RESTRICT part_py = soa_real[RealSoA::uy].dataPtr();
                amrex::ParticleReal* const AMREX_RESTRICT part_pt = soa_real[RealSoA::pt].dataPtr();

                // here we just access the element by its respective type
                std::visit(
                    [=, &ref_part](auto element) {
                        // push beam particles relative to reference particle
                        detail::PushSingleParticle<decltype(element)> const pushSingleParticle(
                            element, aos_ptr, part_px, part_py, part_pt, ref_part);
                        //   loop over beam particles in the box
                        amrex::ParallelFor(np, pushSingleParticle);

                        // push reference particle in global coordinates
                        element(ref_part);
                    },
                    element_variant
                );
            } // end loop over all particle boxes
        } // env mesh-refinement level loop
    }

} // namespace impactx
