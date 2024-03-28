/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Remi Lehe
 * License: BSD-3-Clause-LBNL
 */
#include "ImpactXParticleContainer.H"

#include "initialization/AmrCoreData.H"

#include <ablastr/constant.H>
#include <ablastr/particles/ParticleMoments.H>

#include <AMReX.H>
#include <AMReX_AmrCore.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Particle.H>

#include <algorithm>
#include <stdexcept>


namespace
{
    bool do_omp_dynamic ()
    {
        bool do_dynamic = true;
        amrex::ParmParse const pp_impactx("impactx");
        pp_impactx.query("do_dynamic_scheduling", do_dynamic);
        return do_dynamic;
    }
}

namespace impactx
{
    ParIterSoA::ParIterSoA (ContainerType& pc, int level)
        : amrex::ParIterSoA<RealSoA::nattribs, IntSoA::nattribs>(pc, level,
                   amrex::MFItInfo().SetDynamic(do_omp_dynamic())) {}

    ParIterSoA::ParIterSoA (ContainerType& pc, int level, amrex::MFItInfo& info)
        : amrex::ParIterSoA<RealSoA::nattribs, IntSoA::nattribs>(pc, level,
              info.SetDynamic(do_omp_dynamic())) {}

    ParConstIterSoA::ParConstIterSoA (ContainerType& pc, int level)
        : amrex::ParConstIterSoA<RealSoA::nattribs, IntSoA::nattribs>(pc, level,
              amrex::MFItInfo().SetDynamic(do_omp_dynamic())) {}

    ParConstIterSoA::ParConstIterSoA (ContainerType& pc, int level, amrex::MFItInfo& info)
        : amrex::ParConstIterSoA<RealSoA::nattribs, IntSoA::nattribs>(pc, level,
              info.SetDynamic(do_omp_dynamic())) {}

    ImpactXParticleContainer::ImpactXParticleContainer (initialization::AmrCoreData* amr_core)
        : amrex::ParticleContainerPureSoA<RealSoA::nattribs, IntSoA::nattribs>(amr_core->GetParGDB())
    {
        SetParticleSize();

        // name compile-time attributes
        m_real_soa_names.resize(RealSoA::names_s.size());
        m_int_soa_names.resize(IntSoA::names_s.size());
        std::copy(RealSoA::names_s.begin(), RealSoA::names_s.end(), m_real_soa_names.begin());
        std::copy(IntSoA::names_s.begin(), IntSoA::names_s.end(), m_int_soa_names.begin());
    }

    void
    ImpactXParticleContainer::AddRealComp (std::string const & name, bool communicate)
    {
        m_real_soa_names.push_back(name);
        amrex::ParticleContainerPureSoA<RealSoA::nattribs, IntSoA::nattribs>::AddRealComp(communicate);
    }

    void
    ImpactXParticleContainer::AddIntComp (std::string const & name, bool communicate)
    {
        m_int_soa_names.push_back(name);
        amrex::ParticleContainerPureSoA<RealSoA::nattribs, IntSoA::nattribs>::AddIntComp(communicate);
    }

    void
    ImpactXParticleContainer::SetLostParticleContainer (ImpactXParticleContainer * lost_pc)
    {
        m_particles_lost = lost_pc;
    }

    ImpactXParticleContainer *
    ImpactXParticleContainer::GetLostParticleContainer ()
    {
        if (m_particles_lost == nullptr)
        {
            throw std::logic_error(
                    "ImpactXParticleContainer::GetLostParticleContainer No lost particle container registered yet.");
        } else {
            return m_particles_lost;
        }
    }

    void ImpactXParticleContainer::SetParticleShape (int order) {
        if (m_particle_shape.has_value())
        {
            throw std::logic_error(
                "ImpactXParticleContainer::SetParticleShape This was already called before and cannot be changed.");
        } else
        {
            if (order < 1 || order > 3) {
                amrex::Abort("algo.particle_shape order can be only 1, 2, or 3");
            }
            m_particle_shape = order;
        }
    }

    void ImpactXParticleContainer::SetParticleShape ()
    {
        amrex::ParmParse const pp_algo("algo");
        int v = 0;
        bool const has_shape = pp_algo.query("particle_shape", v);
        if (!has_shape)
            throw std::runtime_error("particle_shape is not set, cannot initialize grids with guard cells.");
        SetParticleShape(v);
    }

    void
    ImpactXParticleContainer::AddNParticles (
        amrex::Gpu::DeviceVector<amrex::ParticleReal> const & x,
        amrex::Gpu::DeviceVector<amrex::ParticleReal> const & y,
        amrex::Gpu::DeviceVector<amrex::ParticleReal> const & t,
        amrex::Gpu::DeviceVector<amrex::ParticleReal> const & px,
        amrex::Gpu::DeviceVector<amrex::ParticleReal> const & py,
        amrex::Gpu::DeviceVector<amrex::ParticleReal> const & pt,
        amrex::ParticleReal qm,
        amrex::ParticleReal bchchg
    )
    {
        BL_PROFILE("ImpactX::AddNParticles");

        AMREX_ALWAYS_ASSERT(x.size() == y.size());
        AMREX_ALWAYS_ASSERT(x.size() == t.size());
        AMREX_ALWAYS_ASSERT(x.size() == px.size());
        AMREX_ALWAYS_ASSERT(x.size() == py.size());
        AMREX_ALWAYS_ASSERT(x.size() == pt.size());

        // number of particles to add
        int const np = x.size();

        // we add particles to lev 0, tile 0 of the first box assigned to this proc
        int lid = 0, gid = 0, tid = 0;
        {
            const auto& pmap = ParticleDistributionMap(lid).ProcessorMap();
            auto it = std::find(pmap.begin(), pmap.end(), amrex::ParallelDescriptor::MyProc());
            if (it == std::end(pmap)) {
                amrex::Abort("Attempting to add particles to box that does not exist.");
            } else {
                gid = *it;
            }
        }
        auto& particle_tile = DefineAndReturnParticleTile(lid, gid, tid);

        auto old_np = particle_tile.numParticles();
        auto new_np = old_np + np;
        particle_tile.resize(new_np);

        // Update NextID to include particles created in this function
        int pid;
#ifdef AMREX_USE_OMP
#pragma omp critical (add_beam_nextid)
#endif
        {
            pid = ParticleType::NextID();
            ParticleType::NextID(pid+np);
        }
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        static_cast<amrex::Long>(pid) + static_cast<amrex::Long>(np) < amrex::LongParticleIds::LastParticleID,
            "ERROR: overflow on particle id numbers");

        const int cpuid = amrex::ParallelDescriptor::MyProc();

        auto & soa = particle_tile.GetStructOfArrays().GetRealData();
        amrex::ParticleReal * const AMREX_RESTRICT x_arr = soa[RealSoA::x].dataPtr();
        amrex::ParticleReal * const AMREX_RESTRICT y_arr = soa[RealSoA::y].dataPtr();
        amrex::ParticleReal * const AMREX_RESTRICT t_arr = soa[RealSoA::t].dataPtr();
        amrex::ParticleReal * const AMREX_RESTRICT px_arr = soa[RealSoA::px].dataPtr();
        amrex::ParticleReal * const AMREX_RESTRICT py_arr = soa[RealSoA::py].dataPtr();
        amrex::ParticleReal * const AMREX_RESTRICT pt_arr = soa[RealSoA::pt].dataPtr();
        amrex::ParticleReal * const AMREX_RESTRICT qm_arr = soa[RealSoA::qm].dataPtr();
        amrex::ParticleReal * const AMREX_RESTRICT w_arr  = soa[RealSoA::w ].dataPtr();

        uint64_t * const AMREX_RESTRICT idcpu_arr = particle_tile.GetStructOfArrays().GetIdCPUData().dataPtr();

        amrex::ParticleReal const * const AMREX_RESTRICT x_ptr = x.data();
        amrex::ParticleReal const * const AMREX_RESTRICT y_ptr = y.data();
        amrex::ParticleReal const * const AMREX_RESTRICT t_ptr = t.data();
        amrex::ParticleReal const * const AMREX_RESTRICT px_ptr = px.data();
        amrex::ParticleReal const * const AMREX_RESTRICT py_ptr = py.data();
        amrex::ParticleReal const * const AMREX_RESTRICT pt_ptr = pt.data();

        amrex::ParallelFor(np,
        [=] AMREX_GPU_DEVICE (int i) noexcept
        {
            idcpu_arr[old_np+i] = amrex::SetParticleIDandCPU(pid + i, cpuid);

            x_arr[old_np+i] = x_ptr[i];
            y_arr[old_np+i] = y_ptr[i];
            t_arr[old_np+i] = t_ptr[i];

            px_arr[old_np+i] = px_ptr[i];
            py_arr[old_np+i] = py_ptr[i];
            pt_arr[old_np+i] = pt_ptr[i];
            qm_arr[old_np+i] = qm;
            w_arr[old_np+i]  = bchchg/ablastr::constant::SI::q_e/np;
        });

        // safety first: in case passed attribute arrays were temporary, we
        // want to make sure the ParallelFor has ended here
        amrex::Gpu::streamSynchronize();
    }

    void
    ImpactXParticleContainer::SetRefParticle (RefPart const & refpart)
    {
        m_refpart = refpart;
    }

    RefPart &
    ImpactXParticleContainer::GetRefParticle ()
    {
        return m_refpart;
    }

    RefPart const &
    ImpactXParticleContainer::GetRefParticle () const
    {
        return m_refpart;
    }

    void
    ImpactXParticleContainer::SetRefParticleEdge ()
    {
        m_refpart.sedge = m_refpart.s;
    }

    std::tuple<
            amrex::ParticleReal, amrex::ParticleReal,
            amrex::ParticleReal, amrex::ParticleReal,
            amrex::ParticleReal, amrex::ParticleReal>
    ImpactXParticleContainer::MinAndMaxPositions ()
    {
        BL_PROFILE("ImpactXParticleContainer::MinAndMaxPositions");
        return ablastr::particles::MinAndMaxPositions(*this);
    }

    std::tuple<
            amrex::ParticleReal, amrex::ParticleReal,
            amrex::ParticleReal, amrex::ParticleReal,
            amrex::ParticleReal, amrex::ParticleReal>
    ImpactXParticleContainer::MeanAndStdPositions ()
    {
        BL_PROFILE("ImpactXParticleContainer::MeanAndStdPositions");
        return ablastr::particles::MeanAndStdPositions<
            ImpactXParticleContainer, RealSoA::w
        >(*this);
    }

    std::vector<std::string>
    ImpactXParticleContainer::RealSoA_names () const
    {
        return m_real_soa_names;
    }

    std::vector<std::string>
    ImpactXParticleContainer::intSoA_names () const
    {
        return m_int_soa_names;
    }

    CoordSystem
    ImpactXParticleContainer::GetCoordSystem () const
    {
        return m_coordsystem;
    }

    void
    ImpactXParticleContainer::SetCoordSystem (CoordSystem coord_system)
    {
        m_coordsystem = coord_system;
    }
} // namespace impactx
