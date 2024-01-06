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

#include <ablastr/constant.H>
#include <ablastr/particles/ParticleMoments.H>

#include <AMReX.H>
#include <AMReX_AmrCore.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParticleTile.H>

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
    ParIter::ParIter (ContainerType& pc, int level)
        : amrex::ParIter<0, 0, RealSoA::nattribs, IntSoA::nattribs>(pc, level,
                   amrex::MFItInfo().SetDynamic(do_omp_dynamic())) {}

    ParIter::ParIter (ContainerType& pc, int level, amrex::MFItInfo& info)
        : amrex::ParIter<0, 0, RealSoA::nattribs, IntSoA::nattribs>(pc, level,
              info.SetDynamic(do_omp_dynamic())) {}

    ParConstIter::ParConstIter (ContainerType& pc, int level)
        : amrex::ParConstIter<0, 0, RealSoA::nattribs, IntSoA::nattribs>(pc, level,
              amrex::MFItInfo().SetDynamic(do_omp_dynamic())) {}

    ParConstIter::ParConstIter (ContainerType& pc, int level, amrex::MFItInfo& info)
        : amrex::ParConstIter<0, 0, RealSoA::nattribs, IntSoA::nattribs>(pc, level,
              info.SetDynamic(do_omp_dynamic())) {}

    ImpactXParticleContainer::ImpactXParticleContainer (amrex::AmrCore* amr_core)
        : amrex::ParticleContainer<0, 0, RealSoA::nattribs, IntSoA::nattribs>(amr_core->GetParGDB())
    {
        SetParticleSize();
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

        auto& particle_tile = DefineAndReturnParticleTile(0, 0, 0);
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

        auto * AMREX_RESTRICT pstructs = particle_tile.GetArrayOfStructs()().dataPtr();
        auto & soa = particle_tile.GetStructOfArrays().GetRealData();
        amrex::ParticleReal * const AMREX_RESTRICT px_arr = soa[RealSoA::px].dataPtr();
        amrex::ParticleReal * const AMREX_RESTRICT py_arr = soa[RealSoA::py].dataPtr();
        amrex::ParticleReal * const AMREX_RESTRICT pt_arr = soa[RealSoA::pt].dataPtr();
        amrex::ParticleReal * const AMREX_RESTRICT qm_arr = soa[RealSoA::qm].dataPtr();
        amrex::ParticleReal * const AMREX_RESTRICT w_arr  = soa[RealSoA::w ].dataPtr();

        amrex::ParticleReal const * const AMREX_RESTRICT x_ptr = x.data();
        amrex::ParticleReal const * const AMREX_RESTRICT y_ptr = y.data();
        amrex::ParticleReal const * const AMREX_RESTRICT t_ptr = t.data();
        amrex::ParticleReal const * const AMREX_RESTRICT px_ptr = px.data();
        amrex::ParticleReal const * const AMREX_RESTRICT py_ptr = py.data();
        amrex::ParticleReal const * const AMREX_RESTRICT pt_ptr = pt.data();

        amrex::ParallelFor(np,
        [=] AMREX_GPU_DEVICE (int i) noexcept
        {
            ParticleType& p = pstructs[old_np + i];
            p.id() = pid + i;
            p.cpu() = cpuid;
            p.pos(RealAoS::x) = x_ptr[i];
            p.pos(RealAoS::y) = y_ptr[i];
            p.pos(RealAoS::t) = t_ptr[i];

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
    ImpactXParticleContainer::RealAoS_names () const
    {
        return get_RealAoS_names();
    }

    std::vector<std::string>
    ImpactXParticleContainer::RealSoA_names () const
    {
        return get_RealSoA_names(this->NumRealComps());
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

    std::vector<std::string>
    get_RealAoS_names ()
    {
        std::vector<std::string> real_aos_names(RealAoS::names_s.size());

        // compile-time attributes
        std::copy(RealAoS::names_s.begin(), RealAoS::names_s.end(), real_aos_names.begin());

        return real_aos_names;
    }

    std::vector<std::string>
    get_RealSoA_names (int num_real_comps)
    {
        std::vector<std::string> real_soa_names(num_real_comps);

        // compile-time attributes
        std::copy(RealSoA::names_s.begin(), RealSoA::names_s.end(), real_soa_names.begin());

        // runtime attributes
        if (num_real_comps > int(RealSoA::names_s.size()))
        {
            // particles lost record their "s" position where they got lost
            real_soa_names[RealSoA::nattribs] = "s_lost";
        }

        return real_soa_names;
    }
} // namespace impactx
