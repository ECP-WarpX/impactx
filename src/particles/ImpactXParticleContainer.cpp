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


namespace impactx
{
    bool do_omp_dynamic () {
        bool do_dynamic = true;
        amrex::ParmParse pp_impactx("impactx");
        pp_impactx.query("do_dynamic_scheduling", do_dynamic);
        return do_dynamic;
    }

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

    void ImpactXParticleContainer::SetParticleShape (int const order) {
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
        amrex::ParmParse pp_algo("algo");
        int v = 0;
        bool has_shape = pp_algo.query("particle_shape", v);
        if (!has_shape)
            throw std::runtime_error("particle_shape is not set, cannot initialize grids with guard cells.");
        SetParticleShape(v);
    }

    void
    ImpactXParticleContainer::AddNParticles (int lev,
                                             amrex::Vector<amrex::ParticleReal> const & x,
                                             amrex::Vector<amrex::ParticleReal> const & y,
                                             amrex::Vector<amrex::ParticleReal> const & z,
                                             amrex::Vector<amrex::ParticleReal> const & px,
                                             amrex::Vector<amrex::ParticleReal> const & py,
                                             amrex::Vector<amrex::ParticleReal> const & pz,
                                             amrex::ParticleReal const & qm,
                                             amrex::ParticleReal const & bchchg)
    {
        BL_PROFILE("ImpactX::AddNParticles");

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(lev == 0, "AddNParticles: only lev=0 is supported yet.");
        AMREX_ALWAYS_ASSERT(x.size() == y.size());
        AMREX_ALWAYS_ASSERT(x.size() == z.size());
        AMREX_ALWAYS_ASSERT(x.size() == px.size());
        AMREX_ALWAYS_ASSERT(x.size() == py.size());
        AMREX_ALWAYS_ASSERT(x.size() == pz.size());

        // number of particles to add
        int const np = x.size();

        // have to resize here, not in the constructor because grids have not
        // been built when constructor was called.
        reserveData();
        resizeData();

        auto& particle_tile = DefineAndReturnParticleTile(0, 0, 0);

        /* Create a temporary tile to obtain data from simulation. This data
         * is then copied to the permanent tile which is stored on the particle
         * (particle_tile).
         */
        using PinnedTile = amrex::ParticleTile<
            amrex::Particle<NStructReal, NStructInt>,
            NArrayReal, NArrayInt,
            amrex::PinnedArenaAllocator
        >;
        PinnedTile pinned_tile;
        pinned_tile.define(NumRuntimeRealComps(), NumRuntimeIntComps());

        for (int i = 0; i < np; i++)
        {
            ParticleType p;
            p.id() = ParticleType::NextID();
            p.cpu() = amrex::ParallelDescriptor::MyProc();
            p.pos(0) = x[i];
            p.pos(1) = y[i];
            p.pos(2) = z[i];
            // write position, creating cpu id, and particle id
            pinned_tile.push_back(p);
        }

        // write Real attributes (SoA) to particle initialized zero
        DefineAndReturnParticleTile(0, 0, 0);

        pinned_tile.push_back_real(RealSoA::ux, px);
        pinned_tile.push_back_real(RealSoA::uy, py);
        pinned_tile.push_back_real(RealSoA::pt, pz);
        pinned_tile.push_back_real(RealSoA::m_qm, np, qm);
        pinned_tile.push_back_real(RealSoA::w, np, bchchg/ablastr::constant::SI::q_e/np);

        /* Redistributes particles to their respective tiles (spatial bucket
         * sort per box over MPI ranks)
         */
        auto old_np = particle_tile.numParticles();
        auto new_np = old_np + pinned_tile.numParticles();
        particle_tile.resize(new_np);
        amrex::copyParticles(
                particle_tile, pinned_tile, 0, old_np, pinned_tile.numParticles());
    }

    void
    ImpactXParticleContainer::SetRefParticle (RefPart const refpart)
    {
        m_refpart = refpart;
    }    RefPart &
    ImpactXParticleContainer::GetRefParticle ()
    {
        return m_refpart;
    }

    RefPart const
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
} // namespace impactx
