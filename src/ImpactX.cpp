/* Copyright 2021 Axel Huebl, Chad Mitchell, Ji Qiang
 *
 * This file is part of ImpactX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "ImpactX.H"

#include <AMReX_ParmParse.H>

namespace impactx
{
    void
    overwrite_amrex_parser_defaults()
    {
        amrex::ParmParse pp_amrex("amrex");

        // https://amrex-codes.github.io/amrex/docs_html/GPU.html#inputs-parameters
        bool abort_on_out_of_gpu_memory = true; // AMReX' default: false
        pp_amrex.query("abort_on_out_of_gpu_memory", abort_on_out_of_gpu_memory);
        pp_amrex.add("abort_on_out_of_gpu_memory", abort_on_out_of_gpu_memory);
    }

    ImpactX::ImpactX (amrex::Geometry const& geom, amrex::AmrInfo const& amr_info)
        : AmrCore(geom, amr_info),
          mypc(std::make_unique<ImpactXParticleContainer>(this))
    {
    }

    void ImpactX::initData ()
    {
        AmrCore::InitFromScratch(0.0);
        amrex::Print() << "boxArray(0) " << boxArray(0) << std::endl;;

        mypc->AddNParticles(0, {0.0}, {0.2}, {0.4});
        amrex::Print() << "# of particles: " << mypc->TotalNumberOfParticles() << std::endl;
    }

    //! Tag cells for refinement.  TagBoxArray tags is built on level lev grids.
    void ImpactX::ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow)
    {
        // todo
    }

    //! Make a new level from scratch using provided BoxArray and DistributionMapping.
    //! Only used during initialization.
    void ImpactX::MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
                                          const amrex::DistributionMapping& dm)
    {
        // todo data_mf.define(ba, dm, 1, 0);
    }

    //! Make a new level using provided BoxArray and DistributionMapping and fill
    //  with interpolated coarse level data.
    void ImpactX::MakeNewLevelFromCoarse (int lev, amrex::Real time, const amrex::BoxArray& ba,
                                         const amrex::DistributionMapping& dm)
    {
        amrex::Print() << "MakeNewLevelFromCoarse" << std::endl;
        // todo
    }

    //! Remake an existing level using provided BoxArray and DistributionMapping
    //  and fill with existing fine and coarse data.
    void ImpactX::RemakeLevel (int lev, amrex::Real time, const amrex::BoxArray& ba,
                              const amrex::DistributionMapping& dm)
    {
        // todo
    }

    //! Delete level data
    void ImpactX::ClearLevel (int lev)
    {
        // todo
    }
} // namespace impactx
