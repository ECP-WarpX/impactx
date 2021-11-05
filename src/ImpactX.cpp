/* Copyright 2021 Axel Huebl, Chad Mitchell, Ji Qiang
 *
 * This file is part of ImpactX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "ImpactX.H"
#include "particles/ImpactXParticleContainer.H"

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>


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

    ImpactX::ImpactX (amrex::Geometry const& simulation_geometry, amrex::AmrInfo const& amr_info)
        : AmrCore(simulation_geometry, amr_info),
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

    /** Tag cells for refinement.  TagBoxArray tags is built on level lev grids.
     *
     * @todo this function is not (yet) implemented.
     */
    void ImpactX::ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow)
    {
        // todo
        amrex::ignore_unused(lev, tags, time, ngrow);
    }

    /** Make a new level from scratch using provided BoxArray and DistributionMapping.
     *
     * Only used during initialization.
     *
     * @todo this function is not (yet) implemented.
     */
    void ImpactX::MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
                                          const amrex::DistributionMapping& dm)
    {
        // todo data_mf.define(ba, dm, 1, 0);
        amrex::ignore_unused(lev, time, ba, dm);
    }

    /** Make a new level using provided BoxArray and DistributionMapping and fill
     *  with interpolated coarse level data.
     *
     * @todo this function is not (yet) implemented.
     */
    void ImpactX::MakeNewLevelFromCoarse (int lev, amrex::Real time, const amrex::BoxArray& ba,
                                         const amrex::DistributionMapping& dm)
    {
        // todo
        amrex::ignore_unused(lev, time, ba, dm);
    }

    /** Remake an existing level using provided BoxArray and DistributionMapping
     *  and fill with existing fine and coarse data.
     *
     * @todo this function is not (yet) implemented.
     */
    void ImpactX::RemakeLevel (int lev, amrex::Real time, const amrex::BoxArray& ba,
                              const amrex::DistributionMapping& dm)
    {
        // todo
        amrex::ignore_unused(lev, time, ba, dm);
    }

    /** Delete level data
     *
     * @todo this function is not (yet) implemented.
     */
    void ImpactX::ClearLevel (int lev)
    {
        // todo
        amrex::ignore_unused(lev);
    }

    void ImpactX::evolve (int num_steps)
    {
        BL_PROFILE("ImpactX::evolve");

        for (int step = 0; step < num_steps; ++step)
        {
            BL_PROFILE("ImpactX::evolve::step");

            using namespace amrex::literals; // for _rt and _prt

            // push all particles

            // loop over refinement levels
            int const nLevel = 0;
            for (int lev = 0; lev < nLevel; ++lev)
            {
                // get simulation geometry information
                //const amrex::Geometry& gm = this->Geom(lev);
                //const auto prob_lo = gm.ProbLo();

                // loop over all particle boxes
                using ParIt = ImpactXParticleContainer::iterator;
                for (ParIt pti(*mypc, lev); pti.isValid(); ++pti) {
                    //const auto t_lev = pti.GetLevel();
                    //const auto index = pti.GetPairIndex();
                    // ...

                    // preparing access to particle data: AoS
                    using PType = ImpactXParticleContainer::ParticleType;
                    auto& aos = pti.GetArrayOfStructs();
                    PType* AMREX_RESTRICT aos_ptr = aos().dataPtr();

                    // preparing access to particle data: SoA of Reals
                    auto& soa_real = pti.GetStructOfArrays().GetRealData();
                    amrex::ParticleReal* const AMREX_RESTRICT part_ux = soa_real[RealSoA::ux].dataPtr();
                    // ...

                    // loop over particles in the box
                    const int np = pti.numParticles();
                    amrex::ParallelFor( np, [=] AMREX_GPU_DEVICE (long i)
                    {
                        // access AoS data such as positions and cpu/id
                        PType& p = aos_ptr[i];
                        amrex::ParticleReal const x = p.pos(0);
                        //amrex::ParticleReal const y = p.pos(1);
                        //amrex::ParticleReal const z = p.pos(2);

                        // acces SoA Real data
                        amrex::ParticleReal const ux = part_ux[i];

                        // advance position
                        amrex::ParticleReal const dt = 1.0_prt;
                        p.pos(0) = x + ux * dt;
                    });
                } // end loop over all particle boxes
            } // env mesh-refinement level loop

            // do more stuff in the step

        } // end step loop
    }

} // namespace impactx
