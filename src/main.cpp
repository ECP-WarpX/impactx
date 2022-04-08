/* Copyright 2021 Axel Huebl, Chad Mitchell, Ji Qiang
 *
 * This file is part of ImpactX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "ImpactX.H"
#include "initialization/InitParser.H"

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>

#include <memory>

int main(int argc, char* argv[])
{
#if defined(AMREX_USE_MPI)
    AMREX_ALWAYS_ASSERT(MPI_SUCCESS == MPI_Init(&argc, &argv));
#endif

    bool const build_parm_parse = true;
    amrex::Initialize(
        argc,
        argv,
        build_parm_parse,
        MPI_COMM_WORLD,
        impactx::initialization::overwrite_amrex_parser_defaults
    );

    {
        BL_PROFILE_VAR("main()", pmain);

        amrex::AmrInfo amr_info;
        const int nprocs = amrex::ParallelDescriptor::NProcs();
        const amrex::IntVect high_end = amr_info.blocking_factor[0]
            * amrex::IntVect(AMREX_D_DECL(nprocs,1,1)) - amrex::IntVect(1);
        amrex::Box domain(amrex::IntVect(0), high_end); // Domain index space
        amrex::RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.,1.,1.)});// Domain physical size
        amrex::Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};
        amrex::Geometry geom(domain, rb, amrex::CoordSys::cartesian, is_periodic);
        auto impactX = std::make_unique<impactx::ImpactX>(geom, amr_info);

        impactX->initData();
        impactX->initElements();
        impactX->evolve( /* num_steps = */ 1);

        BL_PROFILE_VAR_STOP(pmain);
    }

    amrex::Finalize();
#if defined(AMREX_USE_MPI)
    AMREX_ALWAYS_ASSERT(MPI_SUCCESS == MPI_Finalize());
#endif
}
