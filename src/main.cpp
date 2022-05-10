/* Copyright 2021-2022 Axel Huebl, Chad Mitchell, Ji Qiang
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

#if defined(AMREX_USE_MPI)
#   include <mpi.h>
#endif


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

    BL_PROFILE_VAR("main()", pmain);
    {
        impactx::ImpactX impactX;
        impactX.initGrids();
        impactX.initBeamDistributionFromInputs();
        impactX.initLatticeElementsFromInputs();
        impactX.evolve( /* num_steps = */ 1);
    }
    BL_PROFILE_VAR_STOP(pmain);

    amrex::Finalize();
#if defined(AMREX_USE_MPI)
    AMREX_ALWAYS_ASSERT(MPI_SUCCESS == MPI_Finalize());
#endif
}
