/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell, Ji Qiang
 * License: BSD-3-Clause-LBNL
 */
#include "ImpactX.H"
#include "initialization/InitAMReX.H"

#include <AMReX.H>
#include <AMReX_BLProfiler.H>

#if defined(AMREX_USE_MPI)
#   include <mpi.h>
#endif


int main(int argc, char* argv[])
{
#if defined(AMREX_USE_MPI)
    AMREX_ALWAYS_ASSERT(MPI_SUCCESS == MPI_Init(&argc, &argv));
#endif

    // although ImpactX' init_grids will call this if not done before, we call
    // it here so users can pass command line arguments
    impactx::initialization::default_init_AMReX(argc, argv);

    {
        BL_PROFILE_VAR("main()", pmain);
        impactx::ImpactX impactX;
        impactX.init_grids();
        impactX.initBeamDistributionFromInputs();
        impactX.initLatticeElementsFromInputs();
        impactX.evolve();
        BL_PROFILE_VAR_STOP(pmain);
        impactX.finalize();
    }

#if defined(AMREX_USE_MPI)
    AMREX_ALWAYS_ASSERT(MPI_SUCCESS == MPI_Finalize());
#endif
}
