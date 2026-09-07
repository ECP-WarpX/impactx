/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
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

#include <ablastr/parallelization/MPIInitHelpers.H>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>


int main(int argc, char* argv[])
{
    // initialize MPI with the thread support level that this build was
    // configured for, e.g., for async I/O (a no-op if built without MPI)
    ablastr::parallelization::mpi_init(argc, argv);

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

    ablastr::parallelization::mpi_finalize();
}
