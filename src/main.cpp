/* Copyright 2021 Axel Huebl, Chad Mitchell, Ji Qiang
 *
 * This file is part of ImpactX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "ImpactX.H"

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>


int main(int argc, char* argv[])
{
    //auto mpi_thread_levels = utils::warpx_mpi_init(argc, argv);
    //warpx_amrex_init(argc, argv);
    bool const build_parm_parse = true;
    amrex::Initialize(
        argc,
        argv,
        build_parm_parse,
        MPI_COMM_WORLD,
        impactx::overwrite_amrex_parser_defaults
    );
    //utils::warpx_check_mpi_thread_level(mpi_thread_levels);

#if defined(AMREX_USE_HIP)
    //rocfft_setup();
#endif

    //IMPACTX_PROFILE_VAR("main()", pmain);
    {
        impactx::ImpactX impactX;
        // ...
    }
    //IMPACTX_PROFILE_VAR_STOP(pmain);

#if defined(AMREX_USE_HIP)
    //rocfft_cleanup();
#endif

    amrex::Finalize();
#if defined(AMREX_USE_MPI)
    MPI_Finalize();
#endif
}
