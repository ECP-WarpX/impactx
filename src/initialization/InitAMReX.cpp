/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell, Ji Qiang
 * License: BSD-3-Clause-LBNL
 */
#include "InitAMReX.H"

#include "initialization/InitParser.H"

#include <ablastr/parallelization/MPIInitHelpers.H>

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>

#if defined(AMREX_USE_MPI)
#   include <mpi.h>
#else
#   include <AMReX_ccse-mpi.H>
using amrex::mpidatatypes::MPI_COMM_WORLD;
#endif

#include <cstdlib>


namespace impactx::initialization
{
namespace
{
    /** Did ImpactX initialize AMReX?
     *
     * Only then do we finalize AMReX in ImpactX::finalize(). Otherwise, e.g.,
     * pyAMReX or another simulation code in the same process owns AMReX.
     */
    bool s_amrex_initialized_by_impactx = false;

#if defined(AMREX_USE_MPI)
    /** Finalize MPI at process exit
     *
     * This is only registered if ImpactX initialized MPI itself.
     */
    void
    finalize_MPI_atexit ()
    {
        int is_finalized = 0;
        MPI_Finalized(&is_finalized);
        if (!is_finalized)
        {
            ablastr::parallelization::mpi_finalize();
        }
    }

    /** Initialize MPI, if not already done
     *
     * We initialize MPI here instead of letting amrex::Initialize() do it,
     * because AMReX finalizes MPI in amrex::Finalize() if and only if it
     * initialized MPI itself. MPI cannot be re-initialized after MPI_Finalize,
     * which would make an ImpactX simulation the last one in its process.
     * Initializing MPI here and finalizing it at process exit keeps MPI alive
     * over arbitrarily many amrex::Initialize/amrex::Finalize cycles, e.g., for
     * parameter scans and optimization loops that create a new ImpactX
     * simulation per iteration.
     */
    void
    init_MPI (int argc, char* argv[])
    {
        int is_initialized = 0;
        MPI_Initialized(&is_initialized);
        if (is_initialized)
        {
            // MPI is owned by the user, e.g., by main() or by mpi4py
            return;
        }

        ablastr::parallelization::mpi_init(argc, argv);

        std::atexit(finalize_MPI_atexit);
    }
#endif
} // namespace <anonymous>

    void
    default_init_AMReX (int argc, char* argv[])
    {
        if (amrex::Initialized())
        {
            // AMReX is owned by the user, e.g., by pyAMReX or another
            // simulation code in the same process
            return;
        }

        // amrex::Initialize() only initializes its runtime parameter database
        // (ParmParse) if a command line is passed (argc >= 1). Without
        // ParmParse::Initialize, ParmParse also never registers its
        // ParmParse::Finalize, i.e., its parameter table is never cleared in
        // amrex::Finalize() and the inputs of a simulation leak into the next
        // simulation of the same process. Thus, pass a dummy program name if we
        // have no command line.
        // Note: with argc == 1, AMReX passes no command line arguments on to
        //       ParmParse, it only uses argv[0] as the executable name.
        char program_name[] = "impactx";
        char* dummy_argv[] = {program_name, nullptr};
        if (argc < 1)
        {
            argc = 1;
            argv = dummy_argv;
        }

#if defined(AMREX_USE_MPI)
        init_MPI(argc, argv);
#endif

        bool const build_parm_parse = true;
        amrex::Initialize(
                argc,
                argv,
                build_parm_parse,
                MPI_COMM_WORLD,
                impactx::initialization::overwrite_amrex_parser_defaults
        );

        s_amrex_initialized_by_impactx = true;
    }

    void
    default_init_AMReX ()
    {
        default_init_AMReX(0, nullptr);
    }

    bool
    default_finalize_AMReX ()
    {
        if (!s_amrex_initialized_by_impactx)
        {
            return false;
        }

        // reset first: if AMReX was already finalized by the user, e.g., via
        // pyAMReX, we do not own it anymore either
        s_amrex_initialized_by_impactx = false;

        if (!amrex::Initialized())
        {
            return false;
        }

        amrex::Finalize();

        return true;
    }
} // namespace impactx::initialization
