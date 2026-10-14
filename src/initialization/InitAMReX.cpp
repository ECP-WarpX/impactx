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

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>

#if defined(AMREX_USE_MPI)
#   include <mpi.h>
#else
#   include <AMReX_ccse-mpi.H>
using amrex::mpidatatypes::MPI_COMM_WORLD;
#endif

#include <string>


namespace impactx::initialization
{
    void
    default_init_AMReX (int argc, char* argv[])
    {
        if (!amrex::Initialized())
        {
            bool const build_parm_parse = true;
            amrex::Initialize(
                    argc,
                    argv,
                    build_parm_parse,
                    MPI_COMM_WORLD,
                    impactx::initialization::overwrite_amrex_parser_defaults
            );
        }
    }

    void
    default_init_AMReX ()
    {
        default_init_AMReX(0, nullptr);
    }

    void
    load_inputs_file (std::string const & filename)
    {
        if (amrex::Initialized())
        {
            amrex::ParmParse::addfile(filename);
            return;
        }

        // AMReX and ImpactX read some of their parameters while initializing, e.g.,
        // amrex.omp_threads: as for an inputs file of the executable, the file comes first
        int argc = 0;
        char** argv = nullptr;
        bool const build_parm_parse = true;
        amrex::Initialize(
                argc,
                argv,
                build_parm_parse,
                MPI_COMM_WORLD,
                [&filename]()
                {
                    amrex::ParmParse::addfile(filename);
                    impactx::initialization::overwrite_amrex_parser_defaults();
                }
        );
    }
} // namespace impactx::initialization
