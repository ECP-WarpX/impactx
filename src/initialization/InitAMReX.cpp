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

#include <cerrno>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <functional>
#include <string>
#include <system_error>


namespace impactx::initialization
{
namespace
{
    /** Initialize AMReX with the ImpactX defaults for AMReX parameters
     *
     * @param argc number of command line arguments
     * @param argv command line arguments
     * @param read_inputs reads parameters before AMReX and ImpactX look up their
     *                    runtime parameters, e.g., amrex.omp_threads (optional)
     */
    void
    init_AMReX (int argc, char* argv[], std::function<void()> const & read_inputs)
    {
        bool const build_parm_parse = true;
        amrex::Initialize(
                argc,
                argv,
                build_parm_parse,
                MPI_COMM_WORLD,
                [&read_inputs]()
                {
                    if (read_inputs) { read_inputs(); }
                    impactx::initialization::overwrite_amrex_parser_defaults();
                }
        );
    }

    /** Throw unless the IO process can read the inputs file
     *
     * AMReX checks and reads the file on the IO process only and broadcasts its content,
     * so the result of this check is broadcast, too: all processes throw, or none.
     *
     * @param filename path to the inputs file
     */
    void
    check_inputs_file (std::string const & filename)
    {
        int error = 0;
        if (amrex::ParallelDescriptor::IOProcessor())
        {
            std::error_code ec;
            auto const type = std::filesystem::status(filename, ec).type();
            if (type == std::filesystem::file_type::not_found) {
                error = ENOENT;
            } else if (ec) {
                error = ec.default_error_condition().value();
            } else if (type == std::filesystem::file_type::directory) {
                error = EISDIR;
            } else if (std::FILE* file = std::fopen(filename.c_str(), "r"); file == nullptr) {
                error = errno != 0 ? errno : EACCES;
            } else {
                std::fclose(file);
            }
        }
        amrex::ParallelDescriptor::Bcast(
            &error, 1, amrex::ParallelDescriptor::IOProcessorNumber());

        if (error != 0)
        {
            throw std::filesystem::filesystem_error(
                "Cannot read the inputs file",
                filename,
                std::error_code(error, std::generic_category())
            );
        }
    }
} // namespace

    void
    default_init_AMReX (int argc, char* argv[])
    {
        if (!amrex::Initialized())
        {
            init_AMReX(argc, argv, {});
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
            check_inputs_file(filename);
            amrex::ParmParse::addfile(filename);
            return;
        }

        // An error in the file must neither escape amrex::Initialize, which would leave AMReX
        // half initialized, nor abort the process: keep it, then shut AMReX down and throw it.
        std::exception_ptr error;
        init_AMReX(0, nullptr,
            [&filename, &error]()
            {
                // AMReX reads amrex.throw_exception only after this function
                amrex::system::throw_exception = true;
                try
                {
                    check_inputs_file(filename);
                    amrex::ParmParse::addfile(filename);
                }
                catch (...)
                {
                    error = std::current_exception();
                }
            }
        );

        if (error)
        {
            // parameters of a file that was not read to its end stay unused
            amrex::system::abort_on_unused_inputs = false;
            amrex::Finalize();
            std::rethrow_exception(error);
        }
    }
} // namespace impactx::initialization
