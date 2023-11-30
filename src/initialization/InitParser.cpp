/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell, Ji Qiang
 * License: BSD-3-Clause-LBNL
 */
#include "InitParser.H"

#include <AMReX_ParmParse.H>

namespace impactx::initialization
{
    void
    overwrite_amrex_parser_defaults ()
    {
        amrex::ParmParse pp_amrex("amrex");

        // https://amrex-codes.github.io/amrex/docs_html/GPU.html#inputs-parameters
        bool abort_on_out_of_gpu_memory = true; // AMReX' default: false
        pp_amrex.queryAdd("abort_on_out_of_gpu_memory", abort_on_out_of_gpu_memory);

        bool the_arena_is_managed = false; // AMReX' default: true
        pp_amrex.queryAdd("the_arena_is_managed", the_arena_is_managed);

        // https://amrex-codes.github.io/amrex/docs_html/InputsComputeBackends.html
        std::string omp_threads = "nosmt"; // AMReX' default: system
        pp_amrex.queryAdd("omp_threads", omp_threads);

        // Here we override the default tiling option for particles, which is always
        // "false" in AMReX, to "false" if compiling for GPU execution and "true"
        // if compiling for CPU.
        {
            amrex::ParmParse pp_particles("particles");
#ifdef AMREX_USE_GPU
            bool do_tiling = false; // By default, tiling is off on GPU
#else
            bool do_tiling = true;
#endif
            pp_particles.queryAdd("do_tiling", do_tiling);
        }
    }
} // namespace impactx
