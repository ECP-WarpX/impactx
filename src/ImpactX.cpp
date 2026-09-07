/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell, Ji Qiang, Remi Lehe
 * License: BSD-3-Clause-LBNL
 */
#include "ImpactX.H"
#include "diagnostics/DiagnosticOutput.H"
#include "diagnostics/FilePrefix.H"
#include "elements/mixin/accessors.H"
#include "elements/mixin/dynamicdata.H"
#include "initialization/InitAMReX.H"
#include "initialization/InitAmrCore.H"
#include "particles/ImpactXParticleContainer.H"
#include "particles/Push.H"
#include "particles/wakefields/HandleWakefield.H"

#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>

#include <iostream>
#include <memory>
#include <string>


namespace impactx {
    ImpactX::ImpactX() {
        // todo: if amr.n_cells is provided, overwrite/redefine AmrCore object

        // todo: if charge deposition and/or space charge are requested, require
        //       amr.n_cells from user inputs
    }

    ImpactX::~ImpactX()
    {
        this->finalize();
    }

    void ImpactX::finalize ()
    {
        // loop over all beamline elements & finalize them
        finalize_elements();

        if (m_grids_initialized)
        {
            m_lattice.clear();

            // this one last
            amr_data.reset();

            // only finalize AMReX if we initialized it: another library or the
            // user (e.g., via pyAMReX) might own AMReX in this process
            initialization::default_finalize_AMReX();

            // only finalize once
            m_grids_initialized = false;
        }
    }

    void ImpactX::finalize_elements ()
    {
        // loop over all beamline elements & finalize them
        for (auto & element_variant : m_lattice)
        {
            elements::finalize(element_variant);
        }
    }

    void ImpactX::init_grids ()
    {
        BL_PROFILE("ImpactX::init_grids");

        amr_data = std::make_unique<initialization::AmrCoreData>(initialization::init_amr_core());
        amr_data->track_particles.m_particle_container = std::make_unique<ImpactXParticleContainer>(amr_data.get());
        amr_data->track_particles.m_particles_lost = std::make_unique<ImpactXParticleContainer>(amr_data.get());

        // query input for warning logger variables and set up warning logger accordingly
        init_warning_logger();

        // move old diagnostics out of the way
        bool diag_enable = true;
        amrex::ParmParse("diag").queryAdd("enable", diag_enable);
        if (diag_enable && diagnostics::FilePrefixCanClean()) {
            amrex::UtilCreateCleanDirectory(diagnostics::FilePrefix(), true);
        }

        // the particle container has been set to track the same Geometry as ImpactX

        // this is the earliest point that we need to know the particle shape,
        // so that we can initialize the guard size of our MultiFabs
        amr_data->track_particles.m_particle_container->SetParticleShape();

        // init blocks / grids & MultiFabs
        amr_data->InitFromScratch(0.0);

        // prepare particle containers
        //   have to do this here, not in the constructor because grids have not
        //   been built when constructor was called.
        amr_data->track_particles.m_particle_container->prepare();
        amr_data->track_particles.m_particles_lost->prepare();

        // register shortcut
        amr_data->track_particles.m_particle_container->SetLostParticleContainer(amr_data->track_particles.m_particles_lost.get());

        // print AMReX grid summary
        if (amrex::ParallelDescriptor::IOProcessor()) {
            // verbosity
            amrex::ParmParse pp_impactx("impactx");
            int verbose = 1;
            pp_impactx.queryAddWithParser("verbose", verbose);

            if (verbose > 0) {
                std::cout << "\nGrids Summary:\n";
                amr_data->printGridSummary(std::cout, 0, amr_data->finestLevel());
            }
        }

        // warn on inefficient domain decompositions (number of boxes vs. parallel processes)
        {
            std::string const python_hint(
                " In Python, the equivalent options are sim.max_grid_size and "
                "sim.max_grid_size_x/_y/_z."
            );
            int const nprocs = amrex::ParallelDescriptor::NProcs();
            for (int lev = 0; lev <= amr_data->finestLevel(); ++lev) {
                amrex::Long const nboxes = amr_data->boxArray(lev).size();
                if (nprocs == 1 && nboxes > 1) {
                    ablastr::warn_manager::WMRecordWarning(
                        "ImpactX::init_grids",
                        "The mesh on level " + std::to_string(lev) + " is decomposed into " +
                        std::to_string(nboxes) + " boxes, but this simulation runs on a "
                        "single process. More than one box per process causes significant "
                        "overhead (particle redistribution, communication in field solves, "
                        "extra kernel launches) without any parallelization benefit. "
                        "Set amr.max_grid_size at least as large as amr.n_cell "
                        "(per direction: amr.max_grid_size_x/_y/_z) to create a single box, "
                        "or run with more MPI processes." + python_hint,
                        ablastr::warn_manager::WarnPriority::high
                    );
                } else if (nboxes > nprocs) {
                    ablastr::warn_manager::WMRecordWarning(
                        "ImpactX::init_grids",
                        "The mesh on level " + std::to_string(lev) + " is decomposed into " +
                        std::to_string(nboxes) + " boxes for only " + std::to_string(nprocs) +
                        " MPI processes. More than one box per process (e.g., per GPU) "
                        "causes significant overhead (particle redistribution, communication "
                        "in field solves, extra kernel launches). Increase amr.max_grid_size "
                        "(per direction: amr.max_grid_size_x/_y/_z) so that the number of "
                        "boxes matches the number of processes." + python_hint,
                        ablastr::warn_manager::WarnPriority::high
                    );
                } else if (lev == 0 && nboxes < nprocs) {
                    ablastr::warn_manager::WMRecordWarning(
                        "ImpactX::init_grids",
                        "The mesh on level 0 is decomposed into " + std::to_string(nboxes) +
                        " box(es) for " + std::to_string(nprocs) + " MPI processes. Processes "
                        "without a box do not participate in field solves and hold no "
                        "particles after redistribution. Decrease amr.max_grid_size "
                        "(per direction: amr.max_grid_size_x/_y/_z) so that the number of "
                        "boxes matches the number of processes." + python_hint,
                        ablastr::warn_manager::WarnPriority::high
                    );
                }
            }
        }

        // keep track that init is done
        m_grids_initialized = true;
    }

    void ImpactX::evolve ()
    {
        BL_PROFILE("ImpactX::evolve");

        amrex::ParmParse pp_algo("algo");
        std::string track = "particles";
        pp_algo.queryAdd("track", track);

        if (track == "particles") {
            track_particles();
        }
        else if (track == "envelope") {
            track_envelope();
        }
        else if (track == "reference_orbit") {
            if (!amr_data->track_reference.m_ref.has_value())
            {
                throw std::runtime_error("evolve: Reference particle not set.");
            }
            track_reference(amr_data->track_reference.m_ref.value());
        }
        else {
            throw std::runtime_error("Unknown tracking algorithm: algo.track=" + track);
        }
    }

} // namespace impactx
