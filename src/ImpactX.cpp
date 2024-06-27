/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell, Ji Qiang, Remi Lehe
 * License: BSD-3-Clause-LBNL
 */
#include "ImpactX.H"
#include "initialization/InitAmrCore.H"
#include "initialization/InitDistribution.H"
#include "particles/CollectLost.H"
#include "particles/ImpactXParticleContainer.H"
#include "particles/Push.H"
#include "particles/diagnostics/DiagnosticOutput.H"
#include "particles/spacecharge/ForceFromSelfFields.H"
#include "particles/spacecharge/GatherAndPush.H"
#include "particles/spacecharge/PoissonSolve.H"
#include "particles/transformation/CoordinateTransformation.H"
#include "particles/elements/Sbend.H"

#include <AMReX.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>

#include <iostream>
#include <memory>

#include <vector>
#include "particles/wakefields/ChargeBinning.H"
#include "particles/wakefields/WakeConvolution.H"
#include "particles/wakefields/WakePush.H"

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
        if (m_grids_initialized)
        {
            m_lattice.clear();

            // this one last
            amr_data.reset();

            if (amrex::Initialized())
                amrex::Finalize();

            // only finalize once
            m_grids_initialized = false;
        }
    }

    void ImpactX::init_grids ()
    {
        BL_PROFILE("ImpactX::init_grids");

        amr_data = std::make_unique<initialization::AmrCoreData>(initialization::init_amr_core());
        amr_data->m_particle_container = std::make_unique<ImpactXParticleContainer>(amr_data.get());
        amr_data->m_particles_lost = std::make_unique<ImpactXParticleContainer>(amr_data.get());

        // query input for warning logger variables and set up warning logger accordingly
        init_warning_logger();

        // move old diagnostics out of the way
        bool diag_enable = true;
        amrex::ParmParse("diag").queryAdd("enable", diag_enable);
        if (diag_enable) {
            amrex::UtilCreateCleanDirectory("diags", true);
        }

        // the particle container has been set to track the same Geometry as ImpactX

        // this is the earliest point that we need to know the particle shape,
        // so that we can initialize the guard size of our MultiFabs
        amr_data->m_particle_container->SetParticleShape();

        // init blocks / grids & MultiFabs
        amr_data->InitFromScratch(0.0);

        // alloc particle containers
        //   the lost particles have an extra runtime attribute: s when it was lost
        if (!amr_data->m_particles_lost->HasRealComp("s_lost"))
        {
            bool comm = true;
            amr_data->m_particles_lost->AddRealComp("s_lost", comm);
        }

        //   have to resize here, not in the constructor because grids have not
        //   been built when constructor was called.
        amr_data->m_particle_container->reserveData();
        amr_data->m_particle_container->resizeData();
        amr_data->m_particles_lost->reserveData();
        amr_data->m_particles_lost->resizeData();

        // register shortcut
        amr_data->m_particle_container->SetLostParticleContainer(amr_data->m_particles_lost.get());

        // print AMReX grid summary
        if (amrex::ParallelDescriptor::IOProcessor()) {
            // verbosity
            amrex::ParmParse pp_impactx("impactx");
            int verbose = 1;
            pp_impactx.queryAdd("verbose", verbose);

            if (verbose > 0) {
                std::cout << "\nGrids Summary:\n";
                amr_data->printGridSummary(std::cout, 0, amr_data->finestLevel());
            }
        }

        // keep track that init is done
        m_grids_initialized = true;
    }

    void ImpactX::evolve ()
    {

        BL_PROFILE("ImpactX::evolve");

        validate();

        // verbosity
        amrex::ParmParse pp_impactx("impactx");
        int verbose = 1;
        pp_impactx.queryAdd("verbose", verbose);

        // a global step for diagnostics including space charge slice steps in elements
        //   before we start the evolve loop, we are in "step 0" (initial state)
        int global_step = 0;

        // check typos in inputs after step 1
        bool early_params_checked = false;

        amrex::ParmParse pp_diag("diag");
        bool diag_enable = true;
        pp_diag.queryAdd("enable", diag_enable);
        if (verbose > 0) {
            amrex::Print() << " Diagnostics: " << diag_enable << "\n";
        }

        int file_min_digits = 6;
        if (diag_enable)
        {
            pp_diag.queryAdd("file_min_digits", file_min_digits);

            // print initial reference particle to file
            diagnostics::DiagnosticOutput(*amr_data->m_particle_container,
                                          diagnostics::OutputType::PrintRefParticle,
                                          "diags/ref_particle",
                                          global_step);

            // print the initial values of reduced beam characteristics
            diagnostics::DiagnosticOutput(*amr_data->m_particle_container,
                                          diagnostics::OutputType::PrintReducedBeamCharacteristics,
                                          "diags/reduced_beam_characteristics");

        }

        amrex::ParmParse pp_algo("algo");
        bool space_charge = false;
        pp_algo.query("space_charge", space_charge);
        if (verbose > 0) {
            amrex::Print() << " Space Charge effects: " << space_charge << "\n";
        }

        // periods through the lattice
        int periods = 1;
        amrex::ParmParse("lattice").queryAdd("periods", periods);

        for (int cycle=0; cycle < periods; ++cycle) {
            // loop over all beamline elements
            for (auto &element_variant: m_lattice) {
                // update element edge of the reference particle
                amr_data->m_particle_container->SetRefParticleEdge();

                // number of slices used for the application of space charge
                int nslice = 1;
                amrex::ParticleReal slice_ds; // in meters
                std::visit([&nslice, &slice_ds](auto &&element) {
                    nslice = element.nslice();
                    slice_ds = element.ds() / nslice;
                }, element_variant);

                // sub-steps for space charge within the element
                for (int slice_step = 0; slice_step < nslice; ++slice_step) {
                    BL_PROFILE("ImpactX::evolve::slice_step");
                    global_step++;
                    if (verbose > 0) {
                        amrex::Print() << " ++++ Starting global_step=" << global_step
                                       << " slice_step=" << slice_step << "\n";
                    }

                    //CSR Wakefield Response

                    /*bool element_has_csr = false; // Updates to true for example with bend element
                    double R = 0.0; // Updates for bend element rc

                    //Define lambda function inside of std::visit
                    std::visit([&R, &element_has_csr](auto &&element)
                    {
                        if constexpr (std::is_same_v<std::decay_t<decltype(element)>, Sbend>)
                        {
                            R = element.m_rc;
                            std::cout << "My radius of curvature is:" << R << std::endl;
                            element_has_csr = true;
                        }
                        else if constexpr (std::is_same_v<std::decay_t<decltype(element)>, CFbend>)
                        {
                            R = element.m_rc;
                            std::cout << "My radius of curvature is:" << R << std::endl;
                            element_has_csr = true;
                        }
                        
                        else if constexpr (std::is_same_v<std::decay_t<decltype(element)>, ExactSbend>) //Currently internal calculation for m_rc
                        {
                            R = element.m_rc;
                            std::cout << "My radius of curvature is:" << R << std::endl;
                            element_has_csr = true;
                        }
                        
                    }, element_variant);

                    Enter loop if lattice has bend element
                    if (element_has_csr)
                    {
                        // Measure beam size, extract the min, max of particle positions
                        auto const [x_min, y_min, t_min, x_max, y_max, t_max] =
                            amr_data->m_particle_container->MinAndMaxPositions();

                        using amrex::Real;

                        // Set parameters for charge deposition
                        bool is_unity_particle_weight = true;
                        bool GetNumberDensity = true;

                        int num_bins = 100;  // Set resolution
                        Real bin_min = t_min;
                        Real bin_max = t_max;
                        Real bin_size = (bin_max - bin_min) / num_bins;

                        // Allocate memory for the charge profile
                        Real* dptr_data = new Real[num_bins]();
                        auto& particle_container = *(amr_data->m_particle_container);

                        //Call charge deposition function
                        DepositCharge1D(particle_container, dptr_data, num_bins, bin_min, bin_size, is_unity_particle_weight);

                        // Call charge density derivative function
                        std::vector<double> charge_distribution(dptr_data, dptr_data + num_bins);
                        std::vector<double> slopes(num_bins - 1);
                        DerivativeCharge1D(charge_distribution.data(), slopes.data(), num_bins, bin_size, GetNumberDensity); //Use number derivatives for convolution with CSR

                        // Call wake function

                        // Read in external variable bunch_charge
                        std::cout << "My beam charge is:" << bunch_charge << std::endl;

                        std::vector<double> wake_function(num_bins);
                        for (int i = 0; i < num_bins; ++i)
                        {
                            double s = bin_min + i * bin_size;
                            wake_function[i] = w_l_csr(s, R, bunch_charge);
                        }

                        // Call convolution function
                        std::vector<double> convoluted_wakefield(2 * num_bins - 1);
                        convolve_fft(slopes.data(), wake_function.data(), slopes.size(), wake_function.size(), bin_size, convoluted_wakefield.data(), 1);

                        //Check convolution
                        std::cout << "Convoluted wakefield: ";
                        std::ofstream outfile("convoluted_wakefield.txt");
                        for (int i = 0; i < convoluted_wakefield.size(); ++i)
                        {
                            std::cout << convoluted_wakefield[i] << " ";
                            outfile << convoluted_wakefield[i] << std::endl;
                        }
                        std::cout << std::endl;
                        outfile.close();
                        delete[] dptr_data;

                        // Kick particles with wake
                        impactx::wakepush::WakePush(particle_container, convoluted_wakefield, bin_size);
                    }*/

                    // Space-charge calculation: turn off if there is only 1 particle
                    if (space_charge &&
                        amr_data->m_particle_container->TotalNumberOfParticles(true, false)) {

                        // transform from x',y',t to x,y,z
                        transformation::CoordinateTransformation(
                                *amr_data->m_particle_container,
                                CoordSystem::t);

                        // Note: The following operation assume that
                        // the particles are in x, y, z coordinates.

                        // Resize the mesh, based on `m_particle_container` extent
                        ResizeMesh();

                        // Redistribute particles in the new mesh in x, y, z
                        amr_data->m_particle_container->Redistribute();

                        // charge deposition
                        amr_data->m_particle_container->DepositCharge(amr_data->m_rho, amr_data->refRatio());

                        // poisson solve in x,y,z
                        spacecharge::PoissonSolve(*amr_data->m_particle_container, amr_data->m_rho, amr_data->m_phi, amr_data->refRatio());

                        // calculate force in x,y,z
                        spacecharge::ForceFromSelfFields(amr_data->m_space_charge_field,
                                                         amr_data->m_phi,
                                                         amr_data->Geom());

                        // gather and space-charge push in x,y,z , assuming the space-charge
                        // field is the same before/after transformation
                        // TODO: This is currently using linear order.
                        spacecharge::GatherAndPush(*amr_data->m_particle_container,
                                                   amr_data->m_space_charge_field,
                                                   amr_data->Geom(),
                                                   slice_ds);

                        // transform from x,y,z to x',y',t
                        transformation::CoordinateTransformation(*amr_data->m_particle_container,
                                                                 CoordSystem::s);
                    }

                    // for later: original Impact implementation as an option
                    // Redistribute particles in x',y',t
                    //   TODO: only needed if we want to gather and push space charge
                    //         in x',y',t
                    //   TODO: change geometry beforehand according to transformation
                    //m_particle_container->Redistribute();
                    //
                    // in original Impact, we gather and space-charge push in x',y',t ,
                    // assuming that the distribution did not change

                    // push all particles with external maps
                    Push(*amr_data->m_particle_container, element_variant, global_step);

                    // move "lost" particles to another particle container
                    collect_lost_particles(*amr_data->m_particle_container);

                    // just prints an empty newline at the end of the slice_step
                    if (verbose > 0) {
                        amrex::Print() << "\n";
                    }

                    // slice-step diagnostics
                    bool slice_step_diagnostics = false;
                    pp_diag.queryAdd("slice_step_diagnostics", slice_step_diagnostics);

                    if (diag_enable && slice_step_diagnostics) {
                        // print slice step reference particle to file
                        diagnostics::DiagnosticOutput(*amr_data->m_particle_container,
                                                      diagnostics::OutputType::PrintRefParticle,
                                                      "diags/ref_particle",
                                                      global_step,
                                                      true);

                        // print slice step reduced beam characteristics to file
                        diagnostics::DiagnosticOutput(*amr_data->m_particle_container,
                                                      diagnostics::OutputType::PrintReducedBeamCharacteristics,
                                                      "diags/reduced_beam_characteristics",
                                                      global_step,
                                                      true);

                    }

                    // inputs: unused parameters (e.g. typos) check after step 1 has finished
                    if (!early_params_checked) { early_params_checked = early_param_check(); }

                } // end in-element space-charge slice-step loop

            } // end beamline element loop
        } // end periods though the lattice loop

        if (diag_enable)
        {
            // print final reference particle to file
            diagnostics::DiagnosticOutput(*amr_data->m_particle_container,
                                          diagnostics::OutputType::PrintRefParticle,
                                          "diags/ref_particle_final",
                                          global_step);

            // print the final values of the reduced beam characteristics
            diagnostics::DiagnosticOutput(*amr_data->m_particle_container,
                                          diagnostics::OutputType::PrintReducedBeamCharacteristics,
                                          "diags/reduced_beam_characteristics_final",
                                          global_step);

            // output particles lost in apertures
            if (amr_data->m_particles_lost->TotalNumberOfParticles() > 0)
            {
                std::string openpmd_backend = "default";
                pp_diag.queryAdd("backend", openpmd_backend);

                diagnostics::BeamMonitor output_lost("particles_lost", openpmd_backend, "g");
                output_lost(*amr_data->m_particles_lost, 0);
                output_lost.finalize();
            }
        }

        // loop over all beamline elements & finalize them
        for (auto & element_variant : m_lattice)
        {
            std::visit([](auto&& element){
                element.finalize();
            }, element_variant);
        }
    }
} // namespace impactx
