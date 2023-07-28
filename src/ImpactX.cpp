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
#include "particles/ImpactXParticleContainer.H"
#include "particles/Push.H"
#include "particles/diagnostics/DiagnosticOutput.H"
#include "particles/spacecharge/ForceFromSelfFields.H"
#include "particles/spacecharge/GatherAndPush.H"
#include "particles/spacecharge/PoissonSolve.H"
#include "particles/transformation/CoordinateTransformation.H"

#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>

#include <memory>


namespace impactx
{
    ImpactX::ImpactX ()
        : AmrCore(initialization::init_amr_core()),
          m_particle_container(std::make_unique<ImpactXParticleContainer>(this))
    {
        // todo: if amr.n_cells is provided, overwrite/redefine AmrCore object

        // todo: if charge deposition and/or space charge are requested, require
        //       amr.n_cells from user inputs

        // query input for warning logger variables and set up warning logger accordingly
        init_warning_logger();

        // move old diagnostics out of the way
        amrex::UtilCreateCleanDirectory("diags", true);
    }

    void ImpactX::initGrids ()
    {
        BL_PROFILE("ImpactX::initGrids");

        // n_cells has been set using temporary values earlier. We now know the true value of
        // n_cells, so we recompute the Geometry objects for each level here *if* the user
        // has set n_cells in the inputs file
        {
            amrex::Vector<int> n_cell(AMREX_SPACEDIM);
            amrex::ParmParse pp_amr("amr");
            pp_amr.queryarr("n_cell", n_cell);

            amrex::IntVect lo(amrex::IntVect::TheZeroVector()), hi(n_cell);
            hi -= amrex::IntVect::TheUnitVector();
            amrex::Box index_domain(lo,hi);
            for (int i = 0; i <= max_level; i++)
            {
                geom[i].Domain(index_domain);
                if (i < max_level) {
                    index_domain.refine(ref_ratio[i]);
                }
            }
        }

        // the particle container has been set to track the same Geometry as ImpactX

        // this is the earliest point that we need to know the particle shape,
        // so that we can initialize the guard size of our MultiFabs
        m_particle_container->SetParticleShape();

        // init blocks / grids & MultiFabs
        AmrCore::InitFromScratch(0.0);
        amrex::Print() << "boxArray(0) " << boxArray(0) << std::endl;
    }

    void ImpactX::evolve ()
    {
        BL_PROFILE("ImpactX::evolve");

        validate();

        // a global step for diagnostics including space charge slice steps in elements
        //   before we start the evolve loop, we are in "step 0" (initial state)
        int global_step = 0;

        // check typos in inputs after step 1
        bool early_params_checked = false;

        amrex::ParmParse pp_diag("diag");
        bool diag_enable = true;
        pp_diag.queryAdd("enable", diag_enable);
        amrex::Print() << " Diagnostics: " << diag_enable << "\n";

        int file_min_digits = 6;
        if (diag_enable)
        {
            pp_diag.queryAdd("file_min_digits", file_min_digits);

            // print initial reference particle to file
            diagnostics::DiagnosticOutput(*m_particle_container,
                                          diagnostics::OutputType::PrintRefParticle,
                                          "diags/ref_particle",
                                          global_step);

            // print the initial values of the two invariants H and I
            std::string diag_name = amrex::Concatenate("diags/nonlinear_lens_invariants_", global_step, file_min_digits);
            diagnostics::DiagnosticOutput(*m_particle_container,
                                          diagnostics::OutputType::PrintNonlinearLensInvariants,
                                          diag_name);

            // print the initial values of reduced beam characteristics
            diagnostics::DiagnosticOutput(*m_particle_container,
                                          diagnostics::OutputType::PrintReducedBeamCharacteristics,
                                          "diags/reduced_beam_characteristics");

        }

        amrex::ParmParse pp_algo("algo");
        bool space_charge = true;
        pp_algo.queryAdd("space_charge", space_charge);
        amrex::Print() << " Space Charge effects: " << space_charge << "\n";

        // periods through the lattice
        int periods = 1;
        amrex::ParmParse("lattice").queryAdd("periods", periods);

        for (int cycle=0; cycle < periods; ++cycle) {
            // loop over all beamline elements
            for (auto &element_variant: m_lattice) {
                // update element edge of the reference particle
                m_particle_container->SetRefParticleEdge();

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
                    amrex::Print() << " ++++ Starting global_step=" << global_step
                                   << " slice_step=" << slice_step << "\n";

                    // Space-charge calculation: turn off if there is only 1 particle
                    if (space_charge &&
                        m_particle_container->TotalNumberOfParticles(false, false) > 1) {

                        // transform from x',y',t to x,y,z
                        transformation::CoordinateTransformation(
                                *m_particle_container,
                                transformation::Direction::to_fixed_t);

                        // Note: The following operation assume that
                        // the particles are in x, y, z coordinates.

                        // Resize the mesh, based on `m_particle_container` extent
                        ResizeMesh();

                        // Redistribute particles in the new mesh in x, y, z
                        m_particle_container->Redistribute();

                        // charge deposition
                        m_particle_container->DepositCharge(m_rho, this->refRatio());

                        // poisson solve in x,y,z
                        spacecharge::PoissonSolve(*m_particle_container, m_rho, m_phi);

                        // calculate force in x,y,z
                        spacecharge::ForceFromSelfFields(m_space_charge_field,
                                                         m_phi,
                                                         this->geom);

                        // gather and space-charge push in x,y,z , assuming the space-charge
                        // field is the same before/after transformation
                        // TODO: This is currently using linear order.
                        spacecharge::GatherAndPush(*m_particle_container,
                                                   m_space_charge_field,
                                                   this->geom,
                                                   slice_ds);

                        // transform from x,y,z to x',y',t
                        transformation::CoordinateTransformation(*m_particle_container,
                                                                 transformation::Direction::to_fixed_s);
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
                    Push(*m_particle_container, element_variant, global_step);

                    // just prints an empty newline at the end of the slice_step
                    amrex::Print() << "\n";

                    // slice-step diagnostics
                    bool slice_step_diagnostics = false;
                    pp_diag.queryAdd("slice_step_diagnostics", slice_step_diagnostics);

                    if (diag_enable && slice_step_diagnostics) {
                        // print slice step reference particle to file
                        diagnostics::DiagnosticOutput(*m_particle_container,
                                                      diagnostics::OutputType::PrintRefParticle,
                                                      "diags/ref_particle",
                                                      global_step,
                                                      true);

                        // print slice step reduced beam characteristics to file
                        diagnostics::DiagnosticOutput(*m_particle_container,
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
            diagnostics::DiagnosticOutput(*m_particle_container,
                                          diagnostics::OutputType::PrintRefParticle,
                                          "diags/ref_particle_final",
                                          global_step);

            // print the final values of the two invariants H and I
            diagnostics::DiagnosticOutput(*m_particle_container,
                                          diagnostics::OutputType::PrintNonlinearLensInvariants,
                                          "diags/nonlinear_lens_invariants_final",
                                          global_step);

            // print the final values of the reduced beam characteristics
            diagnostics::DiagnosticOutput(*m_particle_container,
                                          diagnostics::OutputType::PrintReducedBeamCharacteristics,
                                          "diags/reduced_beam_characteristics_final",
                                          global_step);
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
