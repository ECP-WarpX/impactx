/* Copyright 2022 The Regents of the University of California, through Lawrence
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
#include <ablastr/utils/TextMsg.H>

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
    }

    void ImpactX::initGrids ()
    {
        BL_PROFILE("ImpactX::initGrids");

        // this is the earliest point that we need to know the particle shape,
        // so that we can initialize the guard size of our MultiFabs
        m_particle_container->SetParticleShape();

        // init blocks / grids & MultiFabs
        AmrCore::InitFromScratch(0.0);
        amrex::Print() << "boxArray(0) " << boxArray(0) << std::endl;

        // move old diagnostics out of the way
        amrex::UtilCreateCleanDirectory("diags", true);
    }

    void ImpactX::evolve ()
    {
        BL_PROFILE("ImpactX::evolve");

        validate();

        // a global step for diagnostics including space charge slice steps in elements
        //   before we start the evolve loop, we are in "step 0" (initial state)
        int global_step = 0;

        bool early_params_checked = false; // check typos in inputs after step 1

        // count particles - if no particles are found in our particle container, then a lot of
        // AMReX routines over ParIter won't work and we have nothing to do here anyways
        {
            int const nLevelPC = finestLevel();
            amrex::Long nParticles = 0;
            for (int lev = 0; lev <= nLevelPC; ++lev) {
                nParticles += m_particle_container->NumberOfParticlesAtLevel(lev);
            }
            if (nParticles == 0) {
                amrex::Abort("No particles found. Cannot run evolve without a beam.");
                return;
            }
        }

        amrex::ParmParse pp_diag("diag");
        bool diag_enable = true;
        pp_diag.queryAdd("enable", diag_enable);
        amrex::Print() << " Diagnostics: " << diag_enable << "\n";

        int file_min_digits = 6;
        if (diag_enable)
        {
            pp_diag.queryAdd("file_min_digits", file_min_digits);

            // print initial particle distribution to file
            std::string diag_name = amrex::Concatenate("diags/beam_", global_step, file_min_digits);
            diagnostics::DiagnosticOutput(*m_particle_container,
                                          diagnostics::OutputType::PrintParticles,
                                          diag_name);

            // print initial reference particle to file
            diagnostics::DiagnosticOutput(*m_particle_container,
                                          diagnostics::OutputType::PrintRefParticle,
                                          "diags/ref_particle",
                                          global_step);

            // print the initial values of the two invariants H and I
            diag_name = amrex::Concatenate("diags/nonlinear_lens_invariants_", global_step, file_min_digits);
            diagnostics::DiagnosticOutput(*m_particle_container,
                                          diagnostics::OutputType::PrintNonlinearLensInvariants,
                                          diag_name);

        }

        amrex::ParmParse pp_algo("algo");
        bool space_charge = true;
        pp_algo.queryAdd("space_charge", space_charge);
        amrex::Print() << " Space Charge effects: " << space_charge << "\n";

        amrex::ParmParse pp_impactx("impactx");

        //"Synthetic" warning messages may be injected in the Warning Manager via
        // inputfile for debug&testing purposes.
        ablastr::warn_manager::GetWMInstance().debug_read_warnings_from_input(pp_impactx);

        // Set the flag to control if WarpX has to emit a warning message as soon as a warning is recorded
        bool always_warn_immediately = false;
        pp_impactx.query("always_warn_immediately", always_warn_immediately);
        ablastr::warn_manager::GetWMInstance().SetAlwaysWarnImmediately(always_warn_immediately);

        // Set the WarnPriority threshold to decide if ImpactX has to abort when a warning is recorded
        if(std::string str_abort_on_warning_threshold = "";
                pp_impactx.query("abort_on_warning_threshold", str_abort_on_warning_threshold)){
            std::optional<ablastr::warn_manager::WarnPriority> abort_on_warning_threshold = std::nullopt;
            if (str_abort_on_warning_threshold == "high")
                abort_on_warning_threshold = ablastr::warn_manager::WarnPriority::high;
            else if (str_abort_on_warning_threshold == "medium" )
                abort_on_warning_threshold = ablastr::warn_manager::WarnPriority::medium;
            else if (str_abort_on_warning_threshold == "low")
                abort_on_warning_threshold = ablastr::warn_manager::WarnPriority::low;
            else {
                amrex::Abort(ablastr::utils::TextMsg::Err(str_abort_on_warning_threshold
                                          +"is not a valid option for impactx.abort_on_warning_threshold (use: low, medium or high)"));
            }
            ablastr::warn_manager::GetWMInstance().SetAbortThreshold(abort_on_warning_threshold);
        }

        // loop over all beamline elements
        for (auto & element_variant : m_lattice)
        {
            // number of slices used for the application of space charge
            int nslice = 1;
            amrex::ParticleReal slice_ds; // in meters
            std::visit([&nslice, &slice_ds](auto&& element){
                nslice = element.nslice();
                slice_ds = element.ds() / nslice;
            }, element_variant);

            // sub-steps for space charge within the element
            for (int slice_step = 0; slice_step < nslice; ++slice_step)
            {
                BL_PROFILE("ImpactX::evolve::slice_step");
                global_step++;
                amrex::Print() << " ++++ Starting global_step=" << global_step
                               << " slice_step=" << slice_step << "\n";

                // Space-charge calculation: turn off if there is only 1 particle
                if (space_charge &&
                    m_particle_container->TotalNumberOfParticles(false,false) > 1)
                {

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
                Push(*m_particle_container, element_variant);

                // just prints an empty newline at the end of the slice_step
                amrex::Print() << "\n";

                // slice-step diagnostics
                bool slice_step_diagnostics = false;
                pp_diag.queryAdd("slice_step_diagnostics", slice_step_diagnostics);

                if (diag_enable && slice_step_diagnostics)
                {
                    // print slice step particle distribution to file
                    std::string diag_name = amrex::Concatenate("diags/beam_", global_step, file_min_digits);
                    diagnostics::DiagnosticOutput(*m_particle_container,
                                                  diagnostics::OutputType::PrintParticles,
                                                  diag_name,
                                                  global_step);

                    // print slice step reference particle to file
                    diagnostics::DiagnosticOutput(*m_particle_container,
                                                  diagnostics::OutputType::PrintRefParticle,
                                                  "diags/ref_particle",
                                                  global_step,
                                                  true);

                }

                // inputs: unused parameters (e.g. typos) check after step 1 has finished
                if (!early_params_checked) {
                    amrex::Print() << "\n"; // better: conditional \n based on return value
                    amrex::ParmParse().QueryUnusedInputs();

                    //Print the warning list right after the first step.
                    amrex::Print() << ablastr::warn_manager::GetWMInstance()
                                      .PrintGlobalWarnings("FIRST STEP");
                    early_params_checked = true;
                }

            } // end in-element space-charge slice-step loop
        } // end beamline element loop

        if (diag_enable)
        {
            // print final particle distribution to file
            diagnostics::DiagnosticOutput(*m_particle_container,
                                          diagnostics::OutputType::PrintParticles,
                                          "diags/beam_final",
                                          global_step);

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
        }

    }
} // namespace impactx
