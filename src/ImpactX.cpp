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
#include "initialization/InitOneBoxPerRank.H"
#include "particles/ImpactXParticleContainer.H"
#include "particles/Push.H"
#include "particles/transformation/CoordinateTransformation.H"
#include "particles/diagnostics/DiagnosticOutput.H"

#include <AMReX.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>

#include <memory>


namespace impactx
{
    ImpactX::ImpactX ()
        : AmrCore(initialization::one_box_per_rank()),
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

        // a global step for diagnostics including space charge slice steps in elements
        //   before we start the evolve loop, we are in "step 0" (initial state)
        int global_step = 0;

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

        // loop over all beamline elements
        for (auto & element_variant : m_lattice)
        {
            // number of slices used for the application of space charge
            int nslice = 1;
            std::visit([&nslice](auto&& element){ nslice = element.nslice(); }, element_variant);

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
                    transformation::CoordinateTransformation(*m_particle_container,
                                                             transformation::Direction::T2Z);

                    // Note: The following operation assume that
                    // the particles are in x, y, z coordinates.

                    // Resize the mesh, based on `m_particle_container` extent
                    ResizeMesh();

                    // Redistribute particles in the new mesh in x, y, z
                    m_particle_container->Redistribute();

                    // charge deposition
                    m_particle_container->DepositCharge(m_rho, this->refRatio());

                    // poisson solve in x,y,z
                    //   TODO

                    // gather and space-charge push in x,y,z , assuming the space-charge
                    // field is the same before/after transformation
                    //   TODO

                    // transform from x,y,z to x',y',t
                    transformation::CoordinateTransformation(*m_particle_container,
                                                             transformation::Direction::Z2T);
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
