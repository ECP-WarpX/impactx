/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#include "ImpactX.H"
#include "diagnostics/DiagnosticOutput.H"
#include "initialization/Algorithms.H"
#include "initialization/InitAmrCore.H"
#include "particles/ImpactXParticleContainer.H"
#include "particles/Push.H"
#include "tracking/common.H"

#include <AMReX.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <memory>
#include <stdexcept>


namespace impactx
{
    void
    ImpactX::track_reference (RefPart & ref)
    {
        BL_PROFILE("ImpactX::track_reference");

        if (m_lattice->empty())
        {
            throw std::runtime_error(
                "Beamline lattice has zero elements. Not yet initialized?");
        }

        // verbosity
        amrex::ParmParse pp_impactx("impactx");
        int verbose = 1;
        pp_impactx.queryAddWithParser("verbose", verbose);

        // a global step for diagnostics including space charge slice steps in elements
        //   before we start the tracking loop, we are in "step 0" (initial state)
        int & step = m_tracking_state.m_step;
        step = 0;
        m_tracking_state.m_direction = TrackingDirection::Forward;

        // check typos in inputs after step 1
        bool early_params_checked = false;

        // output of init state
        amrex::ParmParse pp_diag("diag");
        bool diag_enable = true;
        pp_diag.queryAdd("enable", diag_enable);
        if (verbose > 0)
        {
            amrex::Print() << " Diagnostics: " << diag_enable << "\n";
        }

        if (diag_enable)
        {
            int file_min_digits = 6;
            pp_diag.queryAddWithParser("file_min_digits", file_min_digits);

            // print initial reference particle to file
            diagnostics::DiagnosticOutput(ref, "ref_particle");

        }

        auto space_charge = get_space_charge_algo();
        if (space_charge != SpaceChargeAlgo::False)
        {
            throw std::runtime_error("Space charge effects cannot be modeled for single particle tracking.");
        }

        amrex::ParmParse const pp_algo("algo");
        bool csr = false;
        pp_algo.query("csr", csr);
        if (csr)
        {
            throw std::runtime_error(
                "Coherent Synchrotron Radiation (CSR) cannot be "
                "modeled for single particle tracking. "
                "Please set the CSR option to false."
            );
        }

        // the external-field transport map ``M`` of one element slice
        auto element_push = [&ref] (
            elements::ElementHandle & element_variant,
            //! (unused) the reference particle has no per-step or per-period element output
            [[maybe_unused]] int step_,
            [[maybe_unused]] int period_
        )
        {
            // push the reference particle with external maps
            push(ref, element_variant);
        };

        // book-keeping and diagnostics, applied once at the end of each slice
        auto slice_diagnostics = [
            this, &ref, verbose, &pp_diag, diag_enable, &early_params_checked
        ] (
            int step_,
            //! (unused) the reference particle output carries no period
            [[maybe_unused]] int period_
        )
        {
            // just prints an empty newline at the end of the slice_step
            if (verbose > 0)
            {
                amrex::Print() << "\n";
            }

            // slice-step diagnostics
            bool slice_step_diagnostics = false;
            pp_diag.queryAdd("slice_step_diagnostics", slice_step_diagnostics);

            if (diag_enable && slice_step_diagnostics)
            {
                // print slice step reference particle to file
                diagnostics::DiagnosticOutput(ref, "ref_particle", step_, true);
            }

            // inputs: unused parameters (e.g. typos) check after step 1 has finished
            if (!early_params_checked) { early_params_checked = early_param_check(); }
        };

        // traverse the lattice, applying the element transport per element slice
        //   Reference-particle tracking rejects space charge and CSR above and models no
        //   other collective effect, so it applies ``M`` alone and never calls the kick
        //   (\see track_lattice).
        track_lattice(
            *m_lattice,
            ref,
            m_tracking_state,
            false, // no collective effects
            false, // nothing to Strang-split
            [this](std::string const & name) { call_hook(name); },
            [](elements::ElementHandle &, amrex::ParticleReal) {}, // never called
            element_push,
            slice_diagnostics
        );

        if (diag_enable)
        {
            // print final reference particle to file
            diagnostics::DiagnosticOutput(ref, "ref_particle_final", step);
        }
    }
} // namespace impactx
