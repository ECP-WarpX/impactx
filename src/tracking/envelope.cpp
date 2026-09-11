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
#include "envelope/spacecharge/EnvelopeSpaceChargePush.H"
#include "initialization/Algorithms.H"
#include "initialization/InitAmrCore.H"
#include "particles/ImpactXParticleContainer.H"
#include "particles/Push.H"
#include "tracking/common.H"

#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>

#include <memory>
#include <stdexcept>


namespace impactx
{
    void
    ImpactX::track_envelope ()
    {
        BL_PROFILE("ImpactX::track_envelope");

        if (m_lattice->empty())
        {
            throw std::runtime_error(
                "Beamline lattice has zero elements. Not yet initialized?");
        }

        using namespace amrex::literals;

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

        // access beam data
        if (!amr_data->track_envelope.m_ref.has_value())
        {
            throw std::runtime_error("track_envelope: Reference particle not set.");
        }
        if (!amr_data->track_envelope.m_env.has_value())
        {
            throw std::runtime_error("track_envelope: Envelope (covariance matrix) not set.");
        }
        auto & ref = amr_data->track_envelope.m_ref.value();
        auto & env = amr_data->track_envelope.m_env.value();
        auto & cm = env.m_env;
        auto & intensity = env.m_beam_intensity;

        // output of init state
        amrex::ParmParse pp_diag("diag");
        bool diag_enable = true;
        pp_diag.queryAdd("enable", diag_enable);
        if (verbose > 0) {
            amrex::Print() << " Diagnostics: " << diag_enable << "\n";
        }

        if (diag_enable)
        {
            int file_min_digits = 6;
            pp_diag.queryAddWithParser("file_min_digits", file_min_digits);

            // print initial reference particle to file
            diagnostics::DiagnosticOutput(ref, "ref_particle");

            // print the initial values of reduced beam characteristics
            diagnostics::DiagnosticOutput(cm, ref, "reduced_beam_characteristics");

        }

        amrex::ParmParse const pp_algo("algo");
        auto space_charge = get_space_charge_algo();
        if (verbose > 0)
        {
            amrex::Print() << " Space Charge effects: " << to_string(space_charge) << "\n";
        }
        if (space_charge == SpaceChargeAlgo::True_3D && intensity == 0_prt) {
            ablastr::warn_manager::WMRecordWarning(
                "algo.space_charge",
                "Space charge calculations are enabled but zero bunch charge was provided. "
                "Skipping space charge calculations.",
                ablastr::warn_manager::WarnPriority::high
            );
        }
        if (space_charge == SpaceChargeAlgo::True_2D && intensity == 0_prt) {
            ablastr::warn_manager::WMRecordWarning(
                "algo.space_charge",
                "Space charge calculations are enabled but zero beam current was provided. "
                "Skipping space charge calculations.",
                ablastr::warn_manager::WarnPriority::high
            );
        }
        // envelope tracking models space charge as the linear field of the rms-equivalent
        //   uniform beam ellipse (2D) or ellipsoid (3D); the other models are particle-only
        if (space_charge == SpaceChargeAlgo::Gauss3D ||
            space_charge == SpaceChargeAlgo::Gauss2p5D ||
            space_charge == SpaceChargeAlgo::True_2p5D)
        {
            throw std::runtime_error(
                to_string(space_charge) + " space charge force calculation is only supported "
                "with particle tracking. For envelope tracking, use: 2D or 3D"
            );
        }

        bool csr = false;
        pp_algo.query("csr", csr);
        if (verbose > 0)
        {
            amrex::Print() << " CSR effects: " << csr << "\n";
        }
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!csr, "CSR effects are not yet implemented for envelope tracking.");

        // whether any collective effect is active: only then is a kick applied per slice
        //   At zero intensity the kick leaves the covariance matrix unchanged, so it is
        //   skipped entirely: that is what the warning above announces, and it keeps the
        //   element transport of such a run unsplit.
        bool const collective_effects =
            (space_charge == SpaceChargeAlgo::True_2D ||
             space_charge == SpaceChargeAlgo::True_3D) &&
            intensity != 0_prt;

        // second-order Strang split of the collective kicks, on by default
        //   Disabling it composes kick and transport to first order instead, which is what
        //   most other codes do: useful to compare against them and to show convergence.
        //   This is the same option as in particle tracking, so that both models compose
        //   their collective kicks to the same order and stay comparable.
        bool strang_split = true;
        pp_algo.query("strang_split", strang_split);
        if (verbose > 0 && collective_effects)
        {
            amrex::Print() << " Strang split: " << strang_split << "\n";
        }

        // the collective effect kick ``K`` of one element slice
        auto collective_kicks = [&ref, &cm, &intensity, space_charge] (
            //! (unused) the kick does not depend on the element
            [[maybe_unused]] elements::ElementHandle & element_variant,
            amrex::ParticleReal kick_ds
        )
        {
            if (space_charge == SpaceChargeAlgo::True_2D)
            {
                // push Covariance Matrix in 2D space charge fields
                envelope::spacecharge::space_charge2D_push(ref, cm, intensity, kick_ds);
            } else if (space_charge == SpaceChargeAlgo::True_3D)
            {
                // push Covariance Matrix in 3D space charge fields
                envelope::spacecharge::space_charge3D_push(ref, cm, intensity, kick_ds);
            }
        };

        // the external-field transport map ``M`` of one element slice
        //   The Strang split around collective effects applies this twice per slice, once
        //   per half-map, so it carries no book-keeping: that lives in @see
        //   slice_diagnostics below.
        auto element_push = [&ref, &cm] (
            elements::ElementHandle & element_variant,
            [[maybe_unused]] int step_,   //! (unused) the envelope has no per-step element output
            [[maybe_unused]] int period_  //! (unused) the envelope has no per-period element output
        )
        {
            elements::visit([&ref, &cm](auto&& element)
            {
                // push reference particle in global coordinates
                {
                    BL_PROFILE("impactx::push::RefPart");
                    element(ref);
                }

                // push Covariance Matrix in external fields
                element(cm, ref);

            }, element_variant);
        };

        // book-keeping and diagnostics, applied once at the end of each slice
        auto slice_diagnostics = [
            this, &ref, &cm, verbose, &pp_diag, diag_enable, &early_params_checked
        ] (
            int step_,
            int period_
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

                // print slice step reduced beam characteristics to file
                diagnostics::DiagnosticOutput(
                    cm, ref, "reduced_beam_characteristics", step_, period_, true);
            }

            // inputs: unused parameters (e.g. typos) check after step 1 has finished
            if (!early_params_checked) { early_params_checked = early_param_check(); }
        };

        // traverse the lattice, applying the collective kick and the
        // element transport per element slice (\see track_lattice)
        track_lattice(
            *m_lattice,
            ref,
            m_tracking_state,
            collective_effects,
            strang_split,
            [this](std::string const & name) { call_hook(name); },
            collective_kicks,
            element_push,
            slice_diagnostics
        );

        if (diag_enable)
        {
            // print final reference particle to file
            diagnostics::DiagnosticOutput(ref, "ref_particle_final", step);

            // print the final values of the reduced beam characteristics
            diagnostics::DiagnosticOutput(
                cm, ref, "reduced_beam_characteristics_final", step, m_tracking_state.m_period);
        }
    }
} // namespace impactx
