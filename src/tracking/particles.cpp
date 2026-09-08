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
#include "particles/CollectLost.H"
#include "particles/ImpactXParticleContainer.H"
#include "particles/ParticleBoundary.H"
#include "particles/Push.H"
#include "particles/spacecharge/HandleSpacecharge.H"
#include "particles/wakefields/HandleISR.H"
#include "particles/wakefields/HandleWakefield.H"
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
    void ImpactX::track_particles ()
    {
        BL_PROFILE("ImpactX::track_particles");

        validate();

        // verbosity
        amrex::ParmParse pp_impactx("impactx");
        int verbose = 1;
        pp_impactx.queryAddWithParser("verbose", verbose);

        // book-keeping:
        //   a global step for diagnostics including space charge slice steps in elements
        //   before we start the tracking loop, we are in "step 0" (initial state)
        int & step = m_tracking_state.m_step;
        step = 0;
        m_tracking_state.m_direction = TrackingDirection::Forward;

        // check typos in inputs after step 1
        bool early_params_checked = false;

        // shortcuts
        auto & pc = amr_data->track_particles.m_particle_container;

        // diagnostics
        amrex::ParmParse pp_diag("diag");
        bool diag_enable = true;
        pp_diag.queryAdd("enable", diag_enable);
        if (verbose > 0) {
            amrex::Print() << " Diagnostics: " << diag_enable << "\n";
        }

        pc->reset_beam_moments_history();

        if (diag_enable)
        {
            int file_min_digits = 6;
            pp_diag.queryAddWithParser("file_min_digits", file_min_digits);

            // print initial reference particle to file
            diagnostics::DiagnosticOutput(amr_data->track_particles.m_particle_container->GetRefParticle(),
                                          "ref_particle",
                                          step);

            // print the initial values of reduced beam characteristics
            diagnostics::DiagnosticOutput(*amr_data->track_particles.m_particle_container,
                                          "reduced_beam_characteristics");

        }

        auto space_charge = get_space_charge_algo();
        if (verbose > 0) {
            amrex::Print() << " Space Charge effects: " << to_string(space_charge) << "\n";
        }

        amrex::ParmParse const pp_algo("algo");
        bool csr = false;
        pp_algo.query("csr", csr);
        bool isr = false;
        pp_algo.query("isr", isr);
        bool short_range_wakefields = false;
        pp_algo.query("short_range_wakefields", short_range_wakefields);

        // whether any collective effect is active: only then is a kick applied per slice
        bool const collective_effects =
            (space_charge != SpaceChargeAlgo::False) || csr || isr || short_range_wakefields;

        // second-order Strang split of the collective kicks, on by default
        //   Disabling it composes kick and transport to first order instead, which is what
        //   most other codes do: useful to compare against them and to show convergence.
        bool strang_split = true;
        pp_algo.query("strang_split", strang_split);

        bool spin = false;
        pp_algo.query("spin", spin);

        if (spin && pc->GetRefParticle().gyromagnetic_anomaly == 0.0) {
            throw std::runtime_error(
                "algo.spin: Spin tracking is enabled, but the gyromagnetic "
                "anomaly of the reference particle is zero. Either disable spin "
                "tracking, set the reference particle species or "
                "set the value of the gyromagnetic anomaly on it."
            );
        }

        if (verbose > 0) {
            amrex::Print() << " CSR effects: " << csr << "\n";
            amrex::Print() << " ISR effects: " << isr << "\n";
            amrex::Print() << " Short-range wakefields: " << short_range_wakefields << "\n";
            amrex::Print() << " Spin tracking: " << spin << "\n";
            if (collective_effects) {
                amrex::Print() << " Strang split: " << strang_split << "\n";
            }
        }

        // collective effect kicks applied per element slice:
        // space charge, coherent and incoherent synchrotron radiation, etc.
        auto collective_kicks = [this, &pc] (
            elements::KnownElements & element_variant,
            amrex::ParticleReal slice_ds
        )
        {
            // Wakefield calculation: call wakefield function to apply wake effects
            particles::wakefields::HandleWakefield(*pc, element_variant, slice_ds);

            // ISR calculation: call ISR function to apply incoherent synchrotron radiation effects
            particles::wakefields::HandleISR(*pc, element_variant, slice_ds);

            // Space-charge calculation
            particles::spacecharge::HandleSpacecharge(amr_data, [this](){ this->ResizeMesh(); }, slice_ds);
        };

        // the external-field transport map ``M`` of one element slice
        //   The Strang split around collective effects applies this twice per slice, once
        //   per half-map, so it carries no book-keeping: that lives in @see
        //   slice_diagnostics below.
        auto element_push = [&pc] (
            elements::KnownElements & element_variant,
            int step_,
            int period_
        )
        {
            // push all particles with external maps
            push(*pc, element_variant, step_, period_);
        };

        // book-keeping and diagnostics, applied once at the end of each slice
        auto slice_diagnostics = [
            this, &pc, verbose, &pp_diag, diag_enable, &early_params_checked
        ] (
            int step_,
            int period_
        )
        {
            // Apply optional particle boundary conditions
            particles::ParticleBoundary(*pc);

            // move "lost" particles to another particle container
            collect_lost_particles(*pc);

            // just prints an empty newline at the end of the slice_step
            if (verbose > 0) {
                amrex::Print() << "\n";
            }

            // slice-step diagnostics
            bool slice_step_diagnostics = false;
            pp_diag.queryAdd("slice_step_diagnostics", slice_step_diagnostics);

            if (pc->store_beam_moments) {
                pc->record_beam_moments();
            }

            if (diag_enable && slice_step_diagnostics) {
                // print slice step reference particle to file
                diagnostics::DiagnosticOutput(pc->GetRefParticle(),
                                              "ref_particle",
                                              step_,
                                              true);

                // print slice step reduced beam characteristics to file
                diagnostics::DiagnosticOutput(*pc,
                                              "reduced_beam_characteristics",
                                              step_,
                                              period_,
                                              true);
            }

            // inputs: unused parameters (e.g. typos) check after step 1 has finished
            if (!early_params_checked) { early_params_checked = early_param_check(); }
        };

        // traverse the lattice, applying the collective kick and the
        // element transport per element slice (\see track_lattice)
        track_lattice(
            m_lattice,
            pc->GetRefParticle(),
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
            diagnostics::DiagnosticOutput(amr_data->track_particles.m_particle_container->GetRefParticle(),
                                          "ref_particle_final",
                                          step);

            // print the final values of the reduced beam characteristics
            diagnostics::DiagnosticOutput(*amr_data->track_particles.m_particle_container,
                                          "reduced_beam_characteristics_final",
                                          step,
                                          m_tracking_state.m_period);

            // output particles lost in apertures
            if (amr_data->track_particles.m_particles_lost->TotalNumberOfParticles() > 0)
            {
                std::string openpmd_backend = "default";
                pp_diag.queryAdd("backend", openpmd_backend);

                elements::diagnostics::BeamMonitor output_lost("particles_lost", openpmd_backend, "g");
                output_lost(*amr_data->track_particles.m_particles_lost, 0, 0);
                output_lost.finalize();
            }
        }
    }
} // namespace impactx
