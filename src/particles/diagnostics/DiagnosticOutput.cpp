/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#include "DiagnosticOutput.H"
#include "NonlinearLensInvariants.H"
#include "ReducedBeamCharacteristics.H"

#include <AMReX_BLProfiler.H> // for BL_PROFILE
#include <AMReX_Extension.H>  // for AMREX_RESTRICT
#include <AMReX_ParmParse.H>  // for ParmParse
#include <AMReX_REAL.H>       // for ParticleReal
#include <AMReX_Print.H>      // for PrintToFile
#include <AMReX_ParticleTile.H>     // for constructor of SoAParticle

#include <limits>
#include <utility>


namespace impactx::diagnostics
{
    void DiagnosticOutput (ImpactXParticleContainer const & pc,
                           OutputType const otype,
                           std::string file_name,
                           int step,
                           bool append)
    {
        BL_PROFILE("impactx::diagnostics::DiagnosticOutput");

        using namespace amrex::literals; // for _rt and _prt

        // keep file open as we add more and more lines
        amrex::AllPrintToFile file_handler(std::move(file_name));
        file_handler.SetPrecision(std::numeric_limits<amrex::ParticleReal>::max_digits10);

        // write file header per MPI RANK
        if (!append) {
            if (otype == OutputType::PrintRefParticle) {
                file_handler << "step s beta gamma beta_gamma x y z t px py pz pt\n";
            } else if (otype == OutputType::PrintReducedBeamCharacteristics) {
                file_handler << "step" << " " << "s" << " "
                             << "x_mean" << " " << "x_min" << " " << "x_max" << " "
                             << "y_mean" << " " << "y_min" << " " << "y_max" << " "
                             << "t_mean" << " " << "t_min" << " " << "t_max" << " "
                             << "sig_x" << " " << "sig_y" << " " << "sig_t" << " "
                             << "px_mean" << " " << "px_min" << " " << "px_max" << " "
                             << "py_mean" << " " << "py_min" << " " << "py_max" << " "
                             << "pt_mean" << " " << "pt_min" << " " << "pt_max" << " "
                             << "sig_px" << " " << "sig_py" << " " << "sig_pt" << " "
                             << "emittance_x" << " " << "emittance_y" << " " << "emittance_t" << " "
                             << "alpha_x" << " " << "alpha_y" << " " << "alpha_t" << " "
                             << "beta_x" << " " << "beta_y" << " " << "beta_t" << " "
                             << "dispersion_x" << " " << "dispersion_px" << " "
                             << "dispersion_y" << " " << "dispersion_py" << " "
                             << "charge_C" << " "
                             << "\n";
            }
        }

        if (otype == OutputType::PrintRefParticle) {
            // preparing to access reference particle data: RefPart
            RefPart const ref_part = pc.GetRefParticle();

            amrex::ParticleReal const s = ref_part.s;
            amrex::ParticleReal const beta = ref_part.beta();
            amrex::ParticleReal const gamma = ref_part.gamma();
            amrex::ParticleReal const beta_gamma = ref_part.beta_gamma();
            amrex::ParticleReal const x = ref_part.x;
            amrex::ParticleReal const y = ref_part.y;
            amrex::ParticleReal const z = ref_part.z;
            amrex::ParticleReal const t = ref_part.t;
            amrex::ParticleReal const px = ref_part.px;
            amrex::ParticleReal const py = ref_part.py;
            amrex::ParticleReal const pz = ref_part.pz;
            amrex::ParticleReal const pt = ref_part.pt;

            // write particle data to file
            file_handler
                    << step << " " << s << " "
                    << beta << " " << gamma << " " << beta_gamma << " "
                    << x << " " << y << " " << z << " " << t << " "
                    << px << " " << py << " " << pz << " " << pt << "\n";
        } // if( otype == OutputType::PrintRefParticle)
        else if (otype == OutputType::PrintReducedBeamCharacteristics) {
            std::unordered_map<std::string, amrex::ParticleReal> const rbc =
                diagnostics::reduced_beam_characteristics(pc);

            amrex::ParticleReal const s = pc.GetRefParticle().s;

            file_handler << step << " " << s << " "
                         << rbc.at("x_mean") << " " << rbc.at("x_min") << " " << rbc.at("x_max") << " "
                         << rbc.at("y_mean") << " " << rbc.at("y_min") << " " << rbc.at("y_max") << " "
                         << rbc.at("t_mean") << " " << rbc.at("t_min") << " " << rbc.at("t_max") << " "
                         << rbc.at("sig_x") << " " << rbc.at("sig_y") << " " << rbc.at("sig_t") << " "
                         << rbc.at("px_mean") << " " << rbc.at("px_min") << " " << rbc.at("px_max") << " "
                         << rbc.at("py_mean") << " " << rbc.at("py_min") << " " << rbc.at("py_max") << " "
                         << rbc.at("pt_mean") << " " << rbc.at("pt_min") << " " << rbc.at("pt_max") << " "
                         << rbc.at("sig_px") << " " << rbc.at("sig_py") << " " << rbc.at("sig_pt") << " "
                         << rbc.at("emittance_x") << " " << rbc.at("emittance_y") << " " << rbc.at("emittance_t") << " "
                         << rbc.at("alpha_x") << " " << rbc.at("alpha_y") << " " << rbc.at("alpha_t") << " "
                         << rbc.at("beta_x") << " " << rbc.at("beta_y") << " " << rbc.at("beta_t") << " "
                         << rbc.at("dispersion_x") << " " << rbc.at("dispersion_px") << " "
                         << rbc.at("dispersion_y") << " " << rbc.at("dispersion_py") << " "
                         << rbc.at("charge_C") << "\n";
        } // if( otype == OutputType::PrintReducedBeamCharacteristics)

        // TODO: add as an option to the monitor element
        if (otype == OutputType::PrintNonlinearLensInvariants) {
            // create a host-side particle buffer
            auto tmp = pc.make_alike<amrex::PinnedArenaAllocator>();

            // copy all particles from device to host
            bool const local = true;
            tmp.copyParticles(pc, local);

            // loop over refinement levels
            int const nLevel = tmp.finestLevel();
            for (int lev = 0; lev <= nLevel; ++lev) {
                // loop over all particle boxes
                using MyPinnedParIter = amrex::ParIterSoA<RealSoA::nattribs,IntSoA::nattribs, amrex::PinnedArenaAllocator>;
                for (MyPinnedParIter pti(tmp, lev); pti.isValid(); ++pti) {
                    const int np = pti.numParticles();

                    // preparing access to particle data: SoA of Reals
                    auto const& soa = pti.GetStructOfArrays();
                    auto const& part_x = soa.GetRealData(RealSoA::x);
                    auto const& part_y = soa.GetRealData(RealSoA::y);
                    auto const& part_px = soa.GetRealData(RealSoA::px);
                    auto const& part_py = soa.GetRealData(RealSoA::py);

                    auto const& part_idcpu = soa.GetIdCPUData().dataPtr();

                    // Parse the diagnostic parameters
                    amrex::ParmParse pp_diag("diag");

                    amrex::ParticleReal alpha = 0.0;
                    pp_diag.queryAdd("alpha", alpha);

                    amrex::ParticleReal beta = 1.0;
                    pp_diag.queryAdd("beta", beta);

                    amrex::ParticleReal tn = 0.4;
                    pp_diag.queryAdd("tn", tn);

                    amrex::ParticleReal cn = 0.01;
                    pp_diag.queryAdd("cn", cn);

                    NonlinearLensInvariants const nonlinear_lens_invariants(alpha, beta, tn, cn);

                    // print out particles (this hack works only on CPU and on GPUs with
                    // unified memory access)
                    for (int i = 0; i < np; ++i) {
                        amrex::ParticleReal const x = part_x[i];
                        amrex::ParticleReal const y = part_y[i];
                        uint64_t const global_id = part_idcpu[i];

                        amrex::ParticleReal const px = part_px[i];
                        amrex::ParticleReal const py = part_py[i];

                        // calculate invariants of motion
                        NonlinearLensInvariants::Data const HI_out =
                                nonlinear_lens_invariants(x, y, px, py);

                        // write particle invariant data to file
                        file_handler
                                << global_id << " "
                                << HI_out.H << " " << HI_out.I << "\n";

                    } // i=0...np
                } // end loop over all particle boxes
            } // end mesh-refinement level loop
        }
    }

} // namespace impactx::diagnostics
