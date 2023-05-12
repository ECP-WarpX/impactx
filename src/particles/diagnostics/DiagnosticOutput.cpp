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

#include <ablastr/particles/IndexHandling.H>

#include <AMReX_BLProfiler.H> // for BL_PROFILE
#include <AMReX_Extension.H>  // for AMREX_RESTRICT
#include <AMReX_ParmParse.H>  // for ParmParse
#include <AMReX_REAL.H>       // for ParticleReal
#include <AMReX_Print.H>      // for PrintToFile


namespace impactx::diagnostics
{
    void DiagnosticOutput (ImpactXParticleContainer const & pc,
                           OutputType const otype,
                           std::string file_name,
                           int const step,
                           bool const append)
    {
        BL_PROFILE("impactx::diagnostics::DiagnosticOutput");

        using namespace amrex::literals; // for _rt and _prt

        // keep file open as we add more and more lines
        amrex::AllPrintToFile file_handler(file_name);

        // write file header per MPI RANK
        if (!append) {
            if (otype == OutputType::PrintParticles) {
                file_handler << "id x y t px py pt\n";
            } else if (otype == OutputType::PrintNonlinearLensInvariants) {
                file_handler << "id H I\n";
            } else if (otype == OutputType::PrintRefParticle) {
                file_handler << "step s x y z t px py pz pt\n";
            } else if (otype == OutputType::PrintReducedBeamCharacteristics) {
                file_handler << "step" << " " << "s_ref_part" << " " << "beta_gamma_ref_part" << " "
                             << "x_mean" << " " << "y_mean" << " " << "t_mean" << " "
                             << "sig_x" << " " << "sig_y" << " " << "sig_t" << " "
                             << "px_mean" << " " << "py_mean" << " " << "pt_mean" << " "
                             << "sig_px" << " " << "sig_py" << " " << "sig_pt" << " "
                             << "emittance_x" << " " << "emittance_y" << " " << "emittance_t" << " "
                             << "alpha_x" << " " << "alpha_y" << " " << "alpha_t" << " "
                             << "beta_x" << " " << "beta_y" << " " << "beta_t" << " "
                             << "charge" << " "
                             << "\n";
            }
        }

        if (otype == OutputType::PrintReducedBeamCharacteristics) {
            std::unordered_map<std::string, amrex::ParticleReal> const rbc = diagnostics::reduced_beam_characteristics(
                    pc);

            file_handler << step << " " << rbc.at("s") << " " << rbc.at("beta_gamma") << " "
                         << rbc.at("x_mean") << " " << rbc.at("y_mean") << " " << rbc.at("t_mean") << " "
                         << rbc.at("sig_x") << " " << rbc.at("sig_y") << " " << rbc.at("sig_t") << " "
                         << rbc.at("px_mean") << " " << rbc.at("py_mean") << " " << rbc.at("pt_mean") << " "
                         << rbc.at("sig_px") << " " << rbc.at("sig_py") << " " << rbc.at("sig_pt") << " "
                         << rbc.at("emittance_x") << " " << rbc.at("emittance_y") << " " << rbc.at("emittance_t") << " "
                         << rbc.at("alpha_x") << " " << rbc.at("alpha_y") << " " << rbc.at("alpha_t") << " "
                         << rbc.at("beta_x") << " " << rbc.at("beta_y") << " " << rbc.at("beta_t") << " "
                         << rbc.at("charge") << "\n";
        } // if( otype == OutputType::PrintReducedBeamCharacteristics)

        // create a host-side particle buffer
        auto tmp = pc.make_alike<amrex::PinnedArenaAllocator>();

        // copy device-to-host
        bool const local = true;
        tmp.copyParticles(pc, local);

        // loop over refinement levels
        int const nLevel = tmp.finestLevel();
        for (int lev = 0; lev <= nLevel; ++lev) {
            // loop over all particle boxes
            using ParIt = typename decltype(tmp)::ParConstIterType;
            for (ParIt pti(tmp, lev); pti.isValid(); ++pti) {
                const int np = pti.numParticles();

                // preparing access to particle data: AoS
                using PType = ImpactXParticleContainer::ParticleType;
                auto const &aos = pti.GetArrayOfStructs();
                PType const *const AMREX_RESTRICT aos_ptr = aos().dataPtr();

                // preparing access to particle data: SoA of Reals
                auto &soa_real = pti.GetStructOfArrays().GetRealData();
                amrex::ParticleReal const *const AMREX_RESTRICT part_px = soa_real[RealSoA::px].dataPtr();
                amrex::ParticleReal const *const AMREX_RESTRICT part_py = soa_real[RealSoA::py].dataPtr();
                amrex::ParticleReal const *const AMREX_RESTRICT part_pt = soa_real[RealSoA::pt].dataPtr();

                if (otype == OutputType::PrintParticles) {
                    // print out particles (this hack works only on CPU and on GPUs with
                    // unified memory access)
                    for (int i = 0; i < np; ++i) {

                        // access AoS data such as positions and cpu/id
                        PType const &p = aos_ptr[i];
                        amrex::ParticleReal const x = p.pos(RealAoS::x);
                        amrex::ParticleReal const y = p.pos(RealAoS::y);
                        amrex::ParticleReal const t = p.pos(RealAoS::t);
                        uint64_t const global_id = ablastr::particles::localIDtoGlobal(p.id(), p.cpu());

                        // access SoA Real data
                        amrex::ParticleReal const px = part_px[i];
                        amrex::ParticleReal const py = part_py[i];
                        amrex::ParticleReal const pt = part_pt[i];

                        // write particle data to file
                        file_handler
                                << global_id << " "
                                << x << " " << y << " " << t << " "
                                << px << " " << py << " " << pt << "\n";
                    } // i=0...np
                } // if( otype == OutputType::PrintParticles)
                else if (otype == OutputType::PrintNonlinearLensInvariants) {

                    using namespace amrex::literals;

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

                        // access AoS data such as positions and cpu/id
                        PType const &p = aos_ptr[i];
                        amrex::ParticleReal const x = p.pos(RealAoS::x);
                        amrex::ParticleReal const y = p.pos(RealAoS::y);
                        uint64_t const global_id = ablastr::particles::localIDtoGlobal(p.id(), p.cpu());

                        // access SoA Real data
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
                } // if( otype == OutputType::PrintInvariants)
                if (otype == OutputType::PrintRefParticle) {
                    // print reference particle to file

                    // preparing to access reference particle data: RefPart
                    RefPart const ref_part = pc.GetRefParticle();

                    amrex::ParticleReal const s = ref_part.s;
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
                            << x << " " << y << " " << z << " " << t << " "
                            << px << " " << py << " " << pz << " " << pt << "\n";
                } // if( otype == OutputType::PrintRefParticle)
            } // end loop over all particle boxes
        } // env mesh-refinement level loop
    }

} // namespace impactx::diagnostics
