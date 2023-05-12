/* Copyright 2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Marco Garten, Yinjian Zhao
 * License: BSD-3-Clause-LBNL
 */

#include "ImpactX.H"
#include "particles/NamedComponentParticleContainer.H"
#include "particles/ImpactXParticleContainer.H"
#include "particles/diagnostics/ReducedBeamCharacteristics.H"

#include <AMReX_ParallelDescriptor.H>   // for ParallelDescriptor
#include <AMReX_REAL.H>                 // for ParticleReal
#include <AMReX_Reduce.H>               // for ReduceOps
#include <AMReX_ParticleReduce.H>       // for ParticleReduce
#include <AMReX_ParallelDescriptor.H>   // for ParallelDescriptor
#include <AMReX_BLProfiler.H>           // for TinyProfiler


namespace impactx
{

    namespace diagnostics{


        //ReducedBeamCharacteristics reduced_beam_characteristics(ImpactXParticleContainer const & pc) {
        std::unordered_map<std::string, amrex::ParticleReal> reduced_beam_characteristics(ImpactXParticleContainer const & pc) {

            BL_PROFILE("impactX::diagnostics::reduced_beam_characteristics");

            // preparing to access reference particle data: RefPart
            RefPart const ref_part = pc.GetRefParticle();
            // reference particle charge in units of positive elementary charge
            amrex::ParticleReal const q_qe = ref_part.charge_qe();

            // preparing access to particle data: AoS and SoA
            using PTDType = typename ImpactXParticleContainer::ParticleTileType::ConstParticleTileDataType;

            amrex::ReduceOps<amrex::ReduceOpSum,amrex::ReduceOpSum,amrex::ReduceOpSum,amrex::ReduceOpSum,amrex::ReduceOpSum,
            amrex::ReduceOpSum,amrex::ReduceOpSum> reduce_ops;

            auto r = amrex::ParticleReduce<amrex::ReduceData<amrex::ParticleReal,amrex::ParticleReal,
            amrex::ParticleReal,amrex::ParticleReal,amrex::ParticleReal,amrex::ParticleReal,amrex::ParticleReal>>(
                    pc,
                    [=] AMREX_GPU_DEVICE(const PTDType& ptd, const int i) noexcept -> amrex::GpuTuple
                            <amrex::ParticleReal,amrex::ParticleReal,amrex::ParticleReal,amrex::ParticleReal,amrex::ParticleReal,
                            amrex::ParticleReal,amrex::ParticleReal>
                    {
                        // access AoS particle position data
                        const amrex::ParticleReal p_pos0 = ptd.m_aos[i].pos(0);
                        const amrex::ParticleReal p_pos1 = ptd.m_aos[i].pos(1);
                        const amrex::ParticleReal p_pos2 = ptd.m_aos[i].pos(2);

                        // access SoA particle momentum data and weighting
                        const amrex::ParticleReal p_w = ptd.m_rdata[RealSoA::w][i];
                        const amrex::ParticleReal p_px = ptd.m_rdata[RealSoA::px][i];
                        const amrex::ParticleReal p_py = ptd.m_rdata[RealSoA::py][i];
                        const amrex::ParticleReal p_ut = ptd.m_rdata[RealSoA::pt][i];

                        // prepare mean position values
                        const amrex::ParticleReal p_x_mean = p_pos0*p_w;
                        const amrex::ParticleReal p_y_mean = p_pos1*p_w;
                        const amrex::ParticleReal p_t_mean = p_pos2*p_w;

                        const amrex::ParticleReal p_px_mean = p_px*p_w;
                        const amrex::ParticleReal p_py_mean = p_py*p_w;
                        const amrex::ParticleReal p_ut_mean = p_ut*p_w;;

                        return {p_w,
                                p_x_mean, p_y_mean, p_t_mean,
                                p_px_mean, p_py_mean, p_ut_mean};
                    },
                    reduce_ops);

            std::vector<amrex::ParticleReal> values_per_rank_1st = {
                    amrex::get<0>(r), // w
                    amrex::get<1>(r), // x_mean
                    amrex::get<2>(r), // y_mean
                    amrex::get<3>(r), // z_mean
                    amrex::get<4>(r), // px_mean
                    amrex::get<5>(r), // py_mean
                    amrex::get<6>(r), // ut_mean
            };

            // reduced sum over mpi ranks (allreduce)
            amrex::ParallelAllReduce::Sum
                    ( values_per_rank_1st.data(), values_per_rank_1st.size(), amrex::ParallelDescriptor::Communicator());

            amrex::ParticleReal w_sum   = values_per_rank_1st.at(0);
            amrex::ParticleReal x_mean  = values_per_rank_1st.at(1) /= w_sum;
            amrex::ParticleReal y_mean  = values_per_rank_1st.at(2) /= w_sum;
            amrex::ParticleReal t_mean  = values_per_rank_1st.at(3) /= w_sum;
            amrex::ParticleReal px_mean = values_per_rank_1st.at(4) /= w_sum;
            amrex::ParticleReal py_mean = values_per_rank_1st.at(5) /= w_sum;
            amrex::ParticleReal ut_mean = values_per_rank_1st.at(6) /= w_sum;


            amrex::ReduceOps<amrex::ReduceOpSum,amrex::ReduceOpSum,amrex::ReduceOpSum,amrex::ReduceOpSum,amrex::ReduceOpSum,
                    amrex::ReduceOpSum,amrex::ReduceOpSum,amrex::ReduceOpSum,amrex::ReduceOpSum,amrex::ReduceOpSum> reduce_ops2;


            auto r2 = amrex::ParticleReduce<amrex::ReduceData<amrex::ParticleReal,amrex::ParticleReal,
                    amrex::ParticleReal,amrex::ParticleReal,amrex::ParticleReal,amrex::ParticleReal,amrex::ParticleReal,amrex::ParticleReal,
                    amrex::ParticleReal,amrex::ParticleReal>>(
                    pc,
                    [=] AMREX_GPU_DEVICE(const PTDType& ptd, const int i) noexcept -> amrex::GpuTuple
                            <amrex::ParticleReal,amrex::ParticleReal,amrex::ParticleReal,amrex::ParticleReal,amrex::ParticleReal,amrex::ParticleReal,
                                    amrex::ParticleReal,amrex::ParticleReal,amrex::ParticleReal,amrex::ParticleReal>
                    {
                        const amrex::ParticleReal p_px = ptd.m_rdata[RealSoA::px][i];
                        const amrex::ParticleReal p_py = ptd.m_rdata[RealSoA::py][i];
                        const amrex::ParticleReal p_ut = ptd.m_rdata[RealSoA::pt][i];
                        const amrex::ParticleReal p_w = ptd.m_rdata[RealSoA::w][i];

                        const amrex::ParticleReal p_pos0 = ptd.m_aos[i].pos(0);
                        const amrex::ParticleReal p_pos1 = ptd.m_aos[i].pos(1);
                        const amrex::ParticleReal p_pos2 = ptd.m_aos[i].pos(2);
                        const amrex::ParticleReal p_x = p_pos0;
                        const amrex::ParticleReal p_y = p_pos1;
                        const amrex::ParticleReal p_t = p_pos2;
                        // prepare mean square for positions
                        const amrex::ParticleReal p_x_ms = (p_x-x_mean)*(p_x-x_mean)*p_w;
                        const amrex::ParticleReal p_y_ms = (p_y-y_mean)*(p_y-y_mean)*p_w;
                        const amrex::ParticleReal p_t_ms = (p_t-t_mean)*(p_t-t_mean)*p_w;
                        // prepare mean square for momenta
                        const amrex::ParticleReal p_px_ms = (p_px-px_mean)*(p_px-px_mean)*p_w;
                        const amrex::ParticleReal p_py_ms = (p_py-py_mean)*(p_py-py_mean)*p_w;
                        const amrex::ParticleReal p_ut_ms = (p_ut-ut_mean)*(p_ut-ut_mean)*p_w;

                        const amrex::ParticleReal p_xpx = (p_x-x_mean)*(p_px-px_mean)*p_w;
                        const amrex::ParticleReal p_ypy = (p_y-y_mean)*(p_py-py_mean)*p_w;
                        const amrex::ParticleReal p_zuz = (p_t-t_mean)*(p_ut-ut_mean)*p_w;

                        const amrex::ParticleReal p_charge = q_qe*p_w;

                        return {p_x_ms, p_y_ms, p_t_ms,
                                p_px_ms, p_py_ms, p_ut_ms,
                                p_xpx, p_ypy, p_zuz,
                                p_charge};
                    },
                    reduce_ops2);

            std::vector<amrex::ParticleReal> values_per_rank_2nd = {
                    amrex::get<0>(r2), // x_ms
                    amrex::get<1>(r2), // y_ms
                    amrex::get<2>(r2), // t_ms
                    amrex::get<3>(r2), // px_ms
                    amrex::get<4>(r2), // py_ms
                    amrex::get<5>(r2), // ut_ms
                    amrex::get<6>(r2), // xpx
                    amrex::get<7>(r2), // ypy
                    amrex::get<8>(r2), // tut
                    amrex::get<9>(r2) // charge
            };

            // reduced sum over mpi ranks (reduce to IO rank)
            amrex::ParallelDescriptor::ReduceRealSum
                    ( values_per_rank_2nd.data(), values_per_rank_2nd.size(), amrex::ParallelDescriptor::IOProcessorNumber());

            amrex::ParticleReal x_ms   = values_per_rank_2nd.at(0) /= w_sum;
            amrex::ParticleReal y_ms   = values_per_rank_2nd.at(1) /= w_sum;
            amrex::ParticleReal t_ms   = values_per_rank_2nd.at(2) /= w_sum;
            amrex::ParticleReal px_ms  = values_per_rank_2nd.at(3) /= w_sum;
            amrex::ParticleReal py_ms  = values_per_rank_2nd.at(4) /= w_sum;
            amrex::ParticleReal ut_ms  = values_per_rank_2nd.at(5) /= w_sum;
            amrex::ParticleReal xpx    = values_per_rank_2nd.at(6) /= w_sum;
            amrex::ParticleReal ypy    = values_per_rank_2nd.at(7) /= w_sum;
            amrex::ParticleReal tut    = values_per_rank_2nd.at(8) /= w_sum;
            amrex::ParticleReal charge = values_per_rank_2nd.at(9);
            // standard deviations of positions
            amrex::ParticleReal sig_x = std::sqrt(x_ms);
            amrex::ParticleReal sig_y = std::sqrt(y_ms);
            amrex::ParticleReal sig_t = std::sqrt(t_ms);
            // standard deviations of momenta
            amrex::ParticleReal sig_px = std::sqrt(px_ms);
            amrex::ParticleReal sig_py = std::sqrt(py_ms);
            amrex::ParticleReal sig_pt = std::sqrt(ut_ms);
            // RMS emittances
            amrex::ParticleReal emittance_x = std::sqrt(x_ms*px_ms-xpx*xpx);
            amrex::ParticleReal emittance_y = std::sqrt(y_ms*py_ms-ypy*ypy);
            amrex::ParticleReal emittance_t = std::sqrt(t_ms*ut_ms-tut*tut);
            // Courant-Snyder (Twiss) beta-function
            amrex::ParticleReal beta_x = x_ms / emittance_x;
            amrex::ParticleReal beta_y = y_ms / emittance_y;
            amrex::ParticleReal beta_t = t_ms / emittance_t;
            // Courant-Snyder (Twiss) alpha
            amrex::ParticleReal alpha_x = - xpx / emittance_x;
            amrex::ParticleReal alpha_y = - ypy / emittance_y;
            amrex::ParticleReal alpha_t = - tut / emittance_t;

            std::unordered_map<std::string, amrex::ParticleReal> data;
            data["s"] = ref_part.s;  // TODO: remove when the output gets rerouted to openPMD
            data["beta_gamma"] = ref_part.beta_gamma();  // TODO: remove when the output gets rerouted to openPMD
            data["x_mean"] = x_mean;
            data["y_mean"] = y_mean;
            data["t_mean"] = t_mean;
            data["sig_x"] = sig_x;
            data["sig_y"] = sig_y;
            data["sig_t"] = sig_t;
            data["px_mean"] = px_mean;
            data["py_mean"] = py_mean;
            data["ut_mean"] = ut_mean;
            data["sig_px"] = sig_px;
            data["sig_py"] = sig_py;
            data["sig_pt"] = sig_pt;
            data["emittance_x"] = emittance_x;
            data["emittance_y"] = emittance_y;
            data["emittance_t"] = emittance_t;
            data["alpha_x"] = alpha_x;
            data["alpha_y"] = alpha_y;
            data["alpha_t"] = alpha_t;
            data["beta_x"] = beta_x;
            data["beta_y"] = beta_y;
            data["beta_t"] = beta_t;
            data["charge"] = charge;  // TODO: remove when the output gets rerouted to openPMD

            return data;

        }

    }
}
