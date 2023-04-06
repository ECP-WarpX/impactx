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
#include <Particles/NamedComponentParticleContainer.H>
#include "particles/ImpactXParticleContainer.H"

#include <AMReX_ParallelDescriptor.H>   // for ParallelDescriptor
#include <AMReX_REAL.H>                 // for ParticleReal
#include <AMReX_Reduce.H>               // for ReduceOps
#include <AMReX_ParticleReduce.H>       // for ParticleReduce
#include <AMReX_ParallelDescriptor.H>   // for ParallelDescriptor
#include <AMReX_BLProfiler.H>           // for TinyProfiler

using namespace amrex;

namespace impactx
{

    namespace diagnostics{


        void compute_beam_relevant(ImpactXParticleContainer const & pc, int const step) {

            BL_PROFILE("impactX::diagnostics::compute_beam_relevant");

            // physical constants and reference quantities
            amrex::ParticleReal constexpr c0_SI = 2.99792458e8;  // TODO move out
            // inverse of speed of light squared
            amrex::Real constexpr inv_c2 = 1.0_rt / (c0_SI * c0_SI);

            // preparing to access reference particle data: RefPart
            RefPart const ref_part = pc.GetRefParticle();

            ParticleReal const m_MeV = ref_part.mass_MeV();
            ParticleReal const q_qe = ref_part.charge_qe();

            // preparing access to particle data: AoS and SoA
            using PTDType = typename ImpactXParticleContainer::ParticleTileType::ConstParticleTileDataType;

            amrex::ReduceOps<ReduceOpSum,ReduceOpSum,ReduceOpSum,ReduceOpSum,ReduceOpSum,
            ReduceOpSum,ReduceOpSum,ReduceOpSum> reduce_ops;

            auto r = amrex::ParticleReduce<amrex::ReduceData<ParticleReal,ParticleReal,
            ParticleReal,ParticleReal,ParticleReal,ParticleReal,ParticleReal,ParticleReal>>(
                    pc,
                    [=] AMREX_GPU_DEVICE(const PTDType& ptd, const int i) noexcept -> amrex::GpuTuple
                            <ParticleReal,ParticleReal,ParticleReal,ParticleReal,ParticleReal,
                            ParticleReal,ParticleReal,ParticleReal>
                    {
                        // access AoS particle position data
                        const ParticleReal p_pos0 = ptd.m_aos[i].pos(0);
                        const ParticleReal p_pos1 = ptd.m_aos[i].pos(1);
                        const ParticleReal p_pos2 = ptd.m_aos[i].pos(2);

                        // access SoA particle momentum data and weighting
                        const ParticleReal p_w = ptd.m_rdata[RealSoA::w][i];
                        const ParticleReal p_ux = ptd.m_rdata[RealSoA::ux][i];
                        const ParticleReal p_uy = ptd.m_rdata[RealSoA::uy][i];
                        const ParticleReal p_ut = ptd.m_rdata[RealSoA::pt][i];

                        const ParticleReal p_us = p_ux*p_ux + p_uy*p_uy + p_ut*p_ut;
                        // prepare mean position values
                        const ParticleReal p_x_mean = p_pos0*p_w;
                        const ParticleReal p_y_mean = p_pos1*p_w;
                        const ParticleReal p_t_mean = p_pos2*p_w;

                        const ParticleReal p_ux_mean = p_ux*p_w;
                        const ParticleReal p_uy_mean = p_uy*p_w;
                        const ParticleReal p_ut_mean = p_ut*p_w;
                        const ParticleReal p_gm_mean = std::sqrt(1.0_rt+p_us*inv_c2)*p_w;

                        return {p_w,
                                p_x_mean, p_y_mean, p_t_mean,
                                p_ux_mean, p_uy_mean, p_ut_mean,
                                p_gm_mean};
                    },
                    reduce_ops);

            std::vector<ParticleReal> values_per_rank_1st = {
                    amrex::get<0>(r), // w
                    amrex::get<1>(r), // x_mean
                    amrex::get<2>(r), // y_mean
                    amrex::get<3>(r), // z_mean
                    amrex::get<4>(r), // ux_mean
                    amrex::get<5>(r), // uy_mean
                    amrex::get<6>(r), // ut_mean
                    amrex::get<7>(r), // gm_mean
            };

            // reduced sum over mpi ranks (allreduce)
            amrex::ParallelAllReduce::Sum
                    ( values_per_rank_1st.data(), values_per_rank_1st.size(), ParallelDescriptor::Communicator());

            ParticleReal w_sum   = values_per_rank_1st.at(0);
            ParticleReal x_mean  = values_per_rank_1st.at(1) /= w_sum;
            ParticleReal y_mean  = values_per_rank_1st.at(2) /= w_sum;
            ParticleReal t_mean  = values_per_rank_1st.at(3) /= w_sum;
            ParticleReal ux_mean = values_per_rank_1st.at(4) /= w_sum;
            ParticleReal uy_mean = values_per_rank_1st.at(5) /= w_sum;
            ParticleReal ut_mean = values_per_rank_1st.at(6) /= w_sum;
            ParticleReal gm_mean = values_per_rank_1st.at(7) /= w_sum;


            amrex::ReduceOps<ReduceOpSum,ReduceOpSum,ReduceOpSum,ReduceOpSum,ReduceOpSum,
                    ReduceOpSum,ReduceOpSum,ReduceOpSum,ReduceOpSum,ReduceOpSum,ReduceOpSum> reduce_ops2;


            auto r2 = amrex::ParticleReduce<amrex::ReduceData<ParticleReal,ParticleReal,
                    ParticleReal,ParticleReal,ParticleReal,ParticleReal,ParticleReal,ParticleReal,
                    ParticleReal,ParticleReal,ParticleReal>>(
                    pc,
                    [=] AMREX_GPU_DEVICE(const PTDType& ptd, const int i) noexcept -> amrex::GpuTuple
                            <ParticleReal,ParticleReal,ParticleReal,ParticleReal,ParticleReal,ParticleReal,
                                    ParticleReal,ParticleReal,ParticleReal,ParticleReal,ParticleReal>
                    {
                        const ParticleReal p_ux = ptd.m_rdata[RealSoA::ux][i];
                        const ParticleReal p_uy = ptd.m_rdata[RealSoA::uy][i];
                        const ParticleReal p_ut = ptd.m_rdata[RealSoA::pt][i];
                        const ParticleReal p_us = p_ux*p_ux + p_uy*p_uy + p_ut*p_ut;
                        const ParticleReal p_gm = std::sqrt(1.0_rt+p_us*inv_c2);
                        const ParticleReal p_w = ptd.m_rdata[RealSoA::w][i];


                        const ParticleReal p_pos0 = ptd.m_aos[i].pos(0);
                        const ParticleReal p_pos1 = ptd.m_aos[i].pos(1);
                        const ParticleReal p_pos2 = ptd.m_aos[i].pos(2);
                        const ParticleReal p_x = p_pos0;
                        const ParticleReal p_y = p_pos1;
                        const ParticleReal p_t = p_pos2;
                        // prepare mean square for positions
                        const ParticleReal p_x_ms = (p_x-x_mean)*(p_x-x_mean)*p_w;
                        const ParticleReal p_y_ms = (p_y-y_mean)*(p_y-y_mean)*p_w;
                        const ParticleReal p_z_ms = (p_t-t_mean)*(p_t-t_mean)*p_w;
                        // prepare mean square for momenta
                        const ParticleReal p_ux_ms = (p_ux-ux_mean)*(p_ux-ux_mean)*p_w;
                        const ParticleReal p_uy_ms = (p_uy-uy_mean)*(p_uy-uy_mean)*p_w;
                        const ParticleReal p_uz_ms = (p_ut-ut_mean)*(p_ut-ut_mean)*p_w;
                        const ParticleReal p_gm_ms = (p_gm-gm_mean)*(p_gm-gm_mean)*p_w;

                        const ParticleReal p_xux = (p_x-x_mean)*(p_ux-ux_mean)*p_w;
                        const ParticleReal p_yuy = (p_y-y_mean)*(p_uy-uy_mean)*p_w;
                        const ParticleReal p_zuz = (p_t-t_mean)*(p_ut-ut_mean)*p_w;

                        const ParticleReal p_charge = q_qe*p_w;

                        return {p_x_ms, p_y_ms, p_z_ms,
                                p_ux_ms, p_uy_ms, p_uz_ms,
                                p_gm_ms,
                                p_xux, p_yuy, p_zuz,
                                p_charge};
                    },
                    reduce_ops2);

            std::vector<ParticleReal> values_per_rank_2nd = {
                    amrex::get<0>(r2), // x_ms
                    amrex::get<1>(r2), // y_ms
                    amrex::get<2>(r2), // z_ms
                    amrex::get<3>(r2), // ux_ms
                    amrex::get<4>(r2), // uy_ms
                    amrex::get<5>(r2), // uz_ms
                    amrex::get<6>(r2), // gm_ms
                    amrex::get<7>(r2), // xux
                    amrex::get<8>(r2), // yuy
                    amrex::get<9>(r2), // zuz
                    amrex::get<10>(r2) // charge
            };

            // reduced sum over mpi ranks (reduce to IO rank)
            ParallelDescriptor::ReduceRealSum
                    ( values_per_rank_2nd.data(), values_per_rank_2nd.size(), ParallelDescriptor::IOProcessorNumber());

            ParticleReal x_ms   = values_per_rank_2nd.at(0) /= w_sum;
            ParticleReal y_ms   = values_per_rank_2nd.at(1) /= w_sum;
            ParticleReal t_ms   = values_per_rank_2nd.at(2) /= w_sum;
            ParticleReal ux_ms  = values_per_rank_2nd.at(3) /= w_sum;
            ParticleReal uy_ms  = values_per_rank_2nd.at(4) /= w_sum;
            ParticleReal ut_ms  = values_per_rank_2nd.at(5) /= w_sum;
            ParticleReal gm_ms  = values_per_rank_2nd.at(6) /= w_sum;
            ParticleReal xux    = values_per_rank_2nd.at(7) /= w_sum;
            ParticleReal yuy    = values_per_rank_2nd.at(8) /= w_sum;
            ParticleReal tut    = values_per_rank_2nd.at(9) /= w_sum;
            ParticleReal charge = values_per_rank_2nd.at(10);
            // standard deviations of positions
            ParticleReal sig_x = std::sqrt(x_ms);
            ParticleReal sig_y = std::sqrt(y_ms);
            ParticleReal sig_t = std::sqrt(t_ms);
            // standard deviations of momenta
            ParticleReal sig_px = std::sqrt(ux_ms);
            ParticleReal sig_py = std::sqrt(uy_ms);
            ParticleReal sig_pt = std::sqrt(ut_ms);
            // RMS emittances
            ParticleReal emittance_x = std::sqrt(x_ms*ux_ms-xux*xux);
            ParticleReal emittance_y = std::sqrt(y_ms*uy_ms-yuy*yuy);
            ParticleReal emittance_z = std::sqrt(t_ms*ut_ms-tut*tut);

            amrex::AllPrint() << "step" << " "
                              << "x_mean" << " " << "y_mean" << " " << "t_mean" << " "
                              << "sig_x" << " " << "sig_y" << " " << "sig_t" << " "
                              << "ux_mean" << " " << "uy_mean" << " " << "ut_mean" << " "
                              << "sig_px" << " " << "sig_py" << " " << "sig_pt" << " "
                              << "emittance_x" << " " << "emittance_y" << " " << "emittance_z" << " "
                              << "charge" << " "
                              << "\n";

            amrex::AllPrint() << step << " "
                              << x_mean << " " << y_mean << " " << t_mean << " "
                              << sig_x << " " << sig_y << " " << sig_t << " "
                              << ux_mean << " " << uy_mean << " " << ut_mean << " "
                              << sig_px << " " << sig_py << " " << sig_pt << " "
                              << emittance_x << " " << emittance_y << " " << emittance_z << " "
                               << charge << " "
                               << "\n";

        }

    }
}
