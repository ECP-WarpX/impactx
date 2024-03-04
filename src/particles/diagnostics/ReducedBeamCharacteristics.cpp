/* Copyright 2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Marco Garten, Chad Mitchell, Yinjian Zhao, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */

#include "ReducedBeamCharacteristics.H"

#include "particles/ImpactXParticleContainer.H"
#include "particles/ReferenceParticle.H"

#include <AMReX_BLProfiler.H>           // for TinyProfiler
#include <AMReX_GpuQualifiers.H>        // for AMREX_GPU_DEVICE
#include <AMReX_REAL.H>                 // for ParticleReal
#include <AMReX_Reduce.H>               // for ReduceOps
#include <AMReX_ParallelDescriptor.H>   // for ParallelDescriptor
#include <AMReX_ParticleReduce.H>       // for ParticleReduce


namespace impactx::diagnostics
{
    std::unordered_map<std::string, amrex::ParticleReal>
    reduced_beam_characteristics (ImpactXParticleContainer const & pc)
    {
        BL_PROFILE("impactx::diagnostics::reduced_beam_characteristics");

        // preparing to access reference particle data: RefPart
        RefPart const ref_part = pc.GetRefParticle();
        // reference particle charge in C
        amrex::ParticleReal const q_C = ref_part.charge;

        // preparing access to particle data: SoA
        using PType = typename ImpactXParticleContainer::SuperParticleType;

        // prepare reduction operations for calculation of mean and min/max values in 6D phase space
        amrex::ReduceOps<
            // preparing mean values
            amrex::ReduceOpSum,  // w
            amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum,  // x, y, t
            amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum,  // px, py, pt
            // preparing min/max values
            amrex::ReduceOpMin, amrex::ReduceOpMax,  // x_min/max
            amrex::ReduceOpMin, amrex::ReduceOpMax,  // y_min/max
            amrex::ReduceOpMin, amrex::ReduceOpMax,  // t_min/max
            amrex::ReduceOpMin, amrex::ReduceOpMax,  // px_min/max
            amrex::ReduceOpMin, amrex::ReduceOpMax,  // py_min/max
            amrex::ReduceOpMin, amrex::ReduceOpMax   // pt_min/max
        > reduce_ops;

        auto r = amrex::ParticleReduce<
            amrex::ReduceData<
                amrex::ParticleReal,
                amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                amrex::ParticleReal, amrex::ParticleReal,
                amrex::ParticleReal, amrex::ParticleReal,
                amrex::ParticleReal, amrex::ParticleReal,
                amrex::ParticleReal, amrex::ParticleReal,
                amrex::ParticleReal, amrex::ParticleReal,
                amrex::ParticleReal, amrex::ParticleReal
            >
        >(
            pc,
            [=] AMREX_GPU_DEVICE (const PType& p) noexcept
            -> amrex::GpuTuple<
                amrex::ParticleReal,
                amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                amrex::ParticleReal, amrex::ParticleReal,
                amrex::ParticleReal, amrex::ParticleReal,
                amrex::ParticleReal, amrex::ParticleReal,
                amrex::ParticleReal, amrex::ParticleReal,
                amrex::ParticleReal, amrex::ParticleReal,
                amrex::ParticleReal, amrex::ParticleReal
            >
            {
                // access particle position data
                const amrex::ParticleReal p_x = p.rdata(RealSoA::x);
                const amrex::ParticleReal p_y = p.rdata(RealSoA::y);
                const amrex::ParticleReal p_t = p.rdata(RealSoA::t);

                // access SoA particle momentum data and weighting
                const amrex::ParticleReal p_w = p.rdata(RealSoA::w);
                const amrex::ParticleReal p_px = p.rdata(RealSoA::px);
                const amrex::ParticleReal p_py = p.rdata(RealSoA::py);
                const amrex::ParticleReal p_pt = p.rdata(RealSoA::pt);

                // prepare mean position values
                const amrex::ParticleReal p_x_mean = p_x * p_w;
                const amrex::ParticleReal p_y_mean = p_y * p_w;
                const amrex::ParticleReal p_t_mean = p_t * p_w;

                const amrex::ParticleReal p_px_mean = p_px * p_w;
                const amrex::ParticleReal p_py_mean = p_py * p_w;
                const amrex::ParticleReal p_pt_mean = p_pt * p_w;

                return {p_w,
                        p_x_mean, p_y_mean, p_t_mean,
                        p_px_mean, p_py_mean, p_pt_mean,
                        p_x, p_x,
                        p_y, p_y,
                        p_t, p_t,
                        p_px, p_px,
                        p_py, p_py,
                        p_pt, p_pt};
            },
            reduce_ops
        );

        std::vector<amrex::ParticleReal> values_per_rank_1st = {
                amrex::get<0>(r), // w
                amrex::get<1>(r), // x_mean
                amrex::get<2>(r), // y_mean
                amrex::get<3>(r), // t_mean
                amrex::get<4>(r), // px_mean
                amrex::get<5>(r), // py_mean
                amrex::get<6>(r), // pt_mean
        };

        // reduced sum over mpi ranks (allreduce)
        amrex::ParallelAllReduce::Sum(
            values_per_rank_1st.data(),
            values_per_rank_1st.size(),
            amrex::ParallelDescriptor::Communicator()
        );

        amrex::ParticleReal const w_sum   = values_per_rank_1st.at(0);
        amrex::ParticleReal const x_mean  = values_per_rank_1st.at(1) /= w_sum;
        amrex::ParticleReal const y_mean  = values_per_rank_1st.at(2) /= w_sum;
        amrex::ParticleReal const t_mean  = values_per_rank_1st.at(3) /= w_sum;
        amrex::ParticleReal const px_mean = values_per_rank_1st.at(4) /= w_sum;
        amrex::ParticleReal const py_mean = values_per_rank_1st.at(5) /= w_sum;
        amrex::ParticleReal const pt_mean = values_per_rank_1st.at(6) /= w_sum;

        std::vector<amrex::ParticleReal> values_per_rank_min = {
                amrex::get<7>(r), // x_min
                amrex::get<9>(r), // y_min
                amrex::get<11>(r), // t_min
                amrex::get<13>(r), // px_min
                amrex::get<15>(r), // py_min
                amrex::get<17>(r), // pt_min
        };

        std::vector<amrex::ParticleReal> values_per_rank_max = {
                amrex::get<8>(r), // x_max
                amrex::get<10>(r), // y_max
                amrex::get<12>(r), // t_max
                amrex::get<14>(r), // px_max
                amrex::get<16>(r), // py_max
                amrex::get<18>(r), // pt_max
        };

        // reduced sum over mpi ranks (allreduce)
        amrex::ParallelAllReduce::Min(
                values_per_rank_min.data(),
                values_per_rank_min.size(),
                amrex::ParallelDescriptor::Communicator()
        );

        // reduced sum over mpi ranks (allreduce)
        amrex::ParallelAllReduce::Max(
                values_per_rank_max.data(),
                values_per_rank_max.size(),
                amrex::ParallelDescriptor::Communicator()
        );

        // prepare reduction operations for calculation of mean square and correlation values
        amrex::ReduceOps<
            amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum,
            amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum,
            amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum,
            amrex::ReduceOpSum
        > reduce_ops2;


        auto r2 = amrex::ParticleReduce<
            amrex::ReduceData<
                amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                amrex::ParticleReal
            >
        >(
            pc,
            [=] AMREX_GPU_DEVICE(const PType& p) noexcept
            -> amrex::GpuTuple<
                amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                amrex::ParticleReal
            >
            {
                // access SoA particle momentum data and weighting
                const amrex::ParticleReal p_w = p.rdata(RealSoA::w);
                const amrex::ParticleReal p_px = p.rdata(RealSoA::px);
                const amrex::ParticleReal p_py = p.rdata(RealSoA::py);
                const amrex::ParticleReal p_pt = p.rdata(RealSoA::pt);
                // access position data
                const amrex::ParticleReal p_x = p.rdata(RealSoA::x);
                const amrex::ParticleReal p_y = p.rdata(RealSoA::y);
                const amrex::ParticleReal p_t = p.rdata(RealSoA::t);
                // prepare mean square for positions
                const amrex::ParticleReal p_x_ms = (p_x-x_mean)*(p_x-x_mean)*p_w;
                const amrex::ParticleReal p_y_ms = (p_y-y_mean)*(p_y-y_mean)*p_w;
                const amrex::ParticleReal p_t_ms = (p_t-t_mean)*(p_t-t_mean)*p_w;
                // prepare mean square for momenta
                const amrex::ParticleReal p_px_ms = (p_px-px_mean)*(p_px-px_mean)*p_w;
                const amrex::ParticleReal p_py_ms = (p_py-py_mean)*(p_py-py_mean)*p_w;
                const amrex::ParticleReal p_pt_ms = (p_pt-pt_mean)*(p_pt-pt_mean)*p_w;

                const amrex::ParticleReal p_xpx = (p_x-x_mean)*(p_px-px_mean)*p_w;
                const amrex::ParticleReal p_ypy = (p_y-y_mean)*(p_py-py_mean)*p_w;
                const amrex::ParticleReal p_tpt = (p_t-t_mean)*(p_pt-pt_mean)*p_w;

                const amrex::ParticleReal p_charge = q_C*p_w;

                return {p_x_ms, p_y_ms, p_t_ms,
                        p_px_ms, p_py_ms, p_pt_ms,
                        p_xpx, p_ypy, p_tpt,
                        p_charge};
            },
            reduce_ops2
        );

        std::vector<amrex::ParticleReal> values_per_rank_2nd = {
                amrex::get<0>(r2), // x_ms
                amrex::get<1>(r2), // y_ms
                amrex::get<2>(r2), // t_ms
                amrex::get<3>(r2), // px_ms
                amrex::get<4>(r2), // py_ms
                amrex::get<5>(r2), // pt_ms
                amrex::get<6>(r2), // xpx
                amrex::get<7>(r2), // ypy
                amrex::get<8>(r2), // tpt
                amrex::get<9>(r2) // charge
        };

        // reduced sum over mpi ranks (reduce to IO rank)
        amrex::ParallelDescriptor::ReduceRealSum(
            values_per_rank_2nd.data(),
            values_per_rank_2nd.size(),
            amrex::ParallelDescriptor::IOProcessorNumber()
        );

        // minimum values
        amrex::ParticleReal const x_min = values_per_rank_min.at(0);
        amrex::ParticleReal const y_min = values_per_rank_min.at(1);
        amrex::ParticleReal const t_min = values_per_rank_min.at(2);
        amrex::ParticleReal const px_min = values_per_rank_min.at(3);
        amrex::ParticleReal const py_min = values_per_rank_min.at(4);
        amrex::ParticleReal const pt_min = values_per_rank_min.at(5);
        // maximum values
        amrex::ParticleReal const x_max = values_per_rank_max.at(0);
        amrex::ParticleReal const y_max = values_per_rank_max.at(1);
        amrex::ParticleReal const t_max = values_per_rank_max.at(2);
        amrex::ParticleReal const px_max = values_per_rank_max.at(3);
        amrex::ParticleReal const py_max = values_per_rank_max.at(4);
        amrex::ParticleReal const pt_max = values_per_rank_max.at(5);
        // mean square and correlation values
        amrex::ParticleReal const x_ms   = values_per_rank_2nd.at(0) /= w_sum;
        amrex::ParticleReal const y_ms   = values_per_rank_2nd.at(1) /= w_sum;
        amrex::ParticleReal const t_ms   = values_per_rank_2nd.at(2) /= w_sum;
        amrex::ParticleReal const px_ms  = values_per_rank_2nd.at(3) /= w_sum;
        amrex::ParticleReal const py_ms  = values_per_rank_2nd.at(4) /= w_sum;
        amrex::ParticleReal const pt_ms  = values_per_rank_2nd.at(5) /= w_sum;
        amrex::ParticleReal const xpx    = values_per_rank_2nd.at(6) /= w_sum;
        amrex::ParticleReal const ypy    = values_per_rank_2nd.at(7) /= w_sum;
        amrex::ParticleReal const tpt    = values_per_rank_2nd.at(8) /= w_sum;
        amrex::ParticleReal const charge = values_per_rank_2nd.at(9);
        // standard deviations of positions
        amrex::ParticleReal const sig_x = std::sqrt(x_ms);
        amrex::ParticleReal const sig_y = std::sqrt(y_ms);
        amrex::ParticleReal const sig_t = std::sqrt(t_ms);
        // standard deviations of momenta
        amrex::ParticleReal const sig_px = std::sqrt(px_ms);
        amrex::ParticleReal const sig_py = std::sqrt(py_ms);
        amrex::ParticleReal const sig_pt = std::sqrt(pt_ms);
        // RMS emittances
        amrex::ParticleReal const emittance_x = std::sqrt(x_ms*px_ms-xpx*xpx);
        amrex::ParticleReal const emittance_y = std::sqrt(y_ms*py_ms-ypy*ypy);
        amrex::ParticleReal const emittance_t = std::sqrt(t_ms*pt_ms-tpt*tpt);
        // Courant-Snyder (Twiss) beta-function
        amrex::ParticleReal const beta_x = x_ms / emittance_x;
        amrex::ParticleReal const beta_y = y_ms / emittance_y;
        amrex::ParticleReal const beta_t = t_ms / emittance_t;
        // Courant-Snyder (Twiss) alpha
        amrex::ParticleReal const alpha_x = - xpx / emittance_x;
        amrex::ParticleReal const alpha_y = - ypy / emittance_y;
        amrex::ParticleReal const alpha_t = - tpt / emittance_t;

        std::unordered_map<std::string, amrex::ParticleReal> data;
        data["s"] = ref_part.s;  // TODO: remove when the output gets rerouted to openPMD
        data["ref_beta_gamma"] = ref_part.beta_gamma();  // TODO: remove when the output gets rerouted to openPMD
        data["x_mean"] = x_mean;
        data["x_min"] = x_min;
        data["x_max"] = x_max;
        data["y_mean"] = y_mean;
        data["y_min"] = y_min;
        data["y_max"] = y_max;
        data["t_mean"] = t_mean;
        data["t_min"] = t_min;
        data["t_max"] = t_max;
        data["sig_x"] = sig_x;
        data["sig_y"] = sig_y;
        data["sig_t"] = sig_t;
        data["px_mean"] = px_mean;
        data["px_min"] = px_min;
        data["px_max"] = px_max;
        data["py_mean"] = py_mean;
        data["py_min"] = py_min;
        data["py_max"] = py_max;
        data["pt_mean"] = pt_mean;
        data["pt_min"] = pt_min;
        data["pt_max"] = pt_max;
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
        data["charge_C"] = charge;  // TODO: remove when the output gets rerouted to openPMD

        return data;
    }
} // namespace impactx::diagnostics
