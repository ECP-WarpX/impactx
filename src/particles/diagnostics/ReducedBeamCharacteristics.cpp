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
#include <AMReX_TypeList.H>             // for TypeMultiplier


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

        /* The variables below need to be static to work around an MSVC bug
         * https://stackoverflow.com/questions/55136414/constexpr-variable-captured-inside-lambda-loses-its-constexpr-ness
         */
        // numbers of same type reduction operations in first concurrent batch
        static constexpr std::size_t num_red_ops_1_sum = 7;  // summation
        static constexpr std::size_t num_red_ops_1_min = 6;  // minimum
        static constexpr std::size_t num_red_ops_1_max = 6;  // maximum

        // prepare reduction operations for calculation of mean and min/max values in 6D phase space
        amrex::TypeMultiplier<amrex::ReduceOps,
            amrex::ReduceOpSum[num_red_ops_1_sum],  // preparing mean values for w, x, y, t, px, py, pt
            amrex::ReduceOpMin[num_red_ops_1_min],  // preparing min values for x, y, t, px, py, pt
            amrex::ReduceOpMax[num_red_ops_1_max]   // preparing max values for x, y, t, px, py, pt
        > reduce_ops_1;
        using ReducedDataT1 = amrex::TypeMultiplier<amrex::ReduceData, amrex::ParticleReal[num_red_ops_1_sum + num_red_ops_1_min + num_red_ops_1_max]>;

        auto r1 = amrex::ParticleReduce<ReducedDataT1>(
            pc,
            [=] AMREX_GPU_DEVICE(const PType& p) noexcept -> ReducedDataT1::Type
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
                        p_x, p_y, p_t, p_px, p_py, p_pt,
                        p_x, p_y, p_t, p_px, p_py, p_pt};
            },
            reduce_ops_1
        );

        std::vector<amrex::ParticleReal> values_per_rank_1st(num_red_ops_1_sum);

        /* contains in this order:
         * w, x_mean, y_mean, t_mean
         * px_mean, py_mean, pt_mean
         */
        amrex::constexpr_for<0, num_red_ops_1_sum> ([&](auto i) {
            values_per_rank_1st[i] = amrex::get<i>(r1);
        });

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

        std::vector<amrex::ParticleReal> values_per_rank_min(num_red_ops_1_min);

        /* contains in this order:
         * x_min, y_min, t_min
         * px_min, py_min, pt_min
         */
        amrex::constexpr_for<0, num_red_ops_1_min> ([&](auto i) {
            constexpr std::size_t idx = i + num_red_ops_1_sum;
            values_per_rank_min[i] = amrex::get<idx>(r1);
        });

        std::vector<amrex::ParticleReal> values_per_rank_max(num_red_ops_1_max);

        /* contains in this order:
         * x_max, y_max, t_max
         * px_max, py_max, pt_max
         */
        amrex::constexpr_for<0, num_red_ops_1_max> ([&](auto i) {
            constexpr std::size_t idx = i + num_red_ops_1_sum + num_red_ops_1_min;
            values_per_rank_max[i] = amrex::get<idx>(r1);
        });

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

        /* The variable below needs to be static to work around an MSVC bug
         * https://stackoverflow.com/questions/55136414/constexpr-variable-captured-inside-lambda-loses-its-constexpr-ness
         */
        // number of reduction operations in second concurrent batch
        static constexpr std::size_t num_red_ops_2 = 14;
        // prepare reduction operations for calculation of mean square and correlation values
        amrex::TypeMultiplier<amrex::ReduceOps, amrex::ReduceOpSum[num_red_ops_2]> reduce_ops_2;
        using ReducedDataT2 = amrex::TypeMultiplier<amrex::ReduceData, amrex::ParticleReal[num_red_ops_2]>;

        auto r2 = amrex::ParticleReduce<ReducedDataT2>(
                pc,
                [=] AMREX_GPU_DEVICE(const PType& p) noexcept
            -> ReducedDataT2::Type
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
                // prepare position-momentum correlations
                const amrex::ParticleReal p_xpx = (p_x-x_mean)*(p_px-px_mean)*p_w;
                const amrex::ParticleReal p_ypy = (p_y-y_mean)*(p_py-py_mean)*p_w;
                const amrex::ParticleReal p_tpt = (p_t-t_mean)*(p_pt-pt_mean)*p_w;
                // prepare correlations for dispersion
                const amrex::ParticleReal p_xpt = (p_x-x_mean)*(p_pt-pt_mean)*p_w;
                const amrex::ParticleReal p_pxpt = (p_px-px_mean)*(p_pt-pt_mean)*p_w;
                const amrex::ParticleReal p_ypt = (p_y-y_mean)*(p_pt-pt_mean)*p_w;
                const amrex::ParticleReal p_pypt = (p_py-py_mean)*(p_pt-pt_mean)*p_w;


                const amrex::ParticleReal p_charge = q_C*p_w;

                return {p_x_ms, p_y_ms, p_t_ms,
                        p_px_ms, p_py_ms, p_pt_ms,
                        p_xpx, p_ypy, p_tpt,
                        p_xpt, p_pxpt, p_ypt, p_pypt,
                        p_charge};
            },
                reduce_ops_2
        );

        std::vector<amrex::ParticleReal> values_per_rank_2nd(num_red_ops_2);

        /* contains in this order:
         * x_ms, y_ms, t_ms
         * px_ms, py_ms, pt_ms,
         * xpx, ypy, tpt,
         * p_xpt, p_pxpt, p_ypt, p_pypt,
         * charge
         */
        amrex::constexpr_for<0, num_red_ops_2> ([&](auto i) {
            values_per_rank_2nd[i] = amrex::get<i>(r2);
        });

        // reduced sum over mpi ranks (allreduce)
        amrex::ParallelAllReduce::Sum(
            values_per_rank_2nd.data(),
            values_per_rank_2nd.size(),
            amrex::ParallelDescriptor::Communicator()
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
        amrex::ParticleReal const xpt    = values_per_rank_2nd.at(9) /= w_sum;
        amrex::ParticleReal const pxpt   = values_per_rank_2nd.at(10) /= w_sum;
        amrex::ParticleReal const ypt    = values_per_rank_2nd.at(11) /= w_sum;
        amrex::ParticleReal const pypt   = values_per_rank_2nd.at(12) /= w_sum;
        amrex::ParticleReal const charge = values_per_rank_2nd.at(13);
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
        // Dispersion and dispersive beam moments
        amrex::ParticleReal const dispersion_x = - xpt / pt_ms;
        amrex::ParticleReal const dispersion_px = - pxpt / pt_ms;
        amrex::ParticleReal const dispersion_y = - ypt / pt_ms;
        amrex::ParticleReal const dispersion_py = - pypt / pt_ms;
        amrex::ParticleReal const x_msd = x_ms - pt_ms*dispersion_x*dispersion_x;
        amrex::ParticleReal const px_msd = px_ms - pt_ms*dispersion_px*dispersion_px;
        amrex::ParticleReal const xpx_d = xpx - pt_ms*dispersion_x*dispersion_px;
        amrex::ParticleReal const emittance_xd = std::sqrt(x_msd*px_msd-xpx_d*xpx_d);
        amrex::ParticleReal const y_msd = y_ms - pt_ms*dispersion_y*dispersion_y;
        amrex::ParticleReal const py_msd = py_ms - pt_ms*dispersion_py*dispersion_py;
        amrex::ParticleReal const ypy_d = ypy - pt_ms*dispersion_y*dispersion_py;
        amrex::ParticleReal const emittance_yd = std::sqrt(y_msd*py_msd-ypy_d*ypy_d);
        // Courant-Snyder (Twiss) beta-function
        amrex::ParticleReal const beta_x = x_msd / emittance_xd;
        amrex::ParticleReal const beta_y = y_msd / emittance_yd;
        amrex::ParticleReal const beta_t = t_ms / emittance_t;
        // Courant-Snyder (Twiss) alpha
        amrex::ParticleReal const alpha_x = - xpx_d / emittance_xd;
        amrex::ParticleReal const alpha_y = - ypy_d / emittance_yd;
        amrex::ParticleReal const alpha_t = - tpt / emittance_t;

        std::unordered_map<std::string, amrex::ParticleReal> data;
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
        data["dispersion_x"] = dispersion_x;
        data["dispersion_px"] = dispersion_px;
        data["dispersion_y"] = dispersion_y;
        data["dispersion_py"] = dispersion_py;
        data["charge_C"] = charge;

        return data;
    }
} // namespace impactx::diagnostics
