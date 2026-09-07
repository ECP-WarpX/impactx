/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Ji Qiang
 * License: BSD-3-Clause-LBNL
 */
#include "Gauss2p5dPush.H"

#include "diagnostics/ReducedBeamCharacteristics.H"
#include "particles/wakefields/ChargeBinning.H"

#include <AMReX_REAL.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_GpuParallelReduce.H>
#include <AMReX_ParmParse.H>
#include <ablastr/warn_manager/WarnManager.H>

#include <cmath>
#include <string>


namespace impactx::particles::spacecharge
{

    /** Calculate the Integrals for Space Charge Fields
     *
     * Compute integrals Eqs 25, 31,32 used in the space-charge fields from a 2D transverse Gaussian distribution (including 2A in Eq. 25).
     * Input particle locations (x,y) and RMS sizes (sigx,sigy) and return the integrals for SC fields.
     */
    AMREX_GPU_DEVICE
    void potInt (
        amrex::ParticleReal delta,
        int nint,
        amrex::ParticleReal xin,
        amrex::ParticleReal yin,
        amrex::ParticleReal sigx,
        amrex::ParticleReal sigy,
        amrex::ParticleReal& pintex,
        amrex::ParticleReal& pintey,
        amrex::ParticleReal& pintez
    )
    {
        using namespace amrex::literals;

        // enforce odd grid points
        if (nint % 2 == 0) nint += 1;

        amrex::ParticleReal xmin = delta;
        amrex::ParticleReal xmax = 1.0;
        amrex::ParticleReal h = (xmax - xmin) / (nint - 1);

        amrex::ParticleReal xp = xin / sigx;
        amrex::ParticleReal yp = yin / sigy;
        amrex::ParticleReal asp = sigx / sigy;
        amrex::ParticleReal xp2 = xp*xp;
        amrex::ParticleReal yp2 = yp*yp;

        // Trapezoidal rule approximation integral from 0 to delta
        amrex::ParticleReal x = xmin;
        amrex::ParticleReal t2 = asp*asp + (1_prt - asp*asp)*x*x;
        amrex::ParticleReal sqrt2 = std::sqrt(t2);
        amrex::ParticleReal exparg = -(xp*xp + yp*yp / t2) / 2_prt;
        amrex::ParticleReal f1 = x*exparg;

        amrex::ParticleReal sum0 = x * f1 / sqrt2 / 2_prt;
        f1 = x*std::exp(x*x*exparg);
        amrex::ParticleReal sum0ex = x * f1 / sqrt2 / 2_prt;
        amrex::ParticleReal sum0ey = x*f1 / (t2*sqrt2) / 2_prt;

        // Simpson's rule for two end points
        // Ez
        f1 = std::exp(x*x*exparg);
        amrex::ParticleReal f2 = x*sqrt2;
        amrex::ParticleReal fx = (f1 - 1_prt) / f2;
        sum0 = sum0 - fx*h;

        amrex::ParticleReal f1x = x*f1;
        f2 = sqrt2;
        fx = f1x / f2;
        sum0ex = sum0ex - fx*h;

        f2 = t2*sqrt2;
        fx = f1x / f2;
        sum0ey = sum0ey - fx*h;

        // Ez at xmax
        x = xmax;
        t2 = asp*asp + (1_prt - asp*asp)*x*x;
        sqrt2 = std::sqrt(t2);
        exparg = -(xp*xp + yp*yp / t2) / 2_prt;
        f1 = std::exp(x*x*exparg);
        f2 = x*sqrt2;
        fx = (f1 - 1_prt) / f2;
        sum0 = sum0 - fx*h;
        f1x = x*f1;
        f2 = sqrt2;
        fx = f1x / f2;
        sum0ex = sum0ex - fx*h;
        f2 = t2*sqrt2;
        fx = f1x / f2;
        sum0ey = sum0ey - fx*h;

        amrex::ParticleReal fy = 0.0;
        amrex::ParticleReal fz = 0.0;
        amrex::ParticleReal t2sqrt = 0.0;
        amrex::ParticleReal f2y = 0.0;
        amrex::ParticleReal asp2 = asp*asp;
        amrex::ParticleReal asp2m = 1_prt - asp*asp;
        for (int i = 0; i < nint; ++i)
        {
            x = xmin + i*h;
            t2 = asp2 + asp2m * x*x;
            t2sqrt = std::sqrt(t2);
            f2y = t2*t2sqrt;
            f2 = x*t2sqrt;
            f1 = std::exp(-x*x*(xp2 + yp2 / t2) / 2_prt);
            f1x = x*f1;

            // Ez
            fz = (f1 - 1_prt) / f2;
            if (i%2 == 0)
              sum0 += 2_prt * fz * h;
            else
              sum0 += 4_prt * fz * h;

            // Ex
            fx = f1x / t2sqrt;
            if (i%2 == 0)
              sum0ex += 2_prt * fx * h;
            else
              sum0ex += 4_prt * fx * h;

            // Ey
            fy = f1x / f2y;
            if (i%2 == 0)
              sum0ey += 2_prt * fy * h;
            else
              sum0ey += 4_prt * fy * h;
        }

        pintez = 2_prt * asp * sum0 / 3_prt;
        pintex = 2_prt * asp * xp * sum0ex / 3_prt / sigx;
        pintey = 2_prt * asp * yp * sum0ey / 3_prt / sigy;
    }

    void Gauss2p5dPush (
        ImpactXParticleContainer & pc,
        amrex::ParticleReal const slice_ds
    )
    {
        BL_PROFILE("impactx::spacecharge::Gauss2p5dPush");

        using namespace amrex::literals;

        amrex::ParticleReal const charge = pc.GetRefParticle().charge;

        std::unordered_map<std::string, amrex::ParticleReal> const rbc = diagnostics::reduced_beam_characteristics(pc);
        // get RMS sigma of beam
        // Standard deviation of the particle positions (speed of light times time delay for t) (unit: meter)
        amrex::ParticleReal const sigx = rbc.at("sig_x");
        amrex::ParticleReal const sigy = rbc.at("sig_y");
        amrex::ParticleReal const sigt = rbc.at("sig_t");

        // physical constants and reference quantities
        amrex::ParticleReal const c0_SI = 2.99792458e8_prt;  // TODO move out
        amrex::ParticleReal const mc_SI = pc.GetRefParticle().mass * c0_SI;
        amrex::ParticleReal const pz_ref_SI = pc.GetRefParticle().beta_gamma() * mc_SI;
        amrex::ParticleReal const gamma = pc.GetRefParticle().gamma();
        amrex::ParticleReal const beta = pc.GetRefParticle().beta();
        amrex::ParticleReal const beta_gamma = pc.GetRefParticle().beta_gamma();
        amrex::ParticleReal const inv_gamma2 = 1.0_prt / (gamma * gamma);
        amrex::ParticleReal const rfpiepslon = c0_SI * c0_SI * 1.0e-7_prt;

        amrex::ParticleReal const dt = slice_ds / beta / c0_SI;
        amrex::ParticleReal const aspect_ratio = std::sqrt(sigx*sigx + sigy*sigy) / (beta_gamma * sigt);

        if (aspect_ratio > 1_rt) {
            ablastr::warn_manager::WMRecordWarning(
                "Impactx::particles::spacecharge::Gauss2p5dPush",
                "Gauss2p5D model assumes a long bunch, but here "
                "sigma_perp / (beta_gamma * sigmaz) = " +
                std::to_string(aspect_ratio) + " > 1.",
                ablastr::warn_manager::WarnPriority::medium
            );
        }

        int nint = 101;
        amrex::Real delta = 0.01_rt;
        // note: the default value below is the optimal value minimizing the L2-norm
        // of the error in the on-axis longitudinal field Ez for a 3D Gaussian bunch
        amrex::Real long_scale = 1.103_rt * beta_gamma * sigt;
        amrex::ParmParse pp_algo("algo.space_charge");
        pp_algo.queryAddWithParser("gauss_nint", nint);
        pp_algo.queryAddWithParser("gauss_taylor_delta", delta);
        // note: intentionall w/o add because `sigt` is dynamic!
        //       add would ignore the new beam size in later sim steps
        pp_algo.queryWithParser("gauss_long_scale", long_scale);

        int tp5d_bins = 129;
        pp_algo.queryAddWithParser("gauss_charge_z_bins", tp5d_bins);

        bool apply_longitudinal_kick = true;
        pp_algo.queryAdd("apply_longitudinal_kick", apply_longitudinal_kick);

        // Measure beam size, extract the min, max of particle positions
        [[maybe_unused]] auto const [x_min, y_min, t_min, x_max, y_max, t_max] =
            pc.MinAndMaxPositions();

        // Set parameters for charge deposition
        bool const is_unity_particle_weight = false;
        bool const GetNumberDensity = true;

        int const num_bins = tp5d_bins;  // Set resolution
        amrex::Real const bin_min = t_min;
        amrex::Real const bin_max = t_max;
        amrex::Real const bin_size = (bin_max - bin_min) / (num_bins - 1);  // number of evaluation points
        // Allocate memory for the charge profile
        amrex::Gpu::DeviceVector<amrex::Real> charge_distribution(num_bins + 1, 0.0);
        // Call charge deposition function
        impactx::particles::wakefields::DepositCharge1D(pc, charge_distribution, bin_min, bin_size, is_unity_particle_weight);

        // Sum up all partial charge histograms to each MPI process to calculate
        // the global charge slope.
        amrex::ParallelAllReduce::Sum(
            charge_distribution,
            amrex::ParallelDescriptor::Communicator()
        );

        // Call charge density derivative function
        amrex::Gpu::DeviceVector<amrex::Real> slopes(charge_distribution.size() - 1, 0.0);
        impactx::particles::wakefields::DerivativeCharge1D(charge_distribution, slopes, bin_size,GetNumberDensity);

        amrex::Real const * const beam_profile_slope = slopes.data();
        amrex::Real const * const beam_profile = charge_distribution.data();

        // group together constants for the momentum push
        amrex::ParticleReal const push_consts = rfpiepslon * dt * charge * inv_gamma2 / (beta * pz_ref_SI);
        amrex::ParticleReal const chargesign = charge / std::abs(charge);
        amrex::ParticleReal const log2n = -std::log(2.0_prt);
        amrex::ParticleReal const pz_push_const =
            log2n
            + 0.577216_prt
            - 2.0_prt * std::log((sigx + sigy)/long_scale/2.0_prt);

        // loop over refinement levels
        int const nLevel = pc.finestLevel();
        for (int lev = 0; lev <= nLevel; ++lev)
        {
            // loop over all particle boxes
            using ParIt = ImpactXParticleContainer::iterator;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (ParIt pti(pc, lev); pti.isValid(); ++pti)
            {
                const int np = pti.numParticles();

                // preparing access to particle data: SoA of Reals
                auto& soa_real = pti.GetStructOfArrays().GetRealData();
                amrex::ParticleReal* const AMREX_RESTRICT part_x = soa_real[RealSoA::x].dataPtr();
                amrex::ParticleReal* const AMREX_RESTRICT part_y = soa_real[RealSoA::y].dataPtr();
                amrex::ParticleReal* const AMREX_RESTRICT part_z = soa_real[RealSoA::z].dataPtr(); // note: currently for a fixed t
                amrex::ParticleReal* const AMREX_RESTRICT part_px = soa_real[RealSoA::px].dataPtr();
                amrex::ParticleReal* const AMREX_RESTRICT part_py = soa_real[RealSoA::py].dataPtr();
                amrex::ParticleReal* const AMREX_RESTRICT part_pz = soa_real[RealSoA::pz].dataPtr(); // note: currently for a fixed t

                // gather to each particle and push momentum
                amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) {
                    // access SoA Real data
                    amrex::ParticleReal & AMREX_RESTRICT x = part_x[i];
                    amrex::ParticleReal & AMREX_RESTRICT y = part_y[i];
                    amrex::ParticleReal & AMREX_RESTRICT z = part_z[i];
                    amrex::ParticleReal & AMREX_RESTRICT px = part_px[i];
                    amrex::ParticleReal & AMREX_RESTRICT py = part_py[i];
                    amrex::ParticleReal & AMREX_RESTRICT pz = part_pz[i];

                    // field integrals from a 2D Gaussian bunch
                    amrex::ParticleReal eintx, einty, eintz;
                    potInt(delta,nint,x,y,sigx,sigy,eintx,einty,eintz);

                    // Update momentae with the 2.5D SC force
                    int const idx = static_cast<int>((z - bin_min) / bin_size);  // Find index position along z
#if (defined(AMREX_DEBUG) || defined(DEBUG)) && !defined(AMREX_USE_GPU)
                    if (idx < 0 || idx >= num_bins)
                    {
                        std::cerr << "Warning: Index out of range for 2.5D Gaussian SC: " << idx << std::endl;
                    }
#endif
                    // Interpolation
                    amrex::ParticleReal const idxreal = static_cast<amrex::ParticleReal>(idx);
                    amrex::ParticleReal const lambda = (z - bin_min) / bin_size - idxreal - 0.5;
                    amrex::ParticleReal const beam_profile_interp = (1_prt-lambda)*beam_profile[idx] + lambda*beam_profile[idx+1];
                    amrex::ParticleReal const beam_profile_slope_interp = (idx==num_bins)? beam_profile_slope[idx] : (1_prt-lambda)*beam_profile_slope[idx] + lambda*beam_profile_slope[idx+1];


                    // Force evaluation
                    amrex::ParticleReal const Fxy = beam_profile_interp * chargesign;
                    amrex::ParticleReal const Fz = (apply_longitudinal_kick)? beam_profile_slope_interp * charge : 0_prt;

                    // push momentum
                    px += eintx * Fxy * push_consts;
                    py += einty * Fxy * push_consts;
                    pz -= (eintz + pz_push_const) * Fz * push_consts;

                    // push position is done in the lattice elements
                });
            } // end loop over all particle boxes
        } // env mesh-refinement level loop
    }

}  // namespace impactx::particles::spacecharge
