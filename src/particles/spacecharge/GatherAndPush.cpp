/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Marco Garten, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "GatherAndPush.H"

#include "initialization/Algorithms.H"
#include "particles/wakefields/ChargeBinning.H"
#include "particles/spacecharge/Deposit1D.H"

#include <ablastr/particles/NodalFieldGather.H>

#include <AMReX_BLProfiler.H>
#include <AMReX_REAL.H>       // for Real
#include <AMReX_SPACE.H>      // for AMREX_D_DECL
#include <AMReX_GpuContainers.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Reduce.H>

namespace impactx::particles::spacecharge
{
    void GatherAndPush (
        ImpactXParticleContainer & pc,
        std::unordered_map<int, std::unordered_map<std::string, amrex::MultiFab> > const & space_charge_field,
        std::unordered_map<int, amrex::MultiFab> const & space_charge_potential,
        const amrex::Vector<amrex::Geometry>& geom,
        amrex::ParticleReal const slice_ds
    )
    {
        BL_PROFILE("impactx::spacecharge::GatherAndPush");

        using namespace amrex::literals;

        auto space_charge = get_space_charge_algo();

        amrex::ParticleReal const charge = pc.GetRefParticle().charge;

        // Deposit 1D charge density in cases where it is required.
        int num_bins = 100;
        amrex::ParmParse pp_algo("algo.space_charge");
        pp_algo.queryAddWithParser("num_longitudinal_bins", num_bins);
        amrex::Gpu::DeviceVector<amrex::Real> charge_distribution(num_bins + 1, 0.0);
        amrex::Gpu::DeviceVector<amrex::Real> charge_distribution_slope(num_bins, 0.0);
        amrex::Real Qb_abs = 0.0;

        bool apply_longitudinal_kick = true;

        [[maybe_unused]] auto const [x_min, y_min, t_min, x_max, y_max, t_max] =
            pc.MinAndMaxPositions();

        amrex::Real bin_min = t_min;
        amrex::Real bin_max = t_max;
        amrex::Real const bin_size = (bin_max - bin_min) / (num_bins - 1);

        if (space_charge == SpaceChargeAlgo::True_2p5D) {

            pp_algo.queryAdd("apply_longitudinal_kick", apply_longitudinal_kick);
            charge_distribution = Deposit1D( pc, bin_min, bin_max, num_bins);
            bool const GetNumberDensity = true;
            impactx::particles::wakefields::DerivativeCharge1D(charge_distribution, charge_distribution_slope, bin_size, GetNumberDensity);
            Qb_abs = bin_size * amrex::Reduce::Sum(charge_distribution.size(), charge_distribution.data());
        }

        amrex::Real const * const beam_profile = charge_distribution.data();
        amrex::Real const * const beam_profile_slope = charge_distribution_slope.data();

        // loop over refinement levels
        int const nLevel = pc.finestLevel();
        for (int lev = 0; lev <= nLevel; ++lev)
        {
            // get simulation geometry information
            auto const &gm = geom[lev];
            auto const dr = gm.CellSizeArray();
            amrex::GpuArray<amrex::Real, 3> const invdr{AMREX_D_DECL(1_rt/dr[0], 1_rt/dr[1], 1_rt/dr[2])};
            const auto prob_lo = gm.ProbLoArray();

            // loop over all particle boxes
            using ParIt = ImpactXParticleContainer::iterator;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (ParIt pti(pc, lev); pti.isValid(); ++pti) {
                const int np = pti.numParticles();

                // get the device pointer-wrapper Array4 for 3D field access
                auto const scf_arr_x = space_charge_field.at(lev).at("x")[pti].array();
                auto const scf_arr_y = space_charge_field.at(lev).at("y")[pti].array();
                auto const scf_arr_z = space_charge_field.at(lev).at("z")[pti].array();

                // get the device pointer-wrapper Array4 for 3D potential access
                auto const phi_arr = space_charge_potential.at(lev)[pti].const_array();

                // physical constants and reference quantities
                amrex::ParticleReal const c0_SI = 2.99792458e8_prt;  // TODO move out
                amrex::ParticleReal const mc_SI = pc.GetRefParticle().mass * c0_SI;
                amrex::ParticleReal const pz_ref_SI = pc.GetRefParticle().beta_gamma() * mc_SI;
                amrex::ParticleReal const gamma = pc.GetRefParticle().gamma();
                amrex::ParticleReal const beta_gamma = pc.GetRefParticle().beta_gamma();
                amrex::ParticleReal const beta = beta_gamma / gamma;
                amrex::ParticleReal const inv_gamma2 = 1.0_prt / (gamma * gamma);

                amrex::ParticleReal const dt = slice_ds / pc.GetRefParticle().beta() / c0_SI;

                // preparing access to particle data: SoA of Reals
                auto& soa_real = pti.GetStructOfArrays().GetRealData();
                amrex::ParticleReal* const AMREX_RESTRICT part_x = soa_real[RealSoA::x].dataPtr();
                amrex::ParticleReal* const AMREX_RESTRICT part_y = soa_real[RealSoA::y].dataPtr();
                amrex::ParticleReal* const AMREX_RESTRICT part_z = soa_real[RealSoA::z].dataPtr(); // note: currently for a fixed t
                amrex::ParticleReal* const AMREX_RESTRICT part_px = soa_real[RealSoA::px].dataPtr();
                amrex::ParticleReal* const AMREX_RESTRICT part_py = soa_real[RealSoA::py].dataPtr();
                amrex::ParticleReal* const AMREX_RESTRICT part_pz = soa_real[RealSoA::pz].dataPtr(); // note: currently for a fixed t

                // group together constants for the momentum push
                amrex::ParticleReal const push_consts = dt * charge * inv_gamma2 / pz_ref_SI;

                // gather to each particle and push momentum
                if (space_charge == SpaceChargeAlgo::True_2D) {
                    // flatten 3rd dimension
                    auto prob_lo_2D = gm.ProbLoArray();
                    prob_lo_2D[2] = 0.0_rt;

                    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) {
                        // access SoA Real data
                        amrex::ParticleReal & AMREX_RESTRICT x = part_x[i];
                        amrex::ParticleReal & AMREX_RESTRICT y = part_y[i];
                        amrex::ParticleReal z = 0.0_prt;  // flatten 3rd dimension
                        amrex::ParticleReal & AMREX_RESTRICT px = part_px[i];
                        amrex::ParticleReal & AMREX_RESTRICT py = part_py[i];
                        amrex::ParticleReal & AMREX_RESTRICT pz = part_pz[i];

                        // force gather
                        amrex::GpuArray<amrex::Real, 3> const field_interp =
                            ablastr::particles::doGatherVectorFieldNodal<2>(
                                x, y, z,
                                scf_arr_x, scf_arr_y, scf_arr_z,
                                invdr,
                                prob_lo_2D
                            );

                        // push momentum
                        px += field_interp[0] * push_consts * dr[2] / (beta * c0_SI);
                        py += field_interp[1] * push_consts * dr[2] / (beta * c0_SI);
                        pz += 0.0_rt;

                        // push position is done in the lattice elements
                    });
                }
                if (space_charge == SpaceChargeAlgo::True_2p5D) {
                    // flatten 3rd dimension
                    auto prob_lo_2D = gm.ProbLoArray();
                  //  amrex::Array4<const amrex::Real> const phi_arr;
                    prob_lo_2D[2] = 0.0_rt;

                    // group together constants for the momentum push
                    amrex::ParticleReal const charge_abs = std::abs(charge);

                    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) {
                        // access SoA Real data
                        amrex::ParticleReal & AMREX_RESTRICT x = part_x[i];
                        amrex::ParticleReal & AMREX_RESTRICT y = part_y[i];
                        amrex::ParticleReal & AMREX_RESTRICT z = part_z[i];
                        amrex::ParticleReal z_flat = 0.0_prt;  // flatten 3rd dimension
                        amrex::ParticleReal & AMREX_RESTRICT px = part_px[i];
                        amrex::ParticleReal & AMREX_RESTRICT py = part_py[i];
                        amrex::ParticleReal & AMREX_RESTRICT pz = part_pz[i];

                        // force gather
                        amrex::GpuArray<amrex::Real, 3> const field_interp =
                            ablastr::particles::doGatherVectorFieldNodal<2>(
                                x, y, z_flat,
                                scf_arr_x, scf_arr_y, scf_arr_z,
                                invdr,
                                prob_lo_2D
                            );

                        // potential gather
                        amrex::Real potential_interp = 0.0;
                        if (apply_longitudinal_kick) {
                           potential_interp =
                               ablastr::particles::doGatherScalarFieldNodal<2>(
                                  x, y, z_flat,
                                   phi_arr,
                                   invdr,
                                   prob_lo_2D
                               );
                       }

                       // Update momentae with the 2.5D SC force
                       int const idx = static_cast<int>((z - bin_min) / bin_size);  // Find index position along z
                       #if (defined(AMREX_DEBUG) || defined(DEBUG)) && !defined(AMREX_USE_GPU)
                       if (idx < 0 || idx >= num_bins)
                       {
                            std::cerr << "Warning: Index out of range for 2.5D SC: " << idx << std::endl;
                       }
                       #endif

                       // Interpolation
                       amrex::ParticleReal const idxreal = static_cast<amrex::ParticleReal>(idx);
                       amrex::ParticleReal const lambda = (z - bin_min) / bin_size - idxreal;
                       amrex::ParticleReal const beam_profile_interp = (1_prt-lambda)*beam_profile[idx] + lambda*beam_profile[idx+1];
                       amrex::ParticleReal const beam_profile_slope_interp = (idx==num_bins)? beam_profile_slope[idx] : (1_prt-lambda)*beam_profile_slope[idx] + lambda*beam_profile_slope[idx+1];

                       // Force constributions
                       amrex::ParticleReal const Fxy = (Qb_abs==0.0) ? 0.0_prt : beam_profile_interp / Qb_abs;
                       amrex::ParticleReal const Fz = (Qb_abs==0.0) ? 0.0_prt : beam_profile_slope_interp * charge_abs / Qb_abs;

                       // push momentum
                       px += field_interp[0] * Fxy * push_consts * dr[2] / beta;
                       py += field_interp[1] * Fxy * push_consts * dr[2] / beta;
                       pz -= potential_interp * Fz * push_consts * dr[2] / beta;

                    // push position is done in the lattice elements
                    });
                }
                if (space_charge == SpaceChargeAlgo::True_3D) {
                    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) {
                        // access SoA Real data
                        amrex::ParticleReal & AMREX_RESTRICT x = part_x[i];
                        amrex::ParticleReal & AMREX_RESTRICT y = part_y[i];
                        amrex::ParticleReal & AMREX_RESTRICT z = part_z[i];
                        amrex::ParticleReal & AMREX_RESTRICT px = part_px[i];
                        amrex::ParticleReal & AMREX_RESTRICT py = part_py[i];
                        amrex::ParticleReal & AMREX_RESTRICT pz = part_pz[i];

                        // force gather
                        amrex::GpuArray<amrex::Real, 3> const field_interp =
                            ablastr::particles::doGatherVectorFieldNodal(
                                x, y, z,
                                scf_arr_x, scf_arr_y, scf_arr_z,
                                invdr,
                                prob_lo
                            );

                        // push momentum
                        px += field_interp[0] * push_consts;
                        py += field_interp[1] * push_consts;
                        pz += field_interp[2] * push_consts;

                        // push position is done in the lattice elements
                    });
                }
            } // end loop over all particle boxes
        } // env mesh-refinement level loop
    }
}  // namespace impactx::particles::spacecharge
