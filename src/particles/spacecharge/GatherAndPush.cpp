/* Copyright 2022 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Marco Garten, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "GatherAndPush.H"

#include <ablastr/particles/NodalFieldGather.H>

#include <AMReX_BLProfiler.H>
#include <AMReX_REAL.H>       // for Real
#include <AMReX_SPACE.H>      // for AMREX_D_DECL


namespace impactx::spacecharge
{
    void GatherAndPush (
            ImpactXParticleContainer & pc,
            std::unordered_map<int, std::unordered_map<std::string, amrex::MultiFab> > const & space_charge_field,
            const amrex::Vector<amrex::Geometry>& geom,
            amrex::ParticleReal const slice_ds
    )
    {
        BL_PROFILE("impactx::spacecharge::GatherAndPush");

        using namespace amrex::literals;

        amrex::ParticleReal const charge = pc.GetRefParticle().charge;

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
            for (ParIt pti(pc, lev); pti.isValid(); ++pti) {
                const int np = pti.numParticles();

                // get the device pointer-wrapper Array4 for 3D field access
                auto const scf_arr_x = space_charge_field.at(lev).at("x")[pti].array();
                auto const scf_arr_y = space_charge_field.at(lev).at("y")[pti].array();
                auto const scf_arr_z = space_charge_field.at(lev).at("z")[pti].array();

                // physical constants and reference quantities
                amrex::ParticleReal const c0_SI = 2.99792458e8;  // TODO move out
                amrex::ParticleReal const mc_SI = pc.GetRefParticle().mass * c0_SI;
                amrex::ParticleReal const pz_ref_SI = pc.GetRefParticle().beta_gamma() * mc_SI;
                amrex::ParticleReal const gamma = pc.GetRefParticle().gamma();
                amrex::ParticleReal const inv_gamma2 = 1.0_prt / (gamma * gamma);

                amrex::ParticleReal const dt = slice_ds / pc.GetRefParticle().beta() / c0_SI;

                // preparing access to particle data: AoS
                using PType = ImpactXParticleContainer::ParticleType;
                auto const & aos = pti.GetArrayOfStructs();
                PType const * const AMREX_RESTRICT aos_ptr = aos().dataPtr();

                // preparing access to particle data: SoA of Reals
                auto& soa_real = pti.GetStructOfArrays().GetRealData();
                amrex::ParticleReal* const AMREX_RESTRICT part_px = soa_real[RealSoA::ux].dataPtr();
                amrex::ParticleReal* const AMREX_RESTRICT part_py = soa_real[RealSoA::uy].dataPtr();
                amrex::ParticleReal* const AMREX_RESTRICT part_pz = soa_real[RealSoA::pt].dataPtr(); // note: currently in z

                // group together constants for the momentum push
                amrex::ParticleReal const push_consts = dt * charge * inv_gamma2 / pz_ref_SI;

                // gather to each particle and push momentum
                amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) {
                    // access AoS data such as positions and cpu/id
                    PType const & AMREX_RESTRICT p = aos_ptr[i];

                    // access SoA Real data
                    amrex::ParticleReal & AMREX_RESTRICT px = part_px[i];
                    amrex::ParticleReal & AMREX_RESTRICT py = part_py[i];
                    amrex::ParticleReal & AMREX_RESTRICT pz = part_pz[i];

                    // force gather
                    amrex::GpuArray<amrex::Real, 3> const field_interp =
                        ablastr::particles::doGatherVectorFieldNodal (
                            p.pos(0), p.pos(1), p.pos(2),
                            scf_arr_x, scf_arr_y, scf_arr_z,
                            invdr,
                            prob_lo);

                    // push momentum
                    px += field_interp[0] * push_consts;
                    py += field_interp[1] * push_consts;
                    pz += field_interp[2] * push_consts;

                    // push position is done in the lattice elements
                });


            } // end loop over all particle boxes
        } // env mesh-refinement level loop
    }
} // namespace impactx::spacecharge
