/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */

#include "openPMD.H"
#include "particles/diagnostics/NonlinearLensInvariants.H"
#include "particles/ImpactXParticleContainer.H"

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_REAL.H>
#include <AMReX_RealVect.H>
#include <AMReX_ParmParse.H>

#include <string>
#include <utility>


namespace impactx::diagnostics
{
    void
    add_optional_properties (
        std::string const & element_name,
        ImpactXParticleContainer & pc
    )
    {
        // Parse the diagnostic parameters
        amrex::ParmParse pp_element(element_name);

        bool enabled = false;
        pp_element.queryAdd("nonlinear_lens_invariants", enabled);
        if (!enabled)
            return;

        amrex::ParticleReal alpha = 0.0;
        pp_element.queryAdd("alpha", alpha);

        amrex::ParticleReal beta = 1.0;
        pp_element.queryAdd("beta", beta);

        amrex::ParticleReal tn = 0.4;
        pp_element.queryAdd("tn", tn);

        amrex::ParticleReal cn = 0.01;
        pp_element.queryAdd("cn", cn);

        NonlinearLensInvariants const nonlinear_lens_invariants(alpha, beta, tn, cn);

        // profile time spent here
        std::string profile_name = "impactx::Push::" + std::string(BeamMonitor::type) + "::add_optional_properties";
        BL_PROFILE(profile_name);

        // add runtime properties for H and I
        bool comm = false;
        if (!pc.HasRealComp("H")) {
            pc.AddRealComp("H", comm);
        }
        const int soa_idx_H = pc.GetRealCompIndex("H");
        if (!pc.HasRealComp("I")) {
            pc.AddRealComp("I", comm);
        }
        const int soa_idx_I = pc.GetRealCompIndex("I");

        // undocumented AMReX stuff that does not belong in user code
        // can be removed after these are merged
        //   https://github.com/AMReX-Codes/amrex/pull/3615
        //   https://github.com/AMReX-Codes/amrex/pull/3861
        for (int lev = 0; lev <= pc.finestLevel(); ++lev)
        {
            for (ImpactXParticleContainer::iterator pti(pc, lev); pti.isValid(); ++pti)
            {
                const int grid_id = pti.index();
                const int tile_id = pti.LocalTileIndex();
                pc.DefineAndReturnParticleTile(lev, grid_id, tile_id);
                auto np = pti.numParticles();
                if (np > 0) {
                    auto& soa = pti.GetStructOfArrays();
                    soa.resize(np);
                }
            }
        }

        // loop over refinement levels
        int const nLevel = pc.finestLevel();
        for (int lev = 0; lev <= nLevel; ++lev) {
            // loop over all particle boxes
            using ParIt = ImpactXParticleContainer::iterator;
            for (ParIt pti(pc, lev); pti.isValid(); ++pti) {
                const int np = pti.numParticles();

                // preparing access to particle data: SoA of Reals
                auto & soa = pti.GetStructOfArrays();
                auto const * part_x = soa.GetRealData(RealSoA::x).data();
                auto const * part_y = soa.GetRealData(RealSoA::y).data();
                auto const * part_px = soa.GetRealData(RealSoA::px).data();
                auto const * part_py = soa.GetRealData(RealSoA::py).data();

                auto * part_H = soa.GetRealData(soa_idx_H).data();
                auto * part_I = soa.GetRealData(soa_idx_I).data();

                amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (long i) {
                    amrex::ParticleReal const x = part_x[i];
                    amrex::ParticleReal const y = part_y[i];

                    amrex::ParticleReal const px = part_px[i];
                    amrex::ParticleReal const py = part_py[i];

                    // calculate invariants of motion
                    NonlinearLensInvariants::Data const HI_out =
                        nonlinear_lens_invariants(x, y, px, py);

                    // write particle invariant data
                    part_H[i] = HI_out.H;
                    part_I[i] = HI_out.I;
                });
            } // end loop over all particle boxes
        } // end mesh-refinement level loop
    }

} // namespace impactx::diagnostics
