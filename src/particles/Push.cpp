/* Copyright 2021 Axel Huebl
 *
 * This file is part of ImpactX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Push.H"

#include "elements/Drift.H"
#include "elements/Quad.H"

#include <AMReX_Extension.H>  // for AMREX_RESTRICT
#include <AMReX_REAL.H>       // for ParticleReal


namespace impactx
{
    void Push (ImpactXParticleContainer & pc)
    {
        using namespace amrex::literals; // for _rt and _prt

        // loop over refinement levels
        int const nLevel = 1;
        for (int lev = 0; lev < nLevel; ++lev)
        {
            // get simulation geometry information
            //const amrex::Geometry& gm = this->Geom(lev);
            //const auto prob_lo = gm.ProbLo();

            // loop over all particle boxes
            using ParIt = ImpactXParticleContainer::iterator;
            for (ParIt pti(pc, lev); pti.isValid(); ++pti) {
                //const auto t_lev = pti.GetLevel();
                //const auto index = pti.GetPairIndex();
                // ...

                // preparing access to particle data: AoS
                using PType = ImpactXParticleContainer::ParticleType;
                auto& aos = pti.GetArrayOfStructs();
                PType* AMREX_RESTRICT aos_ptr = aos().dataPtr();

                // preparing access to particle data: SoA of Reals
                auto& soa_real = pti.GetStructOfArrays().GetRealData();
                amrex::ParticleReal* const AMREX_RESTRICT part_px = soa_real[RealSoA::ux].dataPtr();
                amrex::ParticleReal* const AMREX_RESTRICT part_py = soa_real[RealSoA::uy].dataPtr();
                amrex::ParticleReal* const AMREX_RESTRICT part_pt = soa_real[RealSoA::pt].dataPtr();
                // ...

                amrex::ParticleReal const ds = 0.1; // Segment length in m.
                amrex::ParticleReal const k = 1.0; // quadrupole strength in 1/m

                // loop over particles in the box
                const int np = pti.numParticles();
                amrex::ParallelFor( np, [=] AMREX_GPU_DEVICE (long i)
                {
                    // access AoS data such as positions and cpu/id
                    PType& p = aos_ptr[i];

                    // access SoA Real data
                    amrex::ParticleReal & px = part_px[i];
                    amrex::ParticleReal & py = part_py[i];
                    amrex::ParticleReal & pt = part_pt[i];

                    Drift drift(ds);
                    drift(p, px, py, pt);

                    Quad quad(ds, k);
                    quad(p, px, py, pt);
                });

                // print out particles (this hack works only on CPU and on GPUs with
                // unified memory access)
                for (int i=0; i < np; ++i)
                {
                    // access AoS data such as positions and cpu/id
                    PType const& p = aos_ptr[i];
                    auto const id = p.id();
                    auto const cpu = p.cpu();
                    auto const pos = p.pos();

                    amrex::AllPrint()
                              << "Particle created at rank=" << cpu
                              << " (pid=" << id << ") is now at: "
                              << pos << "\n";
                };
            } // end loop over all particle boxes
        } // env mesh-refinement level loop
    }

} // namespace impactx
