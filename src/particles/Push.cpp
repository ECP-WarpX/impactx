/* Copyright 2021 Axel Huebl
 *
 * This file is part of ImpactX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Push.H"

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
                
                // loop over particles in the box
                const int np = pti.numParticles();
                amrex::ParallelFor( np, [=] AMREX_GPU_DEVICE (long i)
                {
                    // access AoS data such as positions and cpu/id
                    PType& p = aos_ptr[i];
                    amrex::ParticleReal const x = p.pos(0);
                    amrex::ParticleReal const y = p.pos(1);
                    amrex::ParticleReal const t = p.pos(2);

                    // acces SoA Real data
                    amrex::ParticleReal const px = part_ux[i];
                    amrex::ParticleReal const py = part_uy[i];
                    amrex::ParticleReal const pt = part_pt[i];

                   // intermediate values
                    amrex::ParticleReal const betgam2 = pow(betgam,2);

                    // advance position and momentum (drift)
                    p.pos(0) = x + ds * px;
                    part_ux[i] = px;
                    p.pos(1) = y + ds * py;
                    part_uy[i] = py;
                    p.pos(2) = t + (ds/betgam2) * pt;
                    part_pt[i] = pt;                    

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
