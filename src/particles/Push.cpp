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
                amrex::ParticleReal* const AMREX_RESTRICT part_ux = soa_real[RealSoA::ux].dataPtr();
                // ...

                // loop over particles in the box
                const int np = pti.numParticles();
                amrex::ParallelFor( np, [=] AMREX_GPU_DEVICE (long i)
                {
                    // access AoS data such as positions and cpu/id
                    PType& p = aos_ptr[i];
                    amrex::ParticleReal const x = p.pos(0);
                    //amrex::ParticleReal const y = p.pos(1);
                    //amrex::ParticleReal const z = p.pos(2);

                    // acces SoA Real data
                    amrex::ParticleReal const ux = part_ux[i];

                    // advance position
                    amrex::ParticleReal const dt = 1.0_prt;
                    amrex::ParticleReal const F_x = 0.1_prt;
                    p.pos(0) = x + ux * dt + F_x;
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
