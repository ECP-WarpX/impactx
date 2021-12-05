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
    void Push (ImpactXParticleContainer & pc,
               std::list<KnownElements> const & beamline_elements)
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
                const int np = pti.numParticles();
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

                // loop over all beamline elements
                for (auto & element_variant : beamline_elements) {
                    // here we just access the element by its respective type
                    std::visit([=](auto&& element) {
                        // loop over particles in the box
                        amrex::ParallelFor( np, [=] AMREX_GPU_DEVICE (long i)
                        {
                            // access AoS data such as positions and cpu/id
                            PType& p = aos_ptr[i];

                            // access SoA Real data
                            amrex::ParticleReal & px = part_px[i];
                            amrex::ParticleReal & py = part_py[i];
                            amrex::ParticleReal & pt = part_pt[i];

                            element(p, px, py, pt);
                        });
                    }, element_variant);
                }; // end loop over all beamline elements

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
