/* Copyright 2021 Axel Huebl, Chad Mitchell
 *
 * This file is part of ImpactX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "DiagnosticOutput.H"

#include <AMReX_Extension.H>  // for AMREX_RESTRICT
#include <AMReX_REAL.H>       // for ParticleReal
#include <AMReX_Print.H>      // for PrintToFile

namespace impactx
{
namespace diagnostics {
    void DiagnosticOutput(ImpactXParticleContainer &pc,
                                  OutputType const &otype) {
        using namespace amrex::literals; // for _rt and _prt

        // create a host-side particle buffer
        auto tmp = pc.make_alike<amrex::PinnedArenaAllocator>();

        // copy device-to-host
        bool const local = true;
        tmp.copyParticles(pc, local);

        // loop over refinement levels
        int const nLevel = tmp.maxLevel() + 1;
        for (int lev = 0; lev < nLevel; ++lev)
        {
            // loop over all particle boxes
            using ParIt = typename decltype(tmp)::ParIterType; // you can even try to use ::ParConstIterType
            for (ParIt pti(tmp, lev); pti.isValid(); ++pti) {
                const int np = pti.numParticles();

               // preparing access to particle data: AoS
                using PType = ImpactXParticleContainer::ParticleType;
                auto &aos = pti.GetArrayOfStructs();
                PType *AMREX_RESTRICT aos_ptr = aos().dataPtr();

                // preparing access to particle data: SoA of Reals
                auto &soa_real = pti.GetStructOfArrays().GetRealData();
                amrex::ParticleReal *const AMREX_RESTRICT part_px = soa_real[RealSoA::ux].dataPtr();
                amrex::ParticleReal *const AMREX_RESTRICT part_py = soa_real[RealSoA::uy].dataPtr();
                amrex::ParticleReal *const AMREX_RESTRICT part_pt = soa_real[RealSoA::pt].dataPtr();


                if( otype == OutputType::PrintParticles) {
                // print out particles (this hack works only on CPU and on GPUs with
                // unified memory access)
                for (int i=0; i < np; ++i)
                {

                        // access AoS data such as positions and cpu/id
                        PType const& p = aos_ptr[i];
                        amrex::ParticleReal const x = p.pos(0);
                        amrex::ParticleReal const y = p.pos(1);
                        amrex::ParticleReal const t = p.pos(2);

                        // access SoA Real data
                        amrex::ParticleReal &px = part_px[i];
                        amrex::ParticleReal &py = part_py[i];
                        amrex::ParticleReal &pt = part_pt[i];

                        // write particle data to file
                        amrex::PrintToFile("output_beam.txt") << x << " " << y << " ";
                        amrex::PrintToFile("output_beam.txt") << t << " " << px << " ";
                        amrex::PrintToFile("output_beam.txt") << py << " " << pt << " " << std::endl;

                };
            } // end loop over all particle boxes
         }
        } // env mesh-refinement level loop

    }
} // namespace diagnostics
} // namespace impactx
