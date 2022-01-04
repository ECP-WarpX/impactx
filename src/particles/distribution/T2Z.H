/* Copyright 2021 Chad Mitchell, Axel Huebl
 *
 * This file is part of ImpactX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_T2Z_H
#define IMPACTX_T2Z_H

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_REAL.H>

#include <cmath>


namespace impactx
{
    struct T2Z
    {
        using PType = ImpactXParticleContainer::ParticleType;

        /** A T-to-Z transformation
         *
         * @param pzd: Design value of pz/mc = beta*gamma.
         */
        T2Z( amrex::ParticleReal const pzd )
        : m_pzd(pzd)
        {
        }

        /** This is a z-to-t map, so that a variable of this type can be used like a z-to-t function.
         *
         * @param p Particle AoS data for positions and cpu/id
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
                PType& p,
                amrex::ParticleReal & px,
                amrex::ParticleReal & py,
                amrex::ParticleReal & pt) const {

            using namespace amrex::literals;

            // access AoS data such as positions and cpu/id
            amrex::ParticleReal const x = p.pos(0);
            amrex::ParticleReal const y = p.pos(1);
            amrex::ParticleReal const t = p.pos(2);

            // intermediate values (universal for all elements - needs tracking for RF elements)
            amrex::ParticleReal const ptf;
            amrex::ParticleReal const ptdf;
            amrex::ParticleReal const arg;
            amrex::ParticleReal const argd;

            // arguments for ptd and pt
            argd = 1.0_prt + pow(pzd, 2);
            arg = 1.0_prt + pow(pzd+pt, 2) + pow(px, 2) + pow(py, 2);

            // check for validity of arguments
            if (argd > 0.0_prt) {
               ptdf = -sqrt(argd);
            } else {
               cout << "invalid ptd arg: " << argd << endl;
               ptdf = -1.0_prt;
            }
            if (arg > 0.0_prt) {
               ptf = -sqrt(arg);
            } else {
               cout << "invalid pt arg: " << arg << endl;
               ptf = -1.0_prt;
            }

           // transform position and momentum (z2t)
            p.pos(0) = x - px*t/(m_pzd+pt);
            // px = px;
            p.pos(1) = y - py*t/(m_pzd+pt);
            // py = py;
            p.pos(2) = ptf*t/(m_pzd+pt);  // This now represents t.
            pt = ptf - ptd;   // This now represents pt.

        }

    private:
        amrex::ParticleReal m_pzd;
    };

} // namespace impactx

#endif // IMPACTX_T2Z_H
