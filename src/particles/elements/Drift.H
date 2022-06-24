/* Copyright 2022 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_DRIFT_H
#define IMPACTX_DRIFT_H

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>

#include <cmath>


namespace impactx
{
    struct Drift
    {
        using PType = ImpactXParticleContainer::ParticleType;

        /** A drift
         *
         * @param ds Segment length in m.
         */
        Drift( amrex::ParticleReal const ds )
        : m_ds(ds)
        {
        }

        /** This is a drift functor, so that a variable of this type can be used like a drift function.
         *
         * @param p Particle AoS data for positions and cpu/id
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t
         * @param refpart reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
                PType& AMREX_RESTRICT p,
                amrex::ParticleReal & AMREX_RESTRICT px,
                amrex::ParticleReal & AMREX_RESTRICT py,
                amrex::ParticleReal & AMREX_RESTRICT pt,
                RefPart const refpart) const {

            using namespace amrex::literals; // for _rt and _prt

            // access AoS data such as positions and cpu/id
            amrex::ParticleReal const x = p.pos(0);
            amrex::ParticleReal const y = p.pos(1);
            amrex::ParticleReal const t = p.pos(2);

            // intialize output values of momenta
            amrex::ParticleReal pxout = px;
            amrex::ParticleReal pyout = py;
            amrex::ParticleReal ptout = pt;

            // access reference particle values to find beta*gamma^2
            amrex::ParticleReal const pt_ref = refpart.pt;
            amrex::ParticleReal const betgam2 = pow(pt_ref, 2) - 1.0_prt;

            // advance position and momentum (drift)
            p.pos(0) = x + m_ds * px;
            // pxout = px;
            p.pos(1) = y + m_ds * py;
            // pyout = py;
            p.pos(2) = t + (m_ds/betgam2) * pt;
            // ptout = pt;

            // assign updated momenta
            px = pxout;
            py = pyout;
            pt = ptout;

        }

        void operator() (
                RefPart refpart) const {

            using namespace amrex::literals; // for _rt and _prt

            // assign input reference particle values
            amrex::ParticleReal const x = refpart.x;
            amrex::ParticleReal const px = refpart.px;
            amrex::ParticleReal const y = refpart.y;
            amrex::ParticleReal const py = refpart.py;
            amrex::ParticleReal const z = refpart.z;
            amrex::ParticleReal const pz = refpart.pz;
            amrex::ParticleReal const t = refpart.t;
            amrex::ParticleReal const pt = refpart.pt;

            // assign intermediate parameter
            amrex::ParticleReal const s = m_ds/sqrt(pow(pt,2)-1.0_prt);

            // advance position and momentum (drift)
            refpart.x = x + s*px;
            refpart.y = y + s*py;
            refpart.z = z + s*pz;
            refpart.t = t - s*pt;

        }

    private:
        amrex::ParticleReal m_ds;
    };

} // namespace impactx

#endif // IMPACTX_DRIFT_H
