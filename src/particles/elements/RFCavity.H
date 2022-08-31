/* Copyright 2022 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_RFCAVITY_H
#define IMPACTX_RFCAVITY_H

#include "particles/ImpactXParticleContainer.H"
#include "particles/integrators/Integrators.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>

#include <cmath>


namespace impactx
{
    struct RFCavity
    {
        static constexpr auto name = "RFCavity";
        using PType = ImpactXParticleContainer::ParticleType;

        /** An RF cavity (with solenoid field)
         *
         * @param ds Segment length in m
         * @param Escale Scaling factor for on-axis RF electric field
         * @param freq RF frequency in Hz
         * @param phase RF driven phase in deg
         * @param Bzsol Strength of solenoid field in T
         * @param nslice number of slices used for the application of space charge
         */
        RFCavity( amrex::ParticleReal const ds, amrex::ParticleReal const
         Escale, amrex::ParticleReal const freq, amrex::ParticleReal const
         phase, amrex::ParticleReal const Bzsol, int const nslice )
        : m_ds(ds), m_Escale(Escale), m_freq(freq), m_phase(phase), m_Bzsol(Bzsol), m_nslice(nslice)
        {
        }

        /** This is an RF cavity functor, so that a variable of this type can be used like
         *  an RF cavity function.
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

            // initialize output values of momenta
            amrex::ParticleReal pxout = px;
            amrex::ParticleReal pyout = py;
            amrex::ParticleReal ptout = pt;

            // length of the current slice
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // access reference particle values to find beta*gamma^2
            amrex::ParticleReal const pt_ref = refpart.pt;
            amrex::ParticleReal const betgam2 = pow(pt_ref, 2) - 1.0_prt;

            // The following map is that of a drift.  It must be
            // replaced by the map for an RF cavity.

            // advance position and momentum (drift)
            p.pos(0) = x + slice_ds * px;
            // pxout = px;
            p.pos(1) = y + slice_ds * py;
            // pyout = py;
            p.pos(2) = t + (slice_ds/betgam2) * pt;
            // ptout = pt;

            // assign updated momenta
            px = pxout;
            py = pyout;
            pt = ptout;

        }

        /** This pushes the reference particle.
         *
         * @param[in,out] refpart reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
                RefPart & AMREX_RESTRICT refpart) const {

            using namespace amrex::literals; // for _rt and _prt

            // assign input reference particle values
            amrex::ParticleReal const x = refpart.x;
            amrex::ParticleReal const px = refpart.px;
            amrex::ParticleReal const y = refpart.y;
            amrex::ParticleReal const py = refpart.py;
            amrex::ParticleReal const z = refpart.z;
            amrex::ParticleReal const pz = refpart.pz;
            amrex::ParticleReal const pt = refpart.pt;
            amrex::ParticleReal const s = refpart.s;

            // length of the current slice
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // compute intial value of beta*gamma
            amrex::ParticleReal const bgi = sqrt(pow(pt,2)-1.0_prt);

            // call integrator to advance (t,pt)
            amrex::ParticleReal const zin = s;
            amrex::ParticleReal const zout = s + slice_ds;
            std::string element_type = "rfcavity";
            int const nsteps = 10;

            integrators::symp4_integrate(refpart,zin,zout,nsteps,element_type);
            amrex::ParticleReal const ptf = refpart.pt;

            // advance position (x,y,z)
            refpart.x = x + slice_ds*px/bgi;
            refpart.y = y + slice_ds*py/bgi;
            refpart.z = z + slice_ds*pz/bgi;

            // compute final value of beta*gamma
            amrex::ParticleReal const bgf = sqrt(pow(ptf,2)-1.0_prt);

            // advance momentum (px,py,pz)
            refpart.px = px*bgf/bgi;
            refpart.py = py*bgf/bgi;
            refpart.pz = pz*bgf/bgi;

            // advance integrated path length
            refpart.s = s + slice_ds;

        }

        /** Number of slices used for the application of space charge
         *
         * @return positive integer
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        int nslice () const
        {
            return m_nslice;
        }

    private:
        amrex::ParticleReal m_ds; //! segment length in m
        amrex::ParticleReal m_Escale; //! scaling factor for RF electric field
        amrex::ParticleReal m_freq; //! RF frequency in Hz
        amrex::ParticleReal m_phase; //! RF driven phase in deg
        amrex::ParticleReal m_Bzsol; //! strength of solenoid field in T
        int m_nslice; //! number of slices used for the application of space charge
    };

} // namespace impactx

#endif // IMPACTX_RFCAVITY_H
