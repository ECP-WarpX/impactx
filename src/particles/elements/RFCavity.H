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
         * @param mapsteps number of integration steps per slice used for
         *        map and reference particle push in applied fields
         * @param nslice number of slices used for the application of space charge
         */
        RFCavity( amrex::ParticleReal const ds, amrex::ParticleReal const
         Escale, amrex::ParticleReal const freq, amrex::ParticleReal const
         phase, amrex::ParticleReal const Bzsol, int const mapsteps, int const nslice )
        : m_ds(ds), m_Escale(Escale), m_freq(freq), m_phase(phase), m_Bzsol(Bzsol), m_mapsteps(mapsteps), m_nslice(nslice)
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
            int const nsteps = m_mapsteps;

            integrators::symp4_integrate(refpart,zin,zout,nsteps,*this);
//            integrators::symp2_integrate(refpart,zin,zout,nsteps,*this);
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

        amrex::ParticleReal
        RF_Efield (amrex::ParticleReal & zeval) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // specify constants
            amrex::ParticleReal const pi = 3.141592653589793_prt;
            amrex::ParticleReal const zlen = 1.31879807_prt;
            amrex::ParticleReal const zmid = zlen/2.0_prt;
            int const ncoef = 25;

            // specify Fourier coefficients
            amrex::ParticleReal cos_coef [ncoef] = {};
            amrex::ParticleReal sin_coef [ncoef] = {};

            // recall array indexing begins with 0
            cos_coef[0] = 0.1644024074311037;
            cos_coef[1] = -0.1324009958969339;
            cos_coef[2] = 4.3443060026047219E-002;
            cos_coef[3] = 8.5602654094946495E-002;
            cos_coef[4] = -0.2433578169042885;
            cos_coef[5] = 0.5297150596779437;
            cos_coef[6] = 0.7164884680963959;
            cos_coef[7] = -5.2579522442877296E-003;
            cos_coef[8] = -5.5025369142193678E-002;
            cos_coef[9] = 4.6845673335028933E-002;
            cos_coef[10] = -2.3279346335638568E-002;
            cos_coef[11] = 4.0800777539657775E-003;
            cos_coef[12] = 4.1378326533752169E-003;
            cos_coef[13] = -2.5040533340490805E-003;
            cos_coef[14] = -4.0654981400000964E-003;
            cos_coef[15] = 9.6630592067498289E-003;
            cos_coef[16] = -8.5275895985990214E-003;
            cos_coef[17] = -5.8078747006425020E-002;
            cos_coef[18] = -2.4044337836660403E-002;
            cos_coef[19] = 1.0968240064697212E-002;
            cos_coef[20] = -3.4461179858301418E-003;
            cos_coef[21] = -8.1201564869443749E-004;
            cos_coef[22] = 2.1438992904959380E-003;
            cos_coef[23] = -1.4997753525697276E-003;
            cos_coef[24] = 1.8685171825676386E-004;

            // compute on-axis electric field (z is relative to cavity midpoint)
            amrex::ParticleReal efield = 0.0;
            amrex::ParticleReal efieldp = 0.0;
            amrex::ParticleReal efieldpp = 0.0;
            amrex::ParticleReal const z = zeval - zmid;

            if (abs(z)<=zmid)
            {
               efield = 0.5_prt*cos_coef[0];
               for(int j=1; j < ncoef; ++j)
               {
                 efield = efield + cos_coef[j]*cos(j*2*pi*z/zlen) + 
                      sin_coef[j]*sin(j*2*pi*z/zlen);
                 efieldp = efieldp-j*2*pi*cos_coef[j]*sin(j*2*pi*z/zlen)/zlen +
                      j*2*pi*sin_coef[j]*cos(j*2*pi*z/zlen)/zlen;
                 efieldpp = efieldpp- pow(j*2*pi*cos_coef[j]/zlen,2) *cos(j*2*pi*z/zlen) -
                      pow(j*2*pi*sin_coef[j]/zlen,2) *sin(j*2*pi*z/zlen);
               }
            }
            return efield;
        }

        void map1 (amrex::ParticleReal const tau, RefPart & refpart,
                               amrex::ParticleReal & zeval) const
        {

            using namespace amrex::literals; // for _rt and _prt


               amrex::ParticleReal const t = refpart.t;
               amrex::ParticleReal const pt = refpart.pt;

               if(pt < -1.0_prt){
                  refpart.t = t + tau/sqrt(1.0_prt - pow(pt,-2));
                  refpart.pt = pt;
               }
               else {
                 refpart.t = t;
                 refpart.pt = pt;
               }
               zeval = zeval + tau;

        }

        void map2 (amrex::ParticleReal const tau, RefPart & refpart,
                               amrex::ParticleReal & zeval) const
        {

            using namespace amrex::literals; // for _rt and _prt


               amrex::ParticleReal const t = refpart.t;
               amrex::ParticleReal const pt = refpart.pt;

               // Define parameters and intermediate constants
               amrex::ParticleReal const pi = 3.141592653589793_prt;
               amrex::ParticleReal const c = 2.99792458e8;
               amrex::ParticleReal const k = 2.0_prt*pi*m_freq/c;
               amrex::ParticleReal const phi = m_phase*pi/180.0_prt;
               amrex::ParticleReal const E0 = m_Escale;

               refpart.t = t;
               refpart.pt = pt - tau*E0*RF_Efield(zeval)*cos(k*t+phi);
               zeval = zeval;
               amrex::AllPrintToFile("rfout") << zeval << " " << RF_Efield(zeval) << "\n";
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
        int m_mapsteps; //! number of map integration steps per slice
        int m_nslice; //! number of slices used for the application of space charge
    };

} // namespace impactx

#endif // IMPACTX_RFCAVITY_H
