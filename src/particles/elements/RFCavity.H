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

#include <AMReX_Array.H>
#include <AMReX_Extension.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <tuple>


namespace impactx
{
    struct RFCavity
    {
        static constexpr auto name = "RFCavity";
        using PType = ImpactXParticleContainer::ParticleType;

        /** An RF cavity
         *
         * @param ds Segment length in m
         * @param escale Scaling factor for on-axis RF electric field
         * @param freq RF frequency in Hz
         * @param phase RF driven phase in deg
         * @param mapsteps number of integration steps per slice used for
         *        map and reference particle push in applied fields
         * @param nslice number of slices used for the application of space charge
         */
        RFCavity(
            amrex::ParticleReal const ds,
            amrex::ParticleReal const escale,
            amrex::ParticleReal const freq,
            amrex::ParticleReal const phase,
            int const mapsteps,
            int const nslice
        )
          : m_ds(ds), m_escale(escale), m_freq(freq), m_phase(phase),
            m_mapsteps(mapsteps), m_nslice(nslice)
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
            [[maybe_unused]] RefPart const & refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // access AoS data such as positions and cpu/id
            amrex::ParticleReal const x = p.pos(0);
            amrex::ParticleReal const y = p.pos(1);
            amrex::ParticleReal const t = p.pos(2);

            // initialize output values of momenta
            amrex::ParticleReal pxout = px;
            amrex::ParticleReal pyout = py;
            amrex::ParticleReal ptout = pt;

            // get the linear map
            amrex::Array2D<amrex::ParticleReal, 1, 6, 1, 6> const R = refpart.map;
            
            // symplectic linear map for the RF cavity is computed using the
            // Hamiltonian formalism as described in:
            // https://uspas.fnal.gov/materials/09UNM/ComputationalMethods.pdf.
            // R denotes the transfer matrix in the basis (x,px,y,py,t,pt),
            // so that, e.g., R(3,4) = dyf/dpyi.

            // push particles using the linear map
            p.pos(0) = R(1,1)*x + R(1,2)*px + R(1,3)*y
                     + R(1,4)*py + R(1,5)*t + R(1,6)*pt;
            pxout = R(2,1)*x + R(2,2)*px + R(2,3)*y
                  + R(2,4)*py + R(2,5)*t + R(2,6)*pt;
            p.pos(1) = R(3,1)*x + R(3,2)*px + R(3,3)*y
                     + R(3,4)*py + R(3,5)*t + R(3,6)*pt;
            pyout = R(4,1)*x + R(4,2)*px + R(4,3)*y
                  + R(4,4)*py + R(4,5)*t + R(4,6)*pt;
            p.pos(2) = R(5,1)*x + R(5,2)*px + R(5,3)*y
                     + R(5,4)*py + R(5,5)*t + R(5,6)*pt;
            ptout = R(6,1)*x + R(6,2)*px + R(6,3)*y
                  + R(6,4)*py + R(6,5)*t + R(6,6)*pt;

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
        void operator() (RefPart & AMREX_RESTRICT refpart) const {

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
            amrex::ParticleReal const sedge = refpart.sedge;

            // initialize linear map (deviation) values
            for (int i=1; i<7; i++) {
               for (int j=1; j<7; j++) {
                  if (i == j)
                      refpart.map(i,j) = 1.0_prt;
                  else
                      refpart.map(i,j) = 0.0_prt;
               }
            }

            // length of the current slice
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // compute intial value of beta*gamma
            amrex::ParticleReal const bgi = sqrt(pow(pt,2)-1.0_prt);

            // call integrator to advance (t,pt)
            amrex::ParticleReal const zin = s - sedge;
            amrex::ParticleReal const zout = zin + slice_ds;
            int const nsteps = m_mapsteps;

            integrators::symp2_integrate_split3(refpart,zin,zout,nsteps,*this);
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

            // convert linear map from dynamic to static units
            amrex::ParticleReal scale_in = 1.0_prt;
            amrex::ParticleReal scale_fin = 1.0_prt;

            for (int i=1; i<7; i++) {
               for (int j=1; j<7; j++) {
                   if( i % 2 == 0)
                      scale_fin = bgf;
                   else
                      scale_fin = 1.0_prt;
                   if( j % 2 == 0)
                      scale_in = bgi;
                   else
                      scale_in = 1.0_prt;
                   refpart.map(i,j) = refpart.map(i,j)*scale_in/scale_fin;
               }
            }

            // advance integrated path length
            refpart.s = s + slice_ds;
        }

        std::tuple<amrex::ParticleReal,amrex::ParticleReal,amrex::ParticleReal>
        RF_Efield (amrex::ParticleReal const zeval) const
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
            
            // Fourier coefficients for the on-axis longitudinal
            // electric field Ez of the 9-cell TESLA superconducting
            // cavity (default cavity model):
            https://journals.aps.org/prab/pdf/10.1103/PhysRevSTAB.3.092001.
            
            // May be superseded by user-provided input values (to do).

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
            amrex::ParticleReal efieldint = 0.0;
            amrex::ParticleReal const z = zeval - zmid;

            if (abs(z)<=zmid)
            {
               efield = 0.5_prt*cos_coef[0];
               efieldint = z*efield;
               for(int j=1; j < ncoef; ++j)
               {
                 efield = efield + cos_coef[j]*cos(j*2*pi*z/zlen) +
                      sin_coef[j]*sin(j*2*pi*z/zlen);
                 efieldp = efieldp-j*2*pi*cos_coef[j]*sin(j*2*pi*z/zlen)/zlen +
                      j*2*pi*sin_coef[j]*cos(j*2*pi*z/zlen)/zlen;
                 efieldpp = efieldpp- pow(j*2*pi*cos_coef[j]/zlen,2) *cos(j*2*pi*z/zlen) -
                      pow(j*2*pi*sin_coef[j]/zlen,2) *sin(j*2*pi*z/zlen);
                 efieldint = efieldint + zlen*cos_coef[j]*sin(j*2*pi*z/zlen)/(j*2*pi) -
                      zlen*sin_coef[j]*cos(j*2*pi*z/zlen)/(j*2*pi);
               }
            }
            return std::make_tuple(efield,efieldp,efieldint);
        }

        void map3 (amrex::ParticleReal const tau,
                   RefPart & refpart,
                   [[maybe_unused]] amrex::ParticleReal & zeval) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // push the reference particle
            amrex::ParticleReal const t = refpart.t;
            amrex::ParticleReal const pt = refpart.pt;

            if (pt < -1.0_prt) {
                refpart.t = t + tau/sqrt(1.0_prt - pow(pt,-2));
                refpart.pt = pt;
            }
            else {
                refpart.t = t;
                refpart.pt = pt;
            }

            // push the linear map equations
            amrex::Array2D<amrex::ParticleReal, 1, 6, 1, 6> const R = refpart.map;
            amrex::ParticleReal const betgam = refpart.beta_gamma();

            refpart.map(5,5) = R(5,5) + tau*R(6,5)/pow(betgam,3);
            refpart.map(5,6) = R(5,6) + tau*R(6,6)/pow(betgam,3);
        }

        void map2 (amrex::ParticleReal const tau,
                   RefPart & refpart,
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
            amrex::ParticleReal const E0 = m_escale;

            // push the reference particle
            [[maybe_unused]] auto [ez, ezp, ezint] = RF_Efield(zeval);

            refpart.t = t;
            refpart.pt = pt;

            // push the linear map equations
            amrex::Array2D<amrex::ParticleReal,1,6,1,6> const R = refpart.map;
            amrex::ParticleReal const s = tau/refpart.beta_gamma();
            amrex::ParticleReal const L = E0*ezp*sin(k*t+phi)/(2.0_prt*k);

            refpart.map(1,1) = (1.0_prt-s*L)*R(1,1) + s*R(2,1);
            refpart.map(1,2) = (1.0_prt-s*L)*R(1,2) + s*R(2,2);
            refpart.map(2,1) = -s*pow(L,2)*R(1,1) + (1.0_prt+s*L)*R(2,1);
            refpart.map(2,2) = -s*pow(L,2)*R(1,2) + (1.0_prt+s*L)*R(2,2);

            refpart.map(3,3) = (1.0_prt-s*L)*R(3,3) + s*R(4,3);
            refpart.map(3,4) = (1.0_prt-s*L)*R(3,4) + s*R(4,4);
            refpart.map(4,3) = -s*pow(L,2)*R(3,3) + (1.0_prt+s*L)*R(4,3);
            refpart.map(4,4) = -s*pow(L,2)*R(3,4) + (1.0_prt+s*L)*R(4,4);
        }

        void map1 (amrex::ParticleReal const tau,
                   RefPart & refpart,
                   amrex::ParticleReal & zeval) const
        {
            using namespace amrex::literals; // for _rt and _prt

            amrex::ParticleReal const t = refpart.t;
            amrex::ParticleReal const pt = refpart.pt;
            amrex::ParticleReal const z = zeval;

            // Define parameters and intermediate constants
            amrex::ParticleReal const pi = 3.141592653589793_prt;
            amrex::ParticleReal const c = 2.99792458e8;
            amrex::ParticleReal const k = 2.0_prt*pi*m_freq/c;
            amrex::ParticleReal const phi = m_phase*pi/180.0_prt;
            amrex::ParticleReal const E0 = m_escale;

            // push the reference particle
            [[maybe_unused]] auto [ez, ezp, ezint] = RF_Efield(z);
            zeval = z + tau;
            [[maybe_unused]] auto [ezf, ezpf, ezintf] = RF_Efield(zeval);

            refpart.t = t;
            refpart.pt = pt - E0*(ezintf-ezint)*cos(k*t+phi);

            // push the linear map equations
            amrex::Array2D<amrex::ParticleReal, 1, 6, 1, 6> const R = refpart.map;
            amrex::ParticleReal const M = E0*(ezintf-ezint)*k*sin(k*t+phi);
            amrex::ParticleReal const L = E0*(ezpf-ezp)*sin(k*t+phi)/(2.0_prt*k)+M/2.0_prt;

            refpart.map(2,1) = L*R(1,1) + R(2,1);
            refpart.map(2,2) = L*R(1,2) + R(2,2);

            refpart.map(4,3) = L*R(3,3) + R(4,3);
            refpart.map(4,4) = L*R(3,4) + R(4,4);

            refpart.map(6,5) = M*R(5,5) + R(6,5);
            refpart.map(6,6) = M*R(5,6) + R(6,6);
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

        /** Return the segment length
         *
         * @return value in meters
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::ParticleReal ds () const
        {
            return m_ds;
        }

    private:
        amrex::ParticleReal m_ds; //! segment length in m
        amrex::ParticleReal m_escale; //! scaling factor for RF electric field
        amrex::ParticleReal m_freq; //! RF frequency in Hz
        amrex::ParticleReal m_phase; //! RF driven phase in deg
        int m_mapsteps; //! number of map integration steps per slice
        int m_nslice; //! number of slices used for the application of space charge
    };

} // namespace impactx

#endif // IMPACTX_RFCAVITY_H
