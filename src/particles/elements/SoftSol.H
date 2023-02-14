/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_SOFTSOL_H
#define IMPACTX_SOFTSOL_H

#include "particles/ImpactXParticleContainer.H"
#include "particles/integrators/Integrators.H"
#include "mixin/beamoptic.H"
#include "mixin/thick.H"

#include <ablastr/constant.H>

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_Extension.H>
#include <AMReX_REAL.H>

#include <array>
#include <cmath>
#include <tuple>


namespace impactx
{
    struct Sol_field_data
    {
       //   Default:
       //   Fourier coefficients for the on-axis longitudinal
       //   magnetic field Bz of a soft-edge solenoid based on
       //   the design used in the LCLS-II photoinjector.
       //   Superseded by user-provided input values.

       std::vector<double> default_cos_coef = {
            0.537354892644806,
            0.405516997152622,
            0.204379467431097,
            8.818833164290837E-002,
            2.990035664995536E-002,
            8.103036020221078E-003,
           -3.457107719342389E-004,
           -1.506042059246567E-003,
           -1.684741337807429E-003,
           -7.496351690988198E-004,
           -5.231073743536983E-004,
           -3.074520833708134E-005,
           -8.136367978784272E-005,
            1.041748970155938E-004,
           -2.263763170923048E-005,
            7.351686613338790E-005,
           -3.220999221920476E-005,
            4.158432715586470E-005,
           -3.417704589118791E-005,
            2.861088103131010E-005,
           -2.787178123388915E-005,
            2.282090995044119E-005,
           -2.181003508883439E-005,
            1.963104345689975E-005,
           -1.626901477387445E-005,
            };
       std::vector<double> default_sin_coef = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0
            };
    };

    struct SoftSolenoid
    : public elements::BeamOptic<SoftSolenoid>,
      public elements::Thick
    {
        static constexpr auto name = "SoftSolenoid";
        using PType = ImpactXParticleContainer::ParticleType;

        /** A soft-edge solenoid
         *
         * @param ds Segment length in m
         * @param bscale Scaling factor for on-axis magnetic field Bz in 1/m
         * @param mapsteps number of integration steps per slice used for
         *        map and reference particle push in applied fields
         * @param nslice number of slices used for the application of space charge
         */
        SoftSolenoid(
            amrex::ParticleReal const ds,
            amrex::ParticleReal const bscale,
            int const mapsteps,
            int const nslice,
            std::vector<double> cos_coef,
            std::vector<double> sin_coef
        )
          : Thick(ds, nslice),
            m_bscale(bscale), m_mapsteps(mapsteps), m_cos_coef(cos_coef),
            m_sin_coef(sin_coef)
        {
        }

        /** Push all particles */
        using BeamOptic::operator();

        /** This is a soft-edge solenoid functor, so that a variable of this type can be used
         *  like a soft-edge solenoid function.
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

            // symplectic linear map for a solenoid is computed using the
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
        void operator() (RefPart & AMREX_RESTRICT refpart) const
        {
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
                  auto const default_value = (i == j) ? 1.0_prt : 0.0_prt;
                  refpart.map(i, j) = default_value;
               }
            }

            // length of the current slice
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // compute intial value of beta*gamma
            amrex::ParticleReal const bgi = sqrt(pow(pt, 2) - 1.0_prt);

            // call integrator to advance (t,pt)
            amrex::ParticleReal const zin = s - sedge;
            amrex::ParticleReal const zout = zin + slice_ds;
            int const nsteps = m_mapsteps;

            integrators::symp2_integrate_split3(refpart,zin,zout,nsteps,*this);
            amrex::ParticleReal const ptf = refpart.pt;

            /* print computed linear map:
               for(int i=1; i<7; ++i){
                 for(int j=1; j<7; ++j){
                    amrex::PrintToFile("SolMap.txt") << i << " " <<
                    j << " " << refpart.map(i,j) << "\n";
                 }
               }
            */

            // advance position (x,y,z)
            refpart.x = x + slice_ds*px/bgi;
            refpart.y = y + slice_ds*py/bgi;
            refpart.z = z + slice_ds*pz/bgi;

            // compute final value of beta*gamma
            amrex::ParticleReal const bgf = sqrt(pow(ptf, 2) - 1.0_prt);

            // advance momentum (px,py,pz)
            refpart.px = px*bgf/bgi;
            refpart.py = py*bgf/bgi;
            refpart.pz = pz*bgf/bgi;

            // advance integrated path length
            refpart.s = s + slice_ds;
        }

        /** This evaluates the on-axis magnetic field Bz at a fixed location
         *  z, together with certain required integrals and derivatives.
         *  The field returned is normalized to a peak value of 1.
         *
         * @param zeval Longitudinal on-axis location in m
         */
        std::tuple<amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        Sol_Bfield (amrex::ParticleReal const zeval) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // specify constants
            using ablastr::constant::math::pi;
            amrex::ParticleReal const zlen = m_ds;
            amrex::ParticleReal const zmid = zlen / 2.0_prt;

            // compute on-axis magnetic field (z is relative to cavity midpoint)
            amrex::ParticleReal bfield = 0.0;
            amrex::ParticleReal bfieldp = 0.0;
            amrex::ParticleReal bfieldpp = 0.0;
            amrex::ParticleReal bfieldint = 0.0;
            amrex::ParticleReal const z = zeval - zmid;
            int const ncoef = m_cos_coef.size();

            if (abs(z)<=zmid)
            {
               bfield = 0.5_prt*m_cos_coef[0];
               bfieldint = z*bfield;
               for (int j=1; j < ncoef; ++j)
               {
                 bfield = bfield + m_cos_coef[j]*cos(j*2*pi*z/zlen) +
                      m_sin_coef[j]*sin(j*2*pi*z/zlen);
                 bfieldp = bfieldp-j*2*pi*m_cos_coef[j]*sin(j*2*pi*z/zlen)/zlen +
                      j*2*pi*m_sin_coef[j]*cos(j*2*pi*z/zlen)/zlen;
                 bfieldpp = bfieldpp- pow(j*2*pi*m_cos_coef[j]/zlen,2) *cos(j*2*pi*z/zlen) -
                      pow(j*2*pi*m_sin_coef[j]/zlen,2) *sin(j*2*pi*z/zlen);
                 bfieldint = bfieldint + zlen*m_cos_coef[j]*sin(j*2*pi*z/zlen)/(j*2*pi) -
                      zlen*m_sin_coef[j]*cos(j*2*pi*z/zlen)/(j*2*pi);
               }
            }
            amrex::PrintToFile("Bfield.txt") << z << " " << bfield << "\n";
            return std::make_tuple(bfield, bfieldp, bfieldint);
        }

        /** This pushes the reference particle and the linear map matrix
         *  elements for a solenoid through the symplectic map associated
         *  with H_3 in the Hamiltonian splitting H = H_1 + H_2 + H_3.
         *
         * @param tau Map step size in m
         * @param[in,out] refpart reference particle
         * @param[in,out] zeval Longitudinal on-axis location in m
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void map1 (amrex::ParticleReal const tau,
                   RefPart & refpart,
                   [[maybe_unused]] amrex::ParticleReal & zeval) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // push the reference particle
            amrex::ParticleReal const t = refpart.t;
            amrex::ParticleReal const pt = refpart.pt;
            amrex::ParticleReal const z = zeval;

            if (pt < -1.0_prt) {
                refpart.t = t + tau/sqrt(1.0_prt - pow(pt, -2));
                refpart.pt = pt;
            }
            else {
                refpart.t = t;
                refpart.pt = pt;
            }

            zeval = z + tau;

            // push the linear map equations
            amrex::Array2D<amrex::ParticleReal, 1, 6, 1, 6> const R = refpart.map;
            amrex::ParticleReal const betgam = refpart.beta_gamma();

            refpart.map(1,1) = R(1,1) + tau*R(2,1);
            refpart.map(1,2) = R(1,2) + tau*R(2,2);
            refpart.map(1,3) = R(1,3) + tau*R(2,3);
            refpart.map(1,4) = R(1,4) + tau*R(2,4);

            refpart.map(3,1) = R(3,1) + tau*R(4,1);
            refpart.map(3,2) = R(3,2) + tau*R(4,2);
            refpart.map(3,3) = R(3,3) + tau*R(4,3);
            refpart.map(3,4) = R(3,4) + tau*R(4,4);

            refpart.map(5,5) = R(5,5) + tau*R(6,5)/pow(betgam,2);
            refpart.map(5,6) = R(5,6) + tau*R(6,6)/pow(betgam,2);

        }

        /** This pushes the reference particle and the linear map matrix
         *  elements for a solenoid through the symplectic map associated
         *  with H_2 in the Hamiltonian splitting H = H_1 + H_2 + H_3.
         *
         * @param tau Map step size in m
         * @param[in,out] refpart reference particle
         * @param[in,out] zeval Longitudinal on-axis location in m
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void map2 (amrex::ParticleReal const tau,
                   RefPart & refpart,
                   amrex::ParticleReal & zeval) const
        {
            using namespace amrex::literals; // for _rt and _prt

            amrex::ParticleReal const t = refpart.t;
            amrex::ParticleReal const pt = refpart.pt;

            // Define parameters and intermediate constants
            amrex::ParticleReal const B0 = m_bscale;

            // push the reference particle
            auto [bz, bzp, bzint] = Sol_Bfield(zeval);
            amrex::ignore_unused(bzp, bzint);

            refpart.t = t;
            refpart.pt = pt;

            // push the linear map equations
            amrex::Array2D<amrex::ParticleReal, 1, 6, 1, 6> const R = refpart.map;
            amrex::ParticleReal const alpha = B0*bz/2.0_prt;
            amrex::ParticleReal const alpha2 = pow(alpha,2);

            refpart.map(2,1) = R(2,1) - tau*alpha2*R(1,1);
            refpart.map(2,2) = R(2,2) - tau*alpha2*R(1,2);
            refpart.map(2,3) = R(2,3) - tau*alpha2*R(1,3);
            refpart.map(2,4) = R(2,4) - tau*alpha2*R(1,4);

            refpart.map(4,1) = R(4,1) - tau*alpha2*R(3,1);
            refpart.map(4,2) = R(4,2) - tau*alpha2*R(3,2);
            refpart.map(4,3) = R(4,3) - tau*alpha2*R(3,3);
            refpart.map(4,4) = R(4,4) - tau*alpha2*R(3,4);

        }

        /** This pushes the reference particle and the linear map matrix
         *  elements for a solenoid through the symplectic map associated
         *  with H_1 in the Hamiltonian splitting H = H_1 + H_2 + H_3.
         *
         * @param tau Map step size in m
         * @param[in,out] refpart reference particle
         * @param[in,out] zeval Longitudinal on-axis location in m
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void map3 (amrex::ParticleReal const tau,
                   RefPart & refpart,
                   amrex::ParticleReal & zeval) const
        {
            using namespace amrex::literals; // for _rt and _prt

            amrex::ParticleReal const t = refpart.t;
            amrex::ParticleReal const pt = refpart.pt;
            amrex::ParticleReal const z = zeval;

            // Define parameters and intermediate constants
            amrex::ParticleReal const B0 = m_bscale;

            // push the reference particle
            auto [bz, bzp, bzint] = Sol_Bfield(z);
            amrex::ignore_unused(bzp, bzint);

            refpart.t = t;
            refpart.pt = pt;

            // push the linear map equations
            amrex::Array2D<amrex::ParticleReal, 1, 6, 1, 6> const R = refpart.map;
            amrex::ParticleReal const theta = tau*B0*bz/2.0_prt;
            amrex::ParticleReal const cs = cos(theta);
            amrex::ParticleReal const sn = sin(theta);

            refpart.map(1,1) = R(1,1)*cs + R(3,1)*sn;
            refpart.map(1,2) = R(1,2)*cs + R(3,2)*sn;
            refpart.map(1,3) = R(1,3)*cs + R(3,3)*sn;
            refpart.map(1,4) = R(1,4)*cs + R(3,4)*sn;

            refpart.map(2,1) = R(2,1)*cs + R(4,1)*sn;
            refpart.map(2,2) = R(2,2)*cs + R(4,2)*sn;
            refpart.map(2,3) = R(2,3)*cs + R(4,3)*sn;
            refpart.map(2,4) = R(2,4)*cs + R(4,4)*sn;

            refpart.map(3,1) = R(3,1)*cs - R(1,1)*sn;
            refpart.map(3,2) = R(3,2)*cs - R(1,2)*sn;
            refpart.map(3,3) = R(3,3)*cs - R(1,3)*sn;
            refpart.map(3,4) = R(3,4)*cs - R(1,4)*sn;

            refpart.map(4,1) = R(4,1)*cs - R(2,1)*sn;
            refpart.map(4,2) = R(4,2)*cs - R(2,2)*sn;
            refpart.map(4,3) = R(4,3)*cs - R(2,3)*sn;
            refpart.map(4,4) = R(4,4)*cs - R(2,4)*sn;

        }

    private:
        amrex::ParticleReal m_bscale; //! scaling factor for solenoid Bz field
        int m_mapsteps; //! number of map integration steps per slice
        std::vector<double> m_cos_coef;
        std::vector<double> m_sin_coef;
    };

} // namespace impactx

#endif // IMPACTX_SOFTSOL_H