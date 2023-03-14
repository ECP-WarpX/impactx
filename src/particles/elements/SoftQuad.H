/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_SOFTQUAD_H
#define IMPACTX_SOFTQUAD_H

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
#include <stdexcept>
#include <tuple>
#include <vector>


namespace impactx
{
   /** Default Fourier coefficients
    *
    * Fourier coefficients for the on-axis longitudinal
    * magnetic field Bz of a soft-edge quadrupole based on
    *
    * Usually superseded by user-provided input values.
    */
    struct Quad_field_data
    {
       amrex::Vector<amrex::ParticleReal> default_cos_coef = {
             0.350807812299706,
             0.323554693720069,
             0.260320578919415,
             0.182848575294969,
             0.106921016050403,
             4.409581845710694E-002,
            -9.416427163897508E-006,
            -2.459452716865687E-002,
            -3.272762575737291E-002,
            -2.936414401076162E-002,
            -1.995780078926890E-002,
            -9.102893342953847E-003,
            -2.456410658713271E-006,
             5.788233017324325E-003,
             8.040408292420691E-003,
             7.480064552867431E-003,
             5.230254569468851E-003,
             2.447614547094685E-003,
            -1.095525090532255E-006,
            -1.614586867387170E-003,
            -2.281365457438345E-003,
            -2.148709081338292E-003,
            -1.522541739363011E-003,
            -7.185505862719508E-004,
            -6.171194824600157E-007,
             4.842109305036943E-004,
             6.874508102002901E-004,
             6.535550288205728E-004,
             4.648795813759210E-004,
             2.216564722797528E-004,
            -4.100982995210341E-007,
            -1.499332112463395E-004,
            -2.151538438342482E-004,
            -2.044590946652016E-004,
            -1.468242784844341E-004
            };

       amrex::Vector<amrex::ParticleReal> default_sin_coef = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0
            };
    };

namespace data
{
    /** Data members we can copy to device with a memcpy.
     */
    struct SoftQuadrupole_device_copyable
    {
        amrex::ParticleReal m_gscale; //! scaling factor for quad field gradient
        int m_mapsteps; //! number of map integration steps per slice

        int m_ncoef = 0; //! number of Fourier coefficients
        amrex::ParticleReal* m_cos_data = nullptr; //! non-owning pointer to device cosine coefficients
        amrex::ParticleReal* m_sin_data = nullptr; //! non-owning pointer to device sine coefficients

        SoftQuadrupole_device_copyable(
            amrex::ParticleReal gscale,
            int mapsteps = 1
        ) : m_gscale(gscale), m_mapsteps(mapsteps)
        {}

        SoftQuadrupole_device_copyable(SoftQuadrupole_device_copyable const &) = default;
        SoftQuadrupole_device_copyable& operator=(SoftQuadrupole_device_copyable const &) = default;
        SoftQuadrupole_device_copyable(SoftQuadrupole_device_copyable &&) = default;
        SoftQuadrupole_device_copyable& operator=(SoftQuadrupole_device_copyable &&) = default;

    };
} // namespace data

    struct SoftQuadrupole
    : public elements::BeamOptic<SoftQuadrupole>,
      public elements::Thick,
      public data::SoftQuadrupole_device_copyable
    {
        static constexpr auto name = "SoftQuadrupole";
        using PType = ImpactXParticleContainer::ParticleType;

        /** A soft-edge quadrupole
         *
         * @param ds Segment length in m
         * @param gscale Scaling factor for on-axis field gradient Bz in 1/m
         * @param cos_coef TODO
         * @param sin_coef TODO
         * @param mapsteps number of integration steps per slice used for
         *        map and reference particle push in applied fields
         * @param nslice number of slices used for the application of space charge
         */
        AMREX_GPU_HOST
        SoftQuadrupole(
            amrex::ParticleReal ds,
            amrex::ParticleReal gscale,
            std::vector<amrex::ParticleReal> cos_coef,
            std::vector<amrex::ParticleReal> sin_coef,
            int mapsteps = 1,
            int nslice = 1
        )
          : Thick(ds, nslice),
            SoftQuadrupole_device_copyable(gscale, mapsteps)
       {
           m_ncoef = cos_coef.size();
           if (m_ncoef !=  int(sin_coef.size()))
                throw std::runtime_error("SoftQuadrupole: cos and sin coefficients must have same length!");

            m_cos_coef.resize(m_ncoef);
            m_sin_coef.resize(m_ncoef);
            amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                                  cos_coef.begin(), cos_coef.end(),
                                  m_cos_coef.begin());
            amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                                  sin_coef.begin(), sin_coef.end(),
                                  m_sin_coef.begin());
            amrex::Gpu::synchronize();

            // low-level objects we can use on device
            m_cos_data = m_cos_coef.data();
            m_sin_data = m_sin_coef.data();
        }

        // copy and move constructors
        SoftQuadrupole(SoftQuadrupole const & other)
          : Thick(other.m_ds, other.m_nslice),
            SoftQuadrupole_device_copyable(other.m_gscale, other.m_mapsteps)
        {
#if !AMREX_DEVICE_COMPILE
            // copy the data container if we copy the host element
            m_cos_coef = other.m_cos_coef;
            m_sin_coef = other.m_sin_coef;
            amrex::Gpu::synchronize();
#endif
            m_ncoef = m_cos_coef.size();
            m_cos_data = m_cos_coef.data();
            m_sin_data = m_sin_coef.data();
        }
        SoftQuadrupole& operator=(SoftQuadrupole const& other)
        {
            if (this == &other)
                return *this;

            Thick::operator=(other);
            SoftQuadrupole_device_copyable::operator=(other);
#if !AMREX_DEVICE_COMPILE
            // copy the data container if we copy the host element
            m_cos_coef = other.m_cos_coef;
            m_sin_coef = other.m_sin_coef;
            amrex::Gpu::synchronize();
#endif
            m_ncoef = m_cos_coef.size();
            m_cos_data = m_cos_coef.data();
            m_sin_data = m_sin_coef.data();
            return *this;
        }
        SoftQuadrupole(SoftQuadrupole && other) = default;
        SoftQuadrupole& operator=(SoftQuadrupole && other) = default;

        AMREX_GPU_HOST
        ~SoftQuadrupole() = default;

        /** Push all particles */
        using BeamOptic::operator();

        /** This is a soft-edge quadrupole functor, so that a variable of this type can be used
         *  like a soft-edge quadrupole function.
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

            // symplectic linear map for a quadrupole is computed using the
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
                    amrex::PrintToFile("QuadMap.txt") << i << " " <<
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
        Quad_Bfield (amrex::ParticleReal const zeval) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // specify constants
            using ablastr::constant::math::pi;
            amrex::ParticleReal const zlen = m_ds;
            amrex::ParticleReal const zmid = zlen / 2.0_prt;

            // compute on-axis magnetic field (z is relative to quadrupole midpoint)
            amrex::ParticleReal bfield = 0.0;
            amrex::ParticleReal bfieldp = 0.0;
            amrex::ParticleReal bfieldint = 0.0;
            amrex::ParticleReal const z = zeval - zmid;

            if (abs(z)<=zmid)
            {
               bfield = 0.5_prt*m_cos_data[0];
               bfieldint = z*bfield;
               for (int j=1; j < m_ncoef; ++j)
               {
                 bfield = bfield + m_cos_data[j]*cos(j*2*pi*z/zlen) +
                      m_sin_data[j]*sin(j*2*pi*z/zlen);
                 bfieldp = bfieldp-j*2*pi*m_cos_data[j]*sin(j*2*pi*z/zlen)/zlen +
                      j*2*pi*m_sin_data[j]*cos(j*2*pi*z/zlen)/zlen;
                 bfieldint = bfieldint + zlen*m_cos_data[j]*sin(j*2*pi*z/zlen)/(j*2*pi) -
                      zlen*m_sin_data[j]*cos(j*2*pi*z/zlen)/(j*2*pi);
               }
            }
            return std::make_tuple(bfield, bfieldp, bfieldint);
        }

        /** This pushes the reference particle and the linear map matrix
         *  elements for a quadrupole through the symplectic map associated
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
         *  elements for a quadrupole through the symplectic map associated
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
            amrex::ParticleReal const G0 = m_gscale;

            // push the reference particle
            auto [bz, bzp, bzint] = Quad_Bfield(zeval);
            amrex::ignore_unused(bzp, bzint);

            refpart.t = t;
            refpart.pt = pt;

            // push the linear map equations
            amrex::Array2D<amrex::ParticleReal, 1, 6, 1, 6> const R = refpart.map;
            amrex::ParticleReal const alpha = G0*bz/2.0_prt;
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
         *  elements for a quadrupole through the symplectic map associated
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
            amrex::ParticleReal const G0 = m_gscale;

            // push the reference particle
            auto [bz, bzp, bzint] = Quad_Bfield(z);
            amrex::ignore_unused(bzp, bzint);

            refpart.t = t;
            refpart.pt = pt;

            // push the linear map equations
            amrex::Array2D<amrex::ParticleReal, 1, 6, 1, 6> const R = refpart.map;
            amrex::ParticleReal const theta = tau*G0*bz/2.0_prt;
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
        // we cannot copy these to device with a memcpy when we copy the element class
        amrex::Gpu::DeviceVector<amrex::ParticleReal> m_cos_coef; //! cosine coefficients in Fourier expansion of on-axis magnetic field Bz
        amrex::Gpu::DeviceVector<amrex::ParticleReal> m_sin_coef; //! sine coefficients in Fourier expansion of on-axis magnetic  field Bz
    };

} // namespace impactx

#endif // IMPACTX_SOFTQUAD_H
