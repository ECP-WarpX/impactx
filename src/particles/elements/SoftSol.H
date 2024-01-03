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
#include "mixin/alignment.H"
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
    * magnetic field Bz of a soft-edge solenoid based on
    * a thin-shell model of the solenoid as a cylinder
    * of current, as described in:
    * https://doi.org/10.1016/j.nima.2022.166706.
    * Parameters correspond to an aspect ratio of
    * cylinder diameter / length = 1/3.
    *
    * Usually superseded by user-provided input values.
    */
    struct Sol_field_data
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

/** Dynamic data for the SoftSolenoid elements
 *
 * Since we copy the element to the device, we cannot store this data on the element itself.
 * But we can store pointers to this data with the element and keep a lookup table here,
 * which we clean up in the end.
 */
namespace SoftSolenoidData
{
    //! last used id for a created soft solenoid
    static inline int next_id = 0;

    //! host: cosine coefficients in Fourier expansion of on-axis magnetic field Bz
    static inline std::map<int, std::vector<amrex::ParticleReal>> h_cos_coef = {};
    //! host: sine coefficients in Fourier expansion of on-axis magnetic field Bz
    static inline std::map<int, std::vector<amrex::ParticleReal>> h_sin_coef = {};

    //! device: cosine coefficients in Fourier expansion of on-axis magnetic field Bz
    static inline std::map<int, amrex::Gpu::DeviceVector<amrex::ParticleReal>> d_cos_coef = {};
    //! device: sine coefficients in Fourier expansion of on-axis magnetic field Bz
    static inline std::map<int, amrex::Gpu::DeviceVector<amrex::ParticleReal>> d_sin_coef = {};

} // namespace SoftSolenoidData

    struct SoftSolenoid
    : public elements::BeamOptic<SoftSolenoid>,
      public elements::Thick,
      public elements::Alignment
    {
        static constexpr auto name = "SoftSolenoid";
        using PType = ImpactXParticleContainer::ParticleType;

        /** A soft-edge solenoid
         *
         * @param ds Segment length in m
         * @param bscale Scaling factor for on-axis magnetic field Bz in 1/m
         * @param cos_coef cosine coefficients in Fourier expansion of on-axis magnetic field Bz
         * @param sin_coef sine coefficients in Fourier expansion of on-axis magnetic field Bz
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param mapsteps number of integration steps per slice used for
         *        map and reference particle push in applied fields
         * @param nslice number of slices used for the application of space charge
         */
        SoftSolenoid (
            amrex::ParticleReal ds,
            amrex::ParticleReal bscale,
            std::vector<amrex::ParticleReal> cos_coef,
            std::vector<amrex::ParticleReal> sin_coef,
            amrex::ParticleReal dx = 0,
            amrex::ParticleReal dy = 0,
            amrex::ParticleReal rotation_degree = 0,
            int mapsteps = 1,
            int nslice = 1
        )
          : Thick(ds, nslice),
            Alignment(dx, dy, rotation_degree),
            m_bscale(bscale), m_mapsteps(mapsteps), m_id(SoftSolenoidData::next_id)
       {
           // next created soft solenoid has another id for its data
           SoftSolenoidData::next_id++;

           // validate sin and cos coefficients are the same length
           m_ncoef = int(cos_coef.size());
           if (m_ncoef != int(sin_coef.size()))
               throw std::runtime_error("SoftSolenoid: cos and sin coefficients must have same length!");

           // host data
           SoftSolenoidData::h_cos_coef[m_id] = cos_coef;
           SoftSolenoidData::h_sin_coef[m_id] = sin_coef;
           m_cos_h_data = SoftSolenoidData::h_cos_coef[m_id].data();
           m_sin_h_data = SoftSolenoidData::h_sin_coef[m_id].data();

           // device data
           SoftSolenoidData::d_cos_coef.emplace(m_id, amrex::Gpu::DeviceVector<amrex::ParticleReal>(m_ncoef));
           SoftSolenoidData::d_sin_coef.emplace(m_id, amrex::Gpu::DeviceVector<amrex::ParticleReal>(m_ncoef));
           amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                                 cos_coef.begin(), cos_coef.end(),
                                 SoftSolenoidData::d_cos_coef[m_id].begin());
           amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                                 sin_coef.begin(), sin_coef.end(),
                                 SoftSolenoidData::d_sin_coef[m_id].begin());
           amrex::Gpu::streamSynchronize();

           // low-level objects we can use on device
           m_cos_d_data = SoftSolenoidData::d_cos_coef[m_id].data();
           m_sin_d_data = SoftSolenoidData::d_sin_coef[m_id].data();
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

            // shift due to alignment errors of the element
            shift_in(p.pos(RealAoS::x), p.pos(RealAoS::y), px, py);

            // access AoS data such as positions and cpu/id
            amrex::ParticleReal const x = p.pos(RealAoS::x);
            amrex::ParticleReal const y = p.pos(RealAoS::y);
            amrex::ParticleReal const t = p.pos(RealAoS::t);

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
            p.pos(RealAoS::x) = R(1,1)*x + R(1,2)*px + R(1,3)*y
                     + R(1,4)*py + R(1,5)*t + R(1,6)*pt;
            pxout = R(2,1)*x + R(2,2)*px + R(2,3)*y
                  + R(2,4)*py + R(2,5)*t + R(2,6)*pt;
            p.pos(RealAoS::y) = R(3,1)*x + R(3,2)*px + R(3,3)*y
                     + R(3,4)*py + R(3,5)*t + R(3,6)*pt;
            pyout = R(4,1)*x + R(4,2)*px + R(4,3)*y
                  + R(4,4)*py + R(4,5)*t + R(4,6)*pt;
            p.pos(RealAoS::t) = R(5,1)*x + R(5,2)*px + R(5,3)*y
                     + R(5,4)*py + R(5,5)*t + R(5,6)*pt;
            ptout = R(6,1)*x + R(6,2)*px + R(6,3)*y
                  + R(6,4)*py + R(6,5)*t + R(6,6)*pt;

            // assign updated momenta
            px = pxout;
            py = pyout;
            pt = ptout;

            // undo shift due to alignment errors of the element
            shift_out(p.pos(RealAoS::x), p.pos(RealAoS::y), px, py);
        }

        /** This pushes the reference particle.
         *
         * @param[in,out] refpart reference particle
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
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

            // pick the right data depending if we are on the host side
            // (reference particle push) or device side (particles):
#if AMREX_DEVICE_COMPILE
            amrex::ParticleReal* cos_data = m_cos_d_data;
            amrex::ParticleReal* sin_data = m_sin_d_data;
#else
            amrex::ParticleReal* cos_data = m_cos_h_data;
            amrex::ParticleReal* sin_data = m_sin_h_data;
#endif

            // specify constants
            using ablastr::constant::math::pi;
            amrex::ParticleReal const zlen = m_ds;
            amrex::ParticleReal const zmid = zlen / 2.0_prt;

            // compute on-axis magnetic field (z is relative to solenoid midpoint)
            amrex::ParticleReal bfield = 0.0;
            amrex::ParticleReal bfieldp = 0.0;
            amrex::ParticleReal bfieldint = 0.0;
            amrex::ParticleReal const z = zeval - zmid;

            if (std::abs(z) <= zmid)
            {
               bfield = 0.5_prt*cos_data[0];
               bfieldint = z*bfield;
               for (int j=1; j < m_ncoef; ++j)
               {
                 bfield = bfield + cos_data[j]*cos(j*2*pi*z/zlen) +
                     sin_data[j]*sin(j*2*pi*z/zlen);
                 bfieldp = bfieldp-j*2*pi*cos_data[j]*sin(j*2*pi*z/zlen)/zlen +
                      j*2*pi*sin_data[j]*cos(j*2*pi*z/zlen)/zlen;
                 bfieldint = bfieldint + zlen*cos_data[j]*sin(j*2*pi*z/zlen)/(j*2*pi) -
                      zlen*sin_data[j]*cos(j*2*pi*z/zlen)/(j*2*pi);
               }
            }
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

        /** Close and deallocate all data and handles.
         */
        void
        finalize ()
        {
            // remove from unique data map
            if (SoftSolenoidData::h_cos_coef.count(m_id) != 0u)
                SoftSolenoidData::h_cos_coef.erase(m_id);
            if (SoftSolenoidData::h_sin_coef.count(m_id) != 0u)
                SoftSolenoidData::h_sin_coef.erase(m_id);

            if (SoftSolenoidData::d_cos_coef.count(m_id) != 0u)
                SoftSolenoidData::d_cos_coef.erase(m_id);
            if (SoftSolenoidData::d_sin_coef.count(m_id) != 0u)
                SoftSolenoidData::d_sin_coef.erase(m_id);
        }

    private:
        amrex::ParticleReal m_bscale; //! scaling factor for solenoid Bz field
        int m_mapsteps; //! number of map integration steps per slice
        int m_id; //! unique soft solenoid id used for data lookup map

        int m_ncoef = 0; //! number of Fourier coefficients
        amrex::ParticleReal* m_cos_h_data = nullptr; //! non-owning pointer to host cosine coefficients
        amrex::ParticleReal* m_sin_h_data = nullptr; //! non-owning pointer to host sine coefficients
        amrex::ParticleReal* m_cos_d_data = nullptr; //! non-owning pointer to device cosine coefficients
        amrex::ParticleReal* m_sin_d_data = nullptr; //! non-owning pointer to device sine coefficients
    };

} // namespace impactx

#endif // IMPACTX_SOFTSOL_H
