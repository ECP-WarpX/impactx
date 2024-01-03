/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
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
     * electric field Ez of the 9-cell TESLA superconducting
     * cavity (default cavity model):
     * https://journals.aps.org/prab/pdf/10.1103/PhysRevSTAB.3.092001.
     *
     * Usually superseded by user-provided input values.
     */
    struct RF_field_data
    {
        amrex::Vector<amrex::ParticleReal> default_cos_coef = {
            0.1644024074311037,
            -0.1324009958969339,
            4.3443060026047219e-002,
            8.5602654094946495e-002,
            -0.2433578169042885,
            0.5297150596779437,
            0.7164884680963959,
            -5.2579522442877296e-003,
            -5.5025369142193678e-002,
            4.6845673335028933e-002,
            -2.3279346335638568e-002,
            4.0800777539657775e-003,
            4.1378326533752169e-003,
            -2.5040533340490805e-003,
            -4.0654981400000964e-003,
            9.6630592067498289e-003,
            -8.5275895985990214e-003,
            -5.8078747006425020e-002,
            -2.4044337836660403e-002,
            1.0968240064697212e-002,
            -3.4461179858301418e-003,
            -8.1201564869443749e-004,
            2.1438992904959380e-003,
            -1.4997753525697276e-003,
            1.8685171825676386e-004
        };

        amrex::Vector<amrex::ParticleReal> default_sin_coef = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0
        };
    };

/** Dynamic data for the RFCavity elements
 *
 * Since we copy the element to the device, we cannot store this data on the element itself.
 * But we can store pointers to this data with the element and keep a lookup table here,
 * which we clean up in the end.
 */
namespace RFCavityData
{
    //! last used id for a created RF cavity
    static inline int next_id = 0;

    //! host: cosine coefficients in Fourier expansion of on-axis electric field Ez
    static inline std::map<int, std::vector<amrex::ParticleReal>> h_cos_coef = {};
    //! host: sine coefficients in Fourier expansion of on-axis electric field Ez
    static inline std::map<int, std::vector<amrex::ParticleReal>> h_sin_coef = {};

    //! device: cosine coefficients in Fourier expansion of on-axis electric field Ez
    static inline std::map<int, amrex::Gpu::DeviceVector<amrex::ParticleReal>> d_cos_coef = {};
    //! device: sine coefficients in Fourier expansion of on-axis electric field Ez
    static inline std::map<int, amrex::Gpu::DeviceVector<amrex::ParticleReal>> d_sin_coef = {};

} // namespace RFCavityData

    struct RFCavity
    : public elements::BeamOptic<RFCavity>,
      public elements::Thick,
      public elements::Alignment
    {
        static constexpr auto name = "RFCavity";
        using PType = ImpactXParticleContainer::ParticleType;

        /** An RF cavity
         *
         * @param ds Segment length in m
         * @param escale Scaling factor for on-axis RF electric field in 1/m
         * @param freq RF frequency in Hz
         * @param phase RF driven phase in deg
         * @param cos_coef TODO
         * @param sin_coef TODO
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param mapsteps number of integration steps per slice used for
         *        map and reference particle push in applied fields
         * @param nslice number of slices used for the application of space charge
         */
        RFCavity (
            amrex::ParticleReal ds,
            amrex::ParticleReal escale,
            amrex::ParticleReal freq,
            amrex::ParticleReal phase,
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
            m_escale(escale), m_freq(freq), m_phase(phase), m_mapsteps(mapsteps)
        {
            // next created RF cavity has another id for its data
            RFCavityData::next_id++;

            // validate sin and cos coefficients are the same length
            m_ncoef = int(cos_coef.size());
            if (m_ncoef != int(sin_coef.size()))
                throw std::runtime_error("RFCavity: cos and sin coefficients must have same length!");

            // host data
            RFCavityData::h_cos_coef[m_id] = cos_coef;
            RFCavityData::h_sin_coef[m_id] = sin_coef;
            m_cos_h_data = RFCavityData::h_cos_coef[m_id].data();
            m_sin_h_data = RFCavityData::h_sin_coef[m_id].data();

            // device data
            RFCavityData::d_cos_coef.emplace(m_id, amrex::Gpu::DeviceVector<amrex::ParticleReal>(m_ncoef));
            RFCavityData::d_sin_coef.emplace(m_id, amrex::Gpu::DeviceVector<amrex::ParticleReal>(m_ncoef));
            amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                                  cos_coef.begin(), cos_coef.end(),
                                  RFCavityData::d_cos_coef[m_id].begin());
            amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                                  sin_coef.begin(), sin_coef.end(),
                                  RFCavityData::d_sin_coef[m_id].begin());
            amrex::Gpu::streamSynchronize();

            // low-level objects we can use on device
            m_cos_d_data = RFCavityData::d_cos_coef[m_id].data();
            m_sin_d_data = RFCavityData::d_sin_coef[m_id].data();
        }

        /** Push all particles */
        using BeamOptic::operator();

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

            // symplectic linear map for the RF cavity is computed using the
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
                  if (i == j)
                      refpart.map(i, j) = 1.0_prt;
                  else
                      refpart.map(i, j) = 0.0_prt;
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
                   refpart.map(i, j) = refpart.map(i, j) * scale_in / scale_fin;
               }
            }

            // advance integrated path length
            refpart.s = s + slice_ds;
        }

        /** This evaluates the on-axis RF electric field at a fixed location
         *  z, together with certain required integrals and derivatives.
         *  The field returned is normalized to a peak value of 1.
         *
         * @param zeval Longitudinal on-axis location in m
         */
        std::tuple<amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        RF_Efield (amrex::ParticleReal const zeval) const
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

            // compute on-axis electric field (z is relative to cavity midpoint)
            amrex::ParticleReal efield = 0.0;
            amrex::ParticleReal efieldp = 0.0;
            amrex::ParticleReal efieldpp = 0.0;
            amrex::ParticleReal efieldint = 0.0;
            amrex::ParticleReal const z = zeval - zmid;

            if (std::abs(z) <= zmid)
            {
               efield = 0.5_prt*cos_data[0];
               efieldint = z*efield;
               for (int j=1; j < m_ncoef; ++j)
               {
                 efield = efield + cos_data[j]*cos(j*2*pi*z/zlen) +
                     sin_data[j]*sin(j*2*pi*z/zlen);
                 efieldp = efieldp-j*2*pi*cos_data[j]*sin(j*2*pi*z/zlen)/zlen +
                      j*2*pi*sin_data[j]*cos(j*2*pi*z/zlen)/zlen;
                 efieldpp = efieldpp- pow(j*2*pi*cos_data[j]/zlen,2) *cos(j*2*pi*z/zlen) -
                      pow(j*2*pi*sin_data[j]/zlen,2) *sin(j*2*pi*z/zlen);
                 efieldint = efieldint + zlen*cos_data[j]*sin(j*2*pi*z/zlen)/(j*2*pi) -
                      zlen*sin_data[j]*cos(j*2*pi*z/zlen)/(j*2*pi);
               }
            }
            return std::make_tuple(efield, efieldp, efieldint);
        }

        /** This pushes the reference particle and the linear map matrix
         *  elements for an RF cavity through the symplectic map associated
         *  with H_3 in the Hamiltonian splitting H = H_1 + H_2 + H_3.
         *
         * @param tau Map step size in m
         * @param[in,out] refpart reference particle
         * @param[in,out] zeval Longitudinal on-axis location in m
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void map3 (amrex::ParticleReal const tau,
                   RefPart & refpart,
                   [[maybe_unused]] amrex::ParticleReal & zeval) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // push the reference particle
            amrex::ParticleReal const t = refpart.t;
            amrex::ParticleReal const pt = refpart.pt;

            if (pt < -1.0_prt) {
                refpart.t = t + tau/sqrt(1.0_prt - pow(pt, -2));
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

        /** This pushes the reference particle and the linear map matrix
         *  elements for an RF cavity through the symplectic map associated
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
            using ablastr::constant::math::pi;
            using ablastr::constant::SI::c;
            amrex::ParticleReal const k = (2.0_prt*pi/c)*m_freq;
            amrex::ParticleReal const phi = m_phase*(pi/180.0_prt);
            amrex::ParticleReal const E0 = m_escale;

            // push the reference particle
            auto [ez, ezp, ezint] = RF_Efield(zeval);
            amrex::ignore_unused(ez, ezint);

            refpart.t = t;
            refpart.pt = pt;

            // push the linear map equations
            amrex::Array2D<amrex::ParticleReal, 1, 6, 1, 6> const R = refpart.map;
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

        /** This pushes the reference particle and the linear map matrix
         *  elements for an RF cavity through the symplectic map associated
         *  with H_1 in the Hamiltonian splitting H = H_1 + H_2 + H_3.
         *
         * @param tau Map step size in m
         * @param[in,out] refpart reference particle
         * @param[in,out] zeval Longitudinal on-axis location in m
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void map1 (amrex::ParticleReal const tau,
                   RefPart & refpart,
                   amrex::ParticleReal & zeval) const
        {
            using namespace amrex::literals; // for _rt and _prt

            amrex::ParticleReal const t = refpart.t;
            amrex::ParticleReal const pt = refpart.pt;
            amrex::ParticleReal const z = zeval;

            // Define parameters and intermediate constants
            using ablastr::constant::math::pi;
            using ablastr::constant::SI::c;
            amrex::ParticleReal const k = (2.0_prt*pi/c)*m_freq;
            amrex::ParticleReal const phi = m_phase*(pi/180.0_prt);
            amrex::ParticleReal const E0 = m_escale;

            // push the reference particle
            auto [ez, ezp, ezint] = RF_Efield(z);
            amrex::ignore_unused(ez);
            zeval = z + tau;
            auto [ezf, ezpf, ezintf] = RF_Efield(zeval);
            amrex::ignore_unused(ezf);

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

        /** Close and deallocate all data and handles.
         */
        void
        finalize ()
        {
            // remove from unique data map
            if (RFCavityData::h_cos_coef.count(m_id) != 0u)
                RFCavityData::h_cos_coef.erase(m_id);
            if (RFCavityData::h_sin_coef.count(m_id) != 0u)
                RFCavityData::h_sin_coef.erase(m_id);

            if (RFCavityData::d_cos_coef.count(m_id) != 0u)
                RFCavityData::d_cos_coef.erase(m_id);
            if (RFCavityData::d_sin_coef.count(m_id) != 0u)
                RFCavityData::d_sin_coef.erase(m_id);
        }

    private:
        amrex::ParticleReal m_escale; //! scaling factor for RF electric field
        amrex::ParticleReal m_freq; //! RF frequency in Hz
        amrex::ParticleReal m_phase; //! RF driven phase in deg
        int m_mapsteps; //! number of map integration steps per slice
        int m_id; //! unique RF cavity id used for data lookup map

        int m_ncoef = 0; //! number of Fourier coefficients
        amrex::ParticleReal* m_cos_h_data = nullptr; //! non-owning pointer to host cosine coefficients
        amrex::ParticleReal* m_sin_h_data = nullptr; //! non-owning pointer to host sine coefficients
        amrex::ParticleReal* m_cos_d_data = nullptr; //! non-owning pointer to device cosine coefficients
        amrex::ParticleReal* m_sin_d_data = nullptr; //! non-owning pointer to device sine coefficients
    };

} // namespace impactx

#endif // IMPACTX_RFCAVITY_H
