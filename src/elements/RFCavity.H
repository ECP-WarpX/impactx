/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
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
#include "mixin/dynamicdata.H"
#include "mixin/lineartransport.H"
#include "mixin/named.H"
#include "mixin/nofinalize.H"
#include "mixin/pipeaperture.H"
#include "mixin/thick.H"
#include "mixin/spintransport.H"

#include <ablastr/constant.H>

#include <AMReX.H>
#include <AMReX_Extension.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>
#include <AMReX_SIMD.H>
#include <AMReX_SmallMatrix.H>
#include <AMReX_TrackedVector.H>

#include <cmath>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace impactx::elements
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

    /** Fourier coefficients for on-axis RF electric field Ez.
     *
     * Host data persists across simulations; device data is uploaded
     * on-demand "lazily" on first use and released when AMReX shuts down.
     */
    struct CavityFourierCoefficients
    {
        amrex::Gpu::TrackedVector<amrex::ParticleReal> cos;
        amrex::Gpu::TrackedVector<amrex::ParticleReal> sin;
    };

    /** On-axis wake function as specified by the user.
     *
     * Host data persists across simulations; device data is uploaded
     * on-demand "lazily" on first use and released when AMReX shuts down.
     */
    struct CavityWakeFunction
    {
        amrex::Gpu::TrackedVector<amrex::ParticleReal> z_data;
        amrex::Gpu::TrackedVector<amrex::ParticleReal> wake_x_data;
        amrex::Gpu::TrackedVector<amrex::ParticleReal> wake_y_data;
        amrex::Gpu::TrackedVector<amrex::ParticleReal> wake_z_data;
    };

    struct RFCavity
    : public mixin::Named,
      public mixin::BeamOptic<RFCavity>,
      public mixin::LinearTransport<RFCavity>,
      public mixin::Thick,
      public mixin::Alignment,
      public mixin::NoFinalize,
      public mixin::PipeAperture,
      public mixin::SpinTransport,
      public amrex::simd::Vectorized<amrex::simd::native_simd_size_particlereal>
    {
        static constexpr auto type = "RFCavity";
        using PType = ImpactXParticleContainer::ParticleType;

        using DynamicData = mixin::GPUDataRegistry<CavityFourierCoefficients>;
        using WakeData = mixin::GPUDataRegistry<CavityWakeFunction>;

        static constexpr int DEFAULT_mapsteps = 10;

        /** An RF cavity
         *
         * @param ds Segment length in m
         * @param escale Scaling factor for on-axis RF electric field in 1/m
         * @param freq RF frequency in Hz
         * @param phase RF driven phase in deg
         * @param cos_coef TODO
         * @param sin_coef TODO
         * @param z_data Optional longitudinal positions of the single-particle wake function in m
         * @param wake_x_data Optional horizontal single-particle wake function data
         * @param wake_y_data Optional vertical single-particle wake function data
         * @param wake_z_data Optional longitudinal single-particle wake function data
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param aperture_x horizontal half-aperture in m
         * @param aperture_y vertical half-aperture in m
         * @param mapsteps number of integration steps per slice used for
         *        map and reference particle push in applied fields
         * @param nslice number of slices used for the application of space charge
         * @param name a user defined and not necessarily unique name of the element
         */
        RFCavity (
            amrex::ParticleReal ds,
            amrex::ParticleReal escale,
            amrex::ParticleReal freq,
            amrex::ParticleReal phase,
            std::vector<amrex::ParticleReal> cos_coef,
            std::vector<amrex::ParticleReal> sin_coef,
            std::optional<std::vector<amrex::ParticleReal>> z_data,
            std::optional<std::vector<amrex::ParticleReal>> wake_x_data,
            std::optional<std::vector<amrex::ParticleReal>> wake_y_data,
            std::optional<std::vector<amrex::ParticleReal>> wake_z_data,
            amrex::ParticleReal dx = DEFAULT_dx,
            amrex::ParticleReal dy = DEFAULT_dy,
            amrex::ParticleReal rotation_degree = DEFAULT_rotation_degree,
            amrex::ParticleReal aperture_x = DEFAULT_aperture_x,
            amrex::ParticleReal aperture_y = DEFAULT_aperture_y,
            int mapsteps = DEFAULT_mapsteps,
            int nslice = DEFAULT_nslice,
            std::optional<std::string> name = DEFAULT_name
        )
          : Named(std::move(name)),
            Thick(ds, nslice),
            Alignment(dx, dy, rotation_degree),
            PipeAperture(aperture_x, aperture_y),
            m_escale(escale), m_freq(freq), m_phase(phase), m_mapsteps(mapsteps),
            m_id(DynamicData::allocate_id()),
            m_wake_id(WakeData::allocate_id())
        {
            if (m_mapsteps <= 0)
                throw std::invalid_argument("RFCavity: mapsteps must be > 0!");

            m_ncoef = int(cos_coef.size());
            if (m_ncoef != int(sin_coef.size()))
                throw std::invalid_argument("RFCavity: cos and sin coefficients must have same length!");

            auto& coef = DynamicData::emplace(
                m_id,
                std::move(cos_coef),
                std::move(sin_coef)
            );
            m_cos_h_data = coef.cos.host_const().data();
            m_sin_h_data = coef.sin.host_const().data();

            // Optional user-specified wake data; absent data is stored as empty.
            WakeData::emplace(
                m_wake_id,
                std::move(z_data).value_or(std::vector<amrex::ParticleReal>{}),
                std::move(wake_x_data).value_or(std::vector<amrex::ParticleReal>{}),
                std::move(wake_y_data).value_or(std::vector<amrex::ParticleReal>{}),
                std::move(wake_z_data).value_or(std::vector<amrex::ParticleReal>{})
            );
        }

        /** Reverse the element (negate ds) */
        void reverse () {
            // Reversing ds traverses the Fourier profile in the opposite z direction,
            // so the odd sine terms change sign implicitly via sin(-kz) = -sin(kz).
            Thick::reverse();
        }

        /** Push all particles */
        using BeamOptic::operator();

        /** Compute and cache the constants for the push.
         *
         * In particular, used to pre-compute and cache variables that are
         * independent of the individually tracked particle.
         *
         * @param refpart reference particle
         */
        void compute_constants (RefPart const & refpart)
        {
            using namespace amrex::literals; // for _rt and _prt

            Alignment::compute_constants(refpart);

            auto const & coef = *DynamicData::get(m_id);
            m_cos_d_data = coef.cos.device_const().data();
            m_sin_d_data = coef.sin.device_const().data();
        }

        /** This is an RF cavity functor, so that a variable of this type can be used like
         *  an RF cavity function.
         *
         * The @see compute_constants method must be called before pushing particles through this operator.
         *
         * @param x particle position in x
         * @param y particle position in y
         * @param t particle position in t
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t
         * @param idcpu particle global index
         * @param refpart reference particle
         */
        template<typename T_Real=amrex::ParticleReal, typename T_IdCpu=uint64_t>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
            T_Real & AMREX_RESTRICT x,
            T_Real & AMREX_RESTRICT y,
            T_Real & AMREX_RESTRICT t,
            T_Real & AMREX_RESTRICT px,
            T_Real & AMREX_RESTRICT py,
            T_Real & AMREX_RESTRICT pt,
            [[maybe_unused]] T_IdCpu const & AMREX_RESTRICT idcpu,
            [[maybe_unused]] RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // get the linear map
            amrex::SmallMatrix<amrex::ParticleReal, 6, 6, amrex::Order::F, 1> const R = m_map;

            // symplectic linear map for the RF cavity is computed using the
            // Hamiltonian formalism as described in:
            // https://uspas.fnal.gov/materials/09UNM/ComputationalMethods.pdf.
            // R denotes the transfer matrix in the basis (x,px,y,py,t,pt),
            // so that, e.g., R(3,4) = dyf/dpyi.
            amrex::SmallVector<T_Real, 6, 1> const v{x, px, y, py, t, pt};

            // push particles using the linear map
            auto const out = R * v;

            // assign updated values
            x  = out[1];
            px = out[2];
            y  = out[3];
            py = out[4];
            t  = out[5];
            pt = out[6];
        }

        /** This pushes the reference particle.
         *
         * @param[in,out] refpart reference particle
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        void operator() (RefPart & AMREX_RESTRICT refpart) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using amrex::Math::powi;

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
            m_map = decltype(m_map)::Identity();

            // initialize the spin-orbit coupling matrix
            m_spin_coupling = {};

            // length of the current slice
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // compute initial value of beta*gamma
            amrex::ParticleReal const bgi = std::sqrt(powi<2>(pt) - 1.0_prt);

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
            amrex::ParticleReal const bgf = std::sqrt(powi<2>(ptf) - 1.0_prt);

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
                   m_map(i, j) = m_map(i, j) * scale_in / scale_fin;
               }
            }

            // convert spin-orbit coupling map from dynamic to static units
            for (int i=1; i<4; i++) {
                for (int j=2; j<7; j+=2) {
                    m_spin_coupling(i, j) *= bgi;
                }
            }

            // advance integrated path length
            refpart.s = s + slice_ds;
        }

        /** This operator pushes the particle spin.
         *
         * @param x particle position in x
         * @param y particle position in y
         * @param t particle position in t
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t
         * @param sx particle spin x-component
         * @param sy particle spin y-component
         * @param sz particle spin z-component
         * @param idcpu particle global index
         * @param refpart reference particle
         */
        template<typename T_Real=amrex::ParticleReal, typename T_IdCpu=uint64_t>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void spin_and_phasespace_push (
            T_Real & AMREX_RESTRICT x,
            T_Real & AMREX_RESTRICT y,
            T_Real & AMREX_RESTRICT t,
            T_Real & AMREX_RESTRICT px,
            T_Real & AMREX_RESTRICT py,
            T_Real & AMREX_RESTRICT pt,
            T_Real & AMREX_RESTRICT sx,
            T_Real & AMREX_RESTRICT sy,
            T_Real & AMREX_RESTRICT sz,
            T_IdCpu const & AMREX_RESTRICT idcpu,
            RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // initialize the three components of the axis-angle vector
            T_Real lambdax = 0_prt;
            T_Real lambday = 0_prt;
            T_Real lambdaz = 0_prt;

            // store the phase space variables in vector form
            amrex::SmallVector<T_Real, 6, 1> const v{x, px, y, py, t, pt};

            // get the spin-orbit coupling matrix
            amrex::SmallMatrix<amrex::ParticleReal, 3, 6, amrex::Order::F, 1> const A = m_spin_coupling;

            // use phase space variables to obtain the angle-axis generator of spin rotation
            auto const out = A * v;

            // update the angle-axis generator
            lambdax = out[1];
            lambday = out[2];
            lambdaz = out[3];

            // push the spin vector using the generator just determined
            rotate_spin(lambdax,lambday,lambdaz,sx,sy,sz);

            // phase space push
            (*this)(x, y, t, px, py, pt, idcpu, refpart);
        }

        /** This pushes the covariance matrix. */
        using LinearTransport::operator();

        /** This function returns the linear transport map.
         *
         * @returns 6x6 transport matrix
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        Map6x6
        transport_map ([[maybe_unused]] RefPart const & AMREX_RESTRICT refpart) const
        {

            Map6x6 R = Map6x6::Identity();
            R = m_map;

            // apply the transverse rotation (roll) alignment error
            return rotate_aligned_map(R);
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
            using amrex::Math::powi;

            // pick the right data depending if we are on the host side
            // (reference particle push) or device side (particles):
#if AMREX_DEVICE_COMPILE
            amrex::ParticleReal const * cos_data = m_cos_d_data;
            amrex::ParticleReal const * sin_data = m_sin_d_data;
#else
            amrex::ParticleReal const * cos_data = m_cos_h_data;
            amrex::ParticleReal const * sin_data = m_sin_h_data;
#endif

            // specify constants
            using ablastr::constant::math::pi;
            amrex::ParticleReal const zlen = std::abs(m_ds);
            amrex::ParticleReal const zmid = zlen * 0.5_prt;

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
                 efield = efield + cos_data[j] * std::cos(j*2*pi*z/zlen) +
                     sin_data[j] * std::sin(j*2*pi*z/zlen);
                 efieldp = efieldp-j*2*pi*cos_data[j] * std::sin(j*2*pi*z/zlen)/zlen +
                      j*2*pi*sin_data[j] * std::cos(j*2*pi*z/zlen)/zlen;
                 efieldpp = efieldpp- powi<2>(j*2*pi*cos_data[j]/zlen)  * std::cos(j*2*pi*z/zlen) -
                     powi<2>(j*2*pi*sin_data[j]/zlen)  * std::sin(j*2*pi*z/zlen);
                 efieldint = efieldint + zlen*cos_data[j] * std::sin(j*2*pi*z/zlen)/(j*2*pi) -
                      zlen*sin_data[j] * std::cos(j*2*pi*z/zlen)/(j*2*pi);
               }
            }
            else  // endpoint of the RF, outsize zlen
            {
               efieldint = std::copysign(zmid, z)*0.5_prt*cos_data[0];
               for (int j=1; j < m_ncoef; ++j)
               {
                 efieldint = efieldint - zlen*sin_data[j] * std::cos(j*pi)/(j*2*pi);
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
            using amrex::Math::powi;

            // push the reference particle
            amrex::ParticleReal const t = refpart.t;
            amrex::ParticleReal const pt = refpart.pt;

            if (pt < -1.0_prt) {
                refpart.t = t + tau/std::sqrt(1.0_prt - powi<-2>(pt));
                refpart.pt = pt;
            }
            else {
                refpart.t = t;
                refpart.pt = pt;
            }

            // push the linear map equations
            amrex::SmallMatrix<amrex::ParticleReal, 6, 6, amrex::Order::F, 1> const R = m_map;
            amrex::SmallMatrix<amrex::ParticleReal, 3, 6, amrex::Order::F, 1> const A = m_spin_coupling;
            amrex::SmallMatrix<amrex::ParticleReal, 3, 6, amrex::Order::F, 1> dA = {};
            amrex::ParticleReal const betgam = refpart.beta_gamma();

            m_map(5,5) = R(5,5) + tau*R(6,5)/powi<3>(betgam);
            m_map(5,6) = R(5,6) + tau*R(6,6)/powi<3>(betgam);

            // BELOW: if spin is needed only:
            // Define parameters and intermediate constants
            using ablastr::constant::math::pi;
            using ablastr::constant::SI::c;
            amrex::ParticleReal const k = (2_prt*pi/c)*m_freq;
            amrex::ParticleReal const phi = m_phase*(pi/180_prt);
            amrex::ParticleReal const E0 = m_escale;
            auto [ez, ezp, ezint] = RF_Efield(zeval);
            amrex::ignore_unused(ezint);

            // Update spin-orbit coupling matrix here
            amrex::ParticleReal G = refpart.gyromagnetic_anomaly;
            amrex::ParticleReal cos_term = std::cos(k*t+phi);
            amrex::ParticleReal sin_term = std::sin(k*t+phi);
            dA(1,3) = E0 * tau / 2_prt * (k * ez * (G*pt - 1_prt) * sin_term / std::sqrt(1_prt + powi<2>(pt)) + cos_term * (G + 1_prt/(1_prt - pt)) * ezp);
            dA(1,4) = E0 * tau / std::sqrt(1_prt + powi<2>(pt)) * cos_term * (G + 1_prt/(1_prt - pt));
            dA(2,1) = -dA(1,3);
            dA(2,2) = -dA(1,4);

            // update the spin-orbit coupling matrix here
            m_spin_coupling = A + dA*R;
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
            using amrex::Math::powi;

            amrex::ParticleReal const t = refpart.t;
            amrex::ParticleReal const pt = refpart.pt;

            // Define parameters and intermediate constants
            using ablastr::constant::math::pi;
            using ablastr::constant::SI::c;
            amrex::ParticleReal const k = (2_prt*pi/c)*m_freq;
            amrex::ParticleReal const phi = m_phase*(pi/180_prt);
            amrex::ParticleReal const E0 = m_escale;

            // push the reference particle
            auto [ez, ezp, ezint] = RF_Efield(zeval);
            amrex::ignore_unused(ez, ezint);

            refpart.t = t;
            refpart.pt = pt;

            // push the linear map equations
            amrex::SmallMatrix<amrex::ParticleReal, 6, 6, amrex::Order::F, 1> const R = m_map;
            amrex::ParticleReal const s = tau/refpart.beta_gamma();
            amrex::ParticleReal const L = E0*ezp * std::sin(k*t+phi)/(2_prt*k);

            m_map(1,1) = (1_prt-s*L)*R(1,1) + s*R(2,1);
            m_map(1,2) = (1_prt-s*L)*R(1,2) + s*R(2,2);
            m_map(2,1) = -s * powi<2>(L)*R(1,1) + (1_prt+s*L)*R(2,1);
            m_map(2,2) = -s * powi<2>(L)*R(1,2) + (1_prt+s*L)*R(2,2);

            m_map(3,3) = (1_prt-s*L)*R(3,3) + s*R(4,3);
            m_map(3,4) = (1_prt-s*L)*R(3,4) + s*R(4,4);
            m_map(4,3) = -s * powi<2>(L)*R(3,3) + (1_prt+s*L)*R(4,3);
            m_map(4,4) = -s * powi<2>(L)*R(3,4) + (1_prt+s*L)*R(4,4);

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
            amrex::ParticleReal const k = (2_prt*pi/c)*m_freq;
            amrex::ParticleReal const phi = m_phase*(pi/180_prt);
            amrex::ParticleReal const E0 = m_escale;

            // push the reference particle
            auto [ez, ezp, ezint] = RF_Efield(z);
            amrex::ignore_unused(ez);
            zeval = z + tau;
            auto [ezf, ezpf, ezintf] = RF_Efield(zeval);
            amrex::ignore_unused(ezf);

            refpart.t = t;
            refpart.pt = pt - E0*(ezintf-ezint) * std::cos(k*t+phi);

            // push the linear map equations
            amrex::SmallMatrix<amrex::ParticleReal, 6, 6, amrex::Order::F, 1> const R = m_map;
            amrex::ParticleReal const M = E0*(ezintf-ezint)*k * std::sin(k*t+phi);
            amrex::ParticleReal const L = E0*(ezpf-ezp) * std::sin(k*t+phi)/(2_prt*k)+M*0.5_prt;

            m_map(2,1) = L*R(1,1) + R(2,1);
            m_map(2,2) = L*R(1,2) + R(2,2);

            m_map(4,3) = L*R(3,3) + R(4,3);
            m_map(4,4) = L*R(3,4) + R(4,4);

            m_map(6,5) = M*R(5,5) + R(6,5);
            m_map(6,6) = M*R(5,6) + R(6,6);

        }

        amrex::ParticleReal m_escale; //! scaling factor for RF electric field
        amrex::ParticleReal m_freq; //! RF frequency in Hz
        amrex::ParticleReal m_phase; //! RF driven phase in deg
        int m_mapsteps; //! number of map integration steps per slice
        int m_id; //! unique RF cavity id used for data lookup map
        int m_wake_id; //! unique id used for user-specified wake data lookup

        int m_ncoef = 0; //! number of Fourier coefficients
        amrex::ParticleReal const * m_cos_h_data = nullptr; //! non-owning pointer to host cosine coefficients
        amrex::ParticleReal const * m_sin_h_data = nullptr; //! non-owning pointer to host sine coefficients
        amrex::ParticleReal const * m_cos_d_data = nullptr; //! non-owning pointer to device cosine coefficients
        amrex::ParticleReal const * m_sin_d_data = nullptr; //! non-owning pointer to device sine coefficients

        // Reference-trajectory linearization around the reference particle.
        // Computed during the reference-particle push and consumed during the
        // particle push.
        // mutable: written by the const reference push before the element is copied
        // into the particle-push functor.
        mutable amrex::SmallMatrix<amrex::ParticleReal, 6, 6, amrex::Order::F, 1> m_map{}; //! linearized map
        mutable amrex::SmallMatrix<amrex::ParticleReal, 3, 6, amrex::Order::F, 1> m_spin_coupling{}; //! linearized spin-orbit coupling matrix
    };

} // namespace impactx

IMPACTX_GPUDATA_EXTERN(impactx::elements::RFCavity)
// second registry of this element: the user-specified wake data.
// IMPACTX_GPUDATA_EXTERN only covers `ElementType::DynamicData`, so spell this one out.
extern template struct impactx::elements::mixin::GPUDataRegistry<impactx::elements::CavityWakeFunction>;
IMPACTX_PUSH_EXTERN_TEMPLATE(impactx::elements::RFCavity)

#endif // IMPACTX_RFCAVITY_H
