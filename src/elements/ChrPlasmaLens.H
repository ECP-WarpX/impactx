/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl, Kyrre Sjobak
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_CHRPLASMALENS_H
#define IMPACTX_CHRPLASMALENS_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/gpuelement.H"
#include "mixin/alignment.H"
#include "mixin/pipeaperture.H"
#include "mixin/beamoptic.H"
#include "mixin/thick.H"
#include "mixin/lineartransport.H"
#include "mixin/spintransport.H"
#include "mixin/nofinalize.H"

#include <AMReX_Extension.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>
#include <AMReX_SIMD.H>

#include <cmath>
#include <stdexcept>

namespace impactx::elements
{
    struct ChrPlasmaLensGPU
    : public mixin::BeamOptic<ChrPlasmaLensGPU>,
      public mixin::LinearTransport<ChrPlasmaLensGPU>,
      public mixin::Thick,
      public mixin::Alignment,
      public mixin::PipeAperture,
      public mixin::SpinTransport,
      public mixin::NoFinalize,
      public amrex::simd::Vectorized<amrex::simd::native_simd_size_particlereal>
    {
        static constexpr auto type = "ChrPlasmaLens";
        using PType = ImpactXParticleContainer::ParticleType;

        static constexpr int DEFAULT_unit = 0;

        /** An active cylindrically symmetric plasma lens with chromatic focusing
         *
         * The Hamiltonian is expanded through second order in the
         * transverse variables (x,px,y,py), with the exact pt dependence retained.
         *
         * @param ds Segment length in m.
         * @param k  plasma lens focusing strength in m^(-2)
         *           = (azimuthal field gradient in T/m) / (rigidity in T-m)
         *        OR azimuthal magnetic field gradient in T/m
         * @param unit  Unit specification
         *           unit = 0 focusing strength in m^(-2)
         *           unit = 1 focusing strength in T/m
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param aperture_x horizontal half-aperture in m
         * @param aperture_y vertical half-aperture in m
         * @param nslice number of slices used for the application of space charge
         */
        ChrPlasmaLensGPU (
            amrex::ParticleReal ds,
            amrex::ParticleReal k,
            int unit = DEFAULT_unit,
            amrex::ParticleReal dx = DEFAULT_dx,
            amrex::ParticleReal dy = DEFAULT_dy,
            amrex::ParticleReal rotation_degree = DEFAULT_rotation_degree,
            amrex::ParticleReal aperture_x = DEFAULT_aperture_x,
            amrex::ParticleReal aperture_y = DEFAULT_aperture_y,
            int nslice = DEFAULT_nslice
        )
        : Thick(ds, nslice),
          Alignment(dx, dy, rotation_degree),
          PipeAperture(aperture_x, aperture_y),
          m_k(k), m_unit(unit)
        {
        }

        /** Reverse the element (negate ds) */
        void reverse () { Thick::reverse(); }

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
            using amrex::Math::powi;

            Alignment::compute_constants(refpart);

            // length of the current slice
            m_slice_ds = m_ds / nslice();

            // access reference particle values to find beta and gamma
            m_beta = refpart.beta();
            amrex::ParticleReal const gamma = refpart.gamma();

            // normalize focusing strength units to MAD-X convention if needed
            m_g = m_unit == 1 ? m_k / refpart.rigidity_Tm() : m_k;

            //For m_g=0 time propagation
            m_const1 = 1_prt / (2_prt * powi<3>(m_beta) * powi<2>(gamma));
        }

        /** This is a chrplasmalens functor, so that a variable of this type can be used like a chrplasmalens function.
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
            T_Real const & AMREX_RESTRICT pt,
            [[maybe_unused]] T_IdCpu const & AMREX_RESTRICT idcpu,
            RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using amrex::Math::powi;
            using namespace std; // for cmath(float)

            // compute particle momentum deviation delta + 1
            T_Real const delta1 = sqrt(1_prt - 2_prt*pt/m_beta + powi<2>(pt));
            T_Real const delta = delta1 - 1_prt;

            // compute phase advance per unit length in s (in rad/m)
            // chromatic dependence on delta is included
            T_Real const omega = sqrt(abs(m_g)/delta1);

            // initialize output values
            T_Real xout = x;
            T_Real yout = y;
            T_Real tout = t;

            // initialize output values of momenta
            T_Real pxout = px;
            T_Real pyout = py;
            // T_Real const ptout = pt;

            // placeholder variables -- canonical position and momentum
            T_Real const q1 = x;
            T_Real const q2 = y;
            T_Real const p1 = px;
            T_Real const p2 = py;

            if (m_g > 0_prt)
            {
                auto const [sin_ods, cos_ods] = amrex::Math::sincos(omega*m_slice_ds);

                // advance transverse position and momentum (focusing)
                xout = cos_ods * x + sin_ods / (omega * delta1) * px;
                pxout = -omega * delta1 * sin_ods * x + cos_ods * px;

                yout = cos_ods * y + sin_ods / (omega * delta1) * py;
                pyout = -omega * delta1 * sin_ods * y + cos_ods * py;

                // the corresponding symplectic update to t (focusing both x and y)
                T_Real const term = pt + delta/m_beta;
                T_Real const t0 = t - term*m_slice_ds/delta1;

                T_Real const w = omega*delta1;
                T_Real const term1 = -(powi<2>(p2) - powi<2>(q2) * powi<2>(w)) * sin(2_prt*m_slice_ds*omega);
                T_Real const term2 = -(powi<2>(p1) - powi<2>(q1) * powi<2>(w)) * sin(2_prt*m_slice_ds*omega);
                T_Real const term3 = -2_prt * q2 * p2 * w * cos(2_prt*m_slice_ds*omega);
                T_Real const term4 = -2_prt * q1 * p1 * w * cos(2_prt*m_slice_ds*omega);
                T_Real const term5 = 2_prt * omega * (q1*p1*delta1 + q2*p2*delta1
                    -(powi<2>(p1) + powi<2>(p2))*m_slice_ds - (powi<2>(q1) + powi<2>(q2)) * powi<2>(w)*m_slice_ds);
                tout = t0 + (-1_prt+m_beta*pt)
                          / (8_prt*m_beta * powi<3>(delta1)*omega)
                          * (term1+term2+term3+term4+term5);

            }
            else if (m_g < 0_prt)
            {
                auto const sinh_ods = sinh(omega*m_slice_ds);
                auto const cosh_ods = cosh(omega*m_slice_ds);

                // advance transverse position and momentum (defocusing)
                xout = cosh_ods * x + sinh_ods / (omega * delta1) * px;
                pxout = +omega * delta1 * sinh_ods * x + cosh_ods * px;

                yout = cosh_ods * y + sinh_ods / (omega * delta1) * py;
                pyout = +omega * delta1 * sinh_ods * y + cosh_ods * py;

                // the corresponding symplectic update to t (defocusing both x and y)
                T_Real const term = pt + delta/m_beta;
                T_Real const t0 = t - term*m_slice_ds/delta1;

                T_Real const w = omega*delta1;
                T_Real const term1 = -(powi<2>(p2) + powi<2>(q2) * powi<2>(w)) * sinh(2_prt*m_slice_ds*omega);
                T_Real const term2 = -(powi<2>(p1) + powi<2>(q1) * powi<2>(w)) * sinh(2_prt*m_slice_ds*omega);
                T_Real const term3 = -2_prt * q2 * p2 * w * cosh(2_prt*m_slice_ds*omega);
                T_Real const term4 = -2_prt * q1 * p1 * w * cosh(2_prt*m_slice_ds*omega);
                T_Real const term5 = 2_prt * omega * (q1*p1*delta1 + q2*p2*delta1
                    -(powi<2>(p1) + powi<2>(p2))*m_slice_ds - (powi<2>(q1) + powi<2>(q2)) * powi<2>(w)*m_slice_ds);
                tout = t0 + (-1_prt+m_beta*pt)
                          / (8_prt*m_beta * powi<3>(delta1)*omega)
                          * (term1+term2+term3+term4+term5);
            }
            else {
                xout = x + px*m_slice_ds / delta1;
                pxout = px;

                yout = y + py*m_slice_ds / delta1;
                pyout = py;

                // the corresponding symplectic update to t (zero strength = drift),
                // avoiding division by 0
                T_Real term = 2_prt * powi<2>(pt) + powi<2>(px) + powi<2>(py);
                term = 2_prt - 4_prt * m_beta * pt + powi<2>(m_beta) * term;
                term = -2_prt + powi<2>(refpart.gamma())*term;
                term = (-1_prt + m_beta * pt) * term;
                term = term * m_const1;
                tout = t - m_slice_ds * (1_prt / m_beta + term / powi<3>(delta1));
            }
            // advance longitudinal position and momentum

            // assign updated position & momenta
            x = xout;
            y = yout;
            t = tout;
            px = pxout;
            py = pyout;
            // pt = ptout;
        }

        /** This pushes the reference particle.
         *
         * @param[in,out] refpart reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
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
            amrex::ParticleReal const t = refpart.t;
            amrex::ParticleReal const pt = refpart.pt;
            amrex::ParticleReal const s = refpart.s;

            // length of the current slice
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // assign intermediate parameter
            amrex::ParticleReal const step = slice_ds / std::sqrt(powi<2>(pt)-1.0_prt);

            // advance position and momentum (straight element)
            refpart.x = x + step*px;
            refpart.y = y + step*py;
            refpart.z = z + step*pz;
            refpart.t = t - step*pt;

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
            T_Real const & AMREX_RESTRICT pt,
            T_Real & AMREX_RESTRICT sx,
            T_Real & AMREX_RESTRICT sy,
            T_Real & AMREX_RESTRICT sz,
            T_IdCpu const & AMREX_RESTRICT idcpu,
            RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using namespace std; // for cmath(float)
            using amrex::Math::powi;

            // initialize the three components of the axis-angle vector
            T_Real lambdax = 0_prt;
            T_Real lambday = 0_prt;
            T_Real lambdaz = 0_prt;

            // compute particle momentum deviation delta + 1
            T_Real const delta1 = sqrt(1_prt - 2_prt*pt/m_beta + powi<2>(pt));
            T_Real const inv_delta1 = 1_prt / delta1;

            // compute particle relativistic gamma
            T_Real const gamma = refpart.gamma() * (1_prt - m_beta*pt);
            T_Real const gyro_const = 1_prt + refpart.gyromagnetic_anomaly * gamma;

            // compute phase advance per unit length in s (in rad/m)
            // chromatic dependence on delta is included
            T_Real const omega = sqrt(abs(m_g)*inv_delta1);

            // compute trigonometric quantities
            auto const [sin_omega_ds, cos_omega_ds] = amrex::Math::sincos(omega*m_slice_ds);
            T_Real const sinh_omega_ds = sinh(omega*m_slice_ds);
            T_Real const cosh_omega_ds = cosh(omega*m_slice_ds);

            if (m_g > 0.0_prt)
            {

                // horizontally focusing quad case
                lambdax = +gyro_const * ( py*inv_delta1*(1_prt - cos_omega_ds) + y*omega*sin_omega_ds );
                lambday = -gyro_const * ( px*inv_delta1*(1_prt - cos_omega_ds) + x*omega*sin_omega_ds );

            } else if (m_g < 0.0_prt)
            {

                // horizontally defocusing quad case
                lambdax = -gyro_const * ( py*inv_delta1*(cosh_omega_ds - 1_prt) + y*omega*sinh_omega_ds );
                lambday = +gyro_const * ( px*inv_delta1*(cosh_omega_ds - 1_prt) + x*omega*sinh_omega_ds );

            } else {
                // treat as a drift
            }

            // push the spin vector using the generator just determined
            rotate_spin(lambdax,lambday,lambdaz,sx,sy,sz);

            // phase space push
            (*this)(x, y, t, px, py, pt, idcpu, refpart);

        }


        /** This pushes the covariance matrix. */
        using LinearTransport::operator();

        /** This function returns the linear transport map.
         *
         * @param[in] refpart reference particle
         * @returns 6x6 transport matrix
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        Map6x6
        transport_map (RefPart const & AMREX_RESTRICT refpart) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using amrex::Math::powi;

            // length of the current slice
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // access reference particle values to find beta*gamma^2
            amrex::ParticleReal const pt_ref = refpart.pt;
            amrex::ParticleReal const betgam2 = powi<2>(pt_ref) - 1_prt;

            // normalize focusing strength units to MAD-X convention if needed
            amrex::ParticleReal const g = m_unit == 1 ? m_k / refpart.rigidity_Tm() : m_k;

            // compute phase advance per unit length in s (in rad/m)
            amrex::ParticleReal const omega = std::sqrt(std::abs(g));

            // initialize linear map matrix elements
            Map6x6 R = Map6x6::Identity();

            if (g > 0.0) {
                R(1,1) = std::cos(omega*slice_ds);
                R(1,2) = std::sin(omega*slice_ds)/omega;
                R(2,1) = -omega*std::sin(omega*slice_ds);
                R(2,2) = std::cos(omega*slice_ds);
                R(3,3) = std::cos(omega*slice_ds);
                R(3,4) = std::sin(omega*slice_ds)/omega;
                R(4,3) = -omega*std::sin(omega*slice_ds);
                R(4,4) = std::cos(omega*slice_ds);
                R(5,6) = slice_ds/betgam2;
            } else if (g < 0.0) {
                R(1,1) = std::cosh(omega*slice_ds);
                R(1,2) = std::sinh(omega*slice_ds)/omega;
                R(2,1) = omega*std::sinh(omega*slice_ds);
                R(2,2) = std::cosh(omega*slice_ds);
                R(3,3) = std::cosh(omega*slice_ds);
                R(3,4) = std::sinh(omega*slice_ds)/omega;
                R(4,3) = omega*std::sinh(omega*slice_ds);
                R(4,4) = std::cosh(omega*slice_ds);
                R(5,6) = slice_ds/betgam2;
            } else {
                R(1,2) = slice_ds;
                R(3,4) = slice_ds;
                R(5,6) = slice_ds / betgam2;
            }

            // apply the transverse rotation (roll) alignment error
            return rotate_aligned_map(R);
        }

        amrex::ParticleReal m_k; //! focusing strength in 1/m^2 (or T/m)
        int m_unit; //! unit specification for focusing strength

      private:
        // constants that are independent of the individually tracked particle,
        // see: compute_constants() to refresh
        amrex::ParticleReal m_slice_ds;  //! m_ds / nslice();
        amrex::ParticleReal m_beta;
        amrex::ParticleReal m_g;
        amrex::ParticleReal m_const1;
    };

    /** ChrPlasmaLens: GPU state plus common host-only metadata */
    using ChrPlasmaLens = Element<ChrPlasmaLensGPU>;

} // namespace impactx

IMPACTX_PUSH_EXTERN_TEMPLATE(impactx::elements::ChrPlasmaLensGPU)

#endif // IMPACTX_CHRPLASMALENS_H
