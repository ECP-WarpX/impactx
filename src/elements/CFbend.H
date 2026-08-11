/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_CFBEND_H
#define IMPACTX_CFBEND_H

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
    struct CFbendGPU
    : public mixin::BeamOptic<CFbendGPU>,
      public mixin::LinearTransport<CFbendGPU>,
      public mixin::Thick,
      public mixin::Alignment,
      public mixin::PipeAperture,
      public mixin::SpinTransport,
      public mixin::NoFinalize
      // At least on Intel AVX512, there is a small overhead to vectorize this element, see
      // https://github.com/BLAST-ImpactX/impactx/pull/1002
      // public amrex::simd::Vectorized<amrex::simd::native_simd_size_particlereal>
    {
        static constexpr auto type = "CFbend";
        using PType = ImpactXParticleContainer::ParticleType;

        /** An combined-function bend, consisting of an ideal sector bend with
         *  an upright quadrupole focusing component.
         *
         * @param ds Segment length in m.
         * @param rc Radius of curvature in m.
         * @param k  Quadrupole strength in m^(-2) (MADX convention)
         *           = (gradient in T/m) / (rigidity in T-m)
         *           k > 0 horizontal focusing
         *           k < 0 horizontal defocusing
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param aperture_x horizontal half-aperture in m
         * @param aperture_y vertical half-aperture in m
         * @param nslice number of slices used for the application of space charge
         */
        CFbendGPU (
            amrex::ParticleReal ds,
            amrex::ParticleReal rc,
            amrex::ParticleReal k,
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
          m_rc(rc), m_k(k)
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
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // find beta*gamma^2, beta
            amrex::ParticleReal const betgam2 = powi<2>(refpart.pt) - 1_prt;
            amrex::ParticleReal const bet = refpart.beta();
            amrex::ParticleReal const ibetgam2 = 1_prt / betgam2;
            amrex::ParticleReal const b2rc2 = powi<2>(bet) * powi<2>(m_rc);

            // update horizontal and longitudinal phase space variables
            amrex::ParticleReal const gx = m_k + powi<-2>(m_rc);
            amrex::ParticleReal const omega_x = std::sqrt(std::abs(gx));

            // update vertical phase space variables
            amrex::ParticleReal const gy = -m_k;
            amrex::ParticleReal const omega_y = std::sqrt(std::abs(gy));

            // trigonometry
            amrex::ParticleReal const cx = gx > 0_prt ? std::cos(omega_x * slice_ds) : std::cosh(omega_x * slice_ds);
            amrex::ParticleReal const sx = gx > 0_prt ? std::sin(omega_x * slice_ds)/omega_x : std::sinh(omega_x * slice_ds)/omega_x;
            amrex::ParticleReal const cy = gy > 0_prt ? std::cos(omega_y * slice_ds) : std::cosh(omega_y * slice_ds);
            amrex::ParticleReal const sy = gy > 0_prt ? std::sin(omega_y * slice_ds)/omega_y : std::sinh(omega_y * slice_ds)/omega_y;

            amrex::ParticleReal const ibrc = 1_prt / (bet * m_rc);
            amrex::ParticleReal const igbrc = 1_prt / (gx * bet * m_rc);
            amrex::ParticleReal const igb2r2 = 1_prt / (gx * b2rc2);

            // elements of the linear transport matrix
            // The matrix elements of R are identical with those found in:
            // F. C. Iselin, Part. Accel. 17, pp. 143-155 (1985), Section 4,
            // except the longitudinal variables used here (t,pt) have the opposite sign.
            m_R11 = cx;
            m_R12 = sx;
            m_R16 = -(1_prt - cx) * igbrc;
            m_R21 = -gx * sx;
            m_R22 = cx;
            m_R26 = -sx * ibrc;
            m_R33 = cy;
            m_R34 = sy;
            m_R43 = -gy * sy;
            m_R44 = cy;
            m_R51 = -m_R26;
            m_R52 = -m_R16;
            m_R56 = slice_ds * ibetgam2 + (sx - slice_ds) * igb2r2;

            // access reference particle values for spin calculation
            amrex::ParticleReal G = refpart.gyromagnetic_anomaly;
            amrex::ParticleReal const gamma = refpart.gamma();
            amrex::ParticleReal const beta = refpart.beta();
            amrex::ParticleReal const gyro_const1 = 1_prt + G * gamma;
            amrex::ParticleReal const gyro_const2 = 1_prt + G * gamma * gamma;

            // trigonometry for spin calculation
            amrex::ParticleReal const h = (m_rc == 0_prt) ? 0_prt : 1_prt/m_rc;
            auto const [sinG, cosG] = amrex::Math::sincos(G * h * gamma * slice_ds);

            // intermediate constants for spin calculation
            amrex::ParticleReal const factor = 1_prt / (-gy + powi<2>(G * h * gamma));

            // nonvanishing components of spin rotation at the design point
            m_lambday = -G * h * gamma * slice_ds;

            // elements of the spin-orbit coupling matrix
            m_A13 = (G*h*gyro_const2*gy*cy*sinG + gy*(powi<2>(G*h)*(gamma - 1_prt)*gamma - gy*gyro_const1 )*sy*cosG) * factor;
            m_A14 = (G*h*gyro_const2*gy*sy*sinG + (powi<2>(G*h)*(gamma - 1_prt)*gamma - gy*gyro_const1)*(1_prt - cosG*cy))*factor;
            m_A21 = -gyro_const1 * gx * sx;
            m_A22 = -gyro_const1 * (1_prt - cx);
            m_A26 = G*h*beta*gamma*slice_ds - h*gyro_const1*sx/(beta);
            m_A33 = (-(-gy*gyro_const1 + powi<2>(G*h)*gamma*(gamma - 1_prt))*gy*sy*sinG -gy*G*h*gyro_const2*(1_prt - cosG*cy))*factor;
            m_A34 = ((-gy*gyro_const1 + powi<2>(G*h)*gamma*(gamma - 1_prt)) * cy*sinG + G*h*gyro_const2*gy*cosG*sy)*factor;
        }

        /** This is a cfbend functor, so that a variable of this type can be used like a cfbend function.
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
         * @param refpart reference particle (unused)
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
            [[maybe_unused]] RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // initialize output values
            T_Real xout = x;
            T_Real yout = y;
            T_Real tout = t;

            // initialize output values of momenta
            T_Real pxout = px;
            T_Real pyout = py;
            // T_Real ptout = pt;

            // update horizontal and longitudinal phase space variables
            //   advance position and momentum: (de)focusing
            x =     m_R11 * xout + m_R12 * px + m_R16 * pt;
            pxout = m_R21 * xout + m_R22 * px + m_R26 * pt;
            t =     m_R51 * xout + m_R52 * px + m_R56 * pt + tout;
            // ptout = pt;

            // update vertical phase space variables
            //   advance position and momentum (de)focusing
            y =     m_R33 * yout + m_R34 * py;
            pyout = m_R43 * yout + m_R44 * py;

            // assign updated momenta
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
            amrex::ParticleReal const theta = slice_ds/m_rc;
            amrex::ParticleReal const B = std::sqrt(powi<2>(pt)-1.0_prt)/m_rc;

            // calculate expensive terms once
            auto const [sin_theta, cos_theta] = amrex::Math::sincos(theta);

            // advance position and momentum (bend)
            refpart.px = px*cos_theta - pz*sin_theta;
            refpart.py = py;
            refpart.pz = pz*cos_theta + px*sin_theta;
            refpart.pt = pt;

            refpart.x = x + (refpart.pz - pz)/B;
            refpart.y = y + (theta/B)*py;
            refpart.z = z - (refpart.px - px)/B;
            refpart.t = t - (theta/B)*pt;

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
         * @param refpart reference particle (unused)
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
            using amrex::Math::powi;

            // initialize the three components of the axis-angle vector
            T_Real lambdax = 0_prt;
            T_Real lambday = 0_prt;
            T_Real lambdaz = 0_prt;

            // set the angle-axis generator based on the phase space variables
            lambdax = m_A13 * y + m_A14 * py;
            lambday = m_A21 * x + m_A22 * px + m_A26 * pt;
            lambdaz = m_A33 * y + m_A34 * py;

            // push the spin vector using the generator just determined
            rotate_spin(lambdax,lambday,lambdaz,sx,sy,sz);

            // axis-angle vector components generating the reference spin map
            lambdax = 0.0_prt;
            lambday = m_lambday;
            lambdaz = 0.0_prt;

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
        transport_map ([[maybe_unused]] RefPart const & AMREX_RESTRICT refpart) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using amrex::Math::powi;

            // length of the current slice
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // find beta*gamma^2, beta
            amrex::ParticleReal const betgam2 = powi<2>(refpart.pt) - 1_prt;
            amrex::ParticleReal const bet = refpart.beta();
            amrex::ParticleReal const ibetgam2 = 1_prt / betgam2;
            amrex::ParticleReal const b2rc2 = powi<2>(bet) * powi<2>(m_rc);

            // update horizontal and longitudinal phase space variables
            amrex::ParticleReal const gx = m_k + powi<-2>(m_rc);
            amrex::ParticleReal const omega_x = std::sqrt(std::abs(gx));

            // update vertical phase space variables
            amrex::ParticleReal const gy = -m_k;
            amrex::ParticleReal const omega_y = std::sqrt(std::abs(gy));

            // trigonometry
            auto const [sinx, cosx] = amrex::Math::sincos(omega_x * slice_ds);
            amrex::ParticleReal const sinhx = std::sinh(omega_x * slice_ds);
            amrex::ParticleReal const coshx = std::cosh(omega_x * slice_ds);
            auto const [siny, cosy] = amrex::Math::sincos(omega_y * slice_ds);
            amrex::ParticleReal const sinhy = std::sinh(omega_y * slice_ds);
            amrex::ParticleReal const coshy = std::cosh(omega_y * slice_ds);

            amrex::ParticleReal const igbrc = 1_prt / (gx * bet * m_rc);
            amrex::ParticleReal const iobrc = 1_prt / (omega_x * bet * m_rc);
            amrex::ParticleReal const igobr = 1_prt / (gx * omega_x * b2rc2);

            // initialize linear map matrix elements
            Map6x6 R = Map6x6::Identity();

            R(1,1) = gx > 0_prt ? cosx : coshx;
            R(1,2) = gx > 0_prt ? sinx / omega_x : sinhx / omega_x;
            R(1,6) = gx > 0_prt ? -(1_prt - cosx) * igbrc : -(1_prt - coshx) * igbrc;
            R(2,1) = gx > 0_prt ? -omega_x * sinx : omega_x * sinhx;
            R(2,2) = gx > 0_prt ? cosx : coshx;
            R(2,6) = gx > 0_prt ? -sinx * iobrc : -sinhx * iobrc;
            R(3,3) = gy > 0_prt ? cosy : coshy;
            R(3,4) = gy > 0_prt ? siny / omega_y : sinhy / omega_y;
            R(4,3) = gy > 0_prt ? -omega_y * siny : omega_y * sinhy;
            R(4,4) = gy > 0_prt ? cosy : coshy;
            R(5,1) = gx > 0_prt ? sinx * iobrc : sinhx * iobrc;
            R(5,2) = gx > 0_prt ? (1_prt - cosx) * igbrc : (1_prt - coshx) * igbrc;
            R(5,6) = gx > 0_prt ?
                slice_ds * ibetgam2 + (sinx  - omega_x * slice_ds) * igobr :
                slice_ds * ibetgam2 + (sinhx - omega_x * slice_ds) * igobr;

            // apply the transverse rotation (roll) alignment error
            return rotate_aligned_map(R);
        }

        amrex::ParticleReal m_rc; //! bend radius in m
        amrex::ParticleReal m_k;  //! quadrupole strength in m^(-2)

    private:
        // constants that are independent of the individually tracked particle,
        // see: compute_constants() to refresh
        amrex::ParticleReal m_R11, m_R12, m_R16, m_R21, m_R22, m_R26, m_R33, m_R34, m_R43, m_R44, m_R51, m_R52, m_R56;
        amrex::ParticleReal m_lambday, m_A13, m_A14, m_A21, m_A22, m_A26, m_A33, m_A34;
    };

    /** CFbend: GPU state plus common host-only metadata */
    using CFbend = Element<CFbendGPU>;

} // namespace impactx

IMPACTX_PUSH_EXTERN_TEMPLATE(impactx::elements::CFbendGPU)

#endif // IMPACTX_CFBEND_H
