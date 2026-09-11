/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_SOL_H
#define IMPACTX_SOL_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/gpuelement.H"
#include "mixin/alignment.H"
#include "mixin/pipeaperture.H"
#include "mixin/beamoptic.H"
#include "mixin/thick.H"
#include "mixin/lineartransport.H"
#include "mixin/nofinalize.H"
#include "mixin/spintransport.H"

#include <AMReX_Extension.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>
#include <AMReX_SIMD.H>

#include <cmath>
#include <stdexcept>

namespace impactx::elements
{
    struct SolGPU
    : public mixin::BeamOptic<SolGPU>,
      public mixin::LinearTransport<SolGPU>,
      public mixin::Thick,
      public mixin::Alignment,
      public mixin::PipeAperture,
      public mixin::SpinTransport,
      public mixin::NoFinalize
      // At least on Intel AVX512, there is a small overhead to vectorize this element, see
      // https://github.com/BLAST-ImpactX/impactx/pull/1002
      // public amrex::simd::Vectorized<amrex::simd::native_simd_size_particlereal>
    {
        static constexpr auto type = "Sol";
        using PType = ImpactXParticleContainer::ParticleType;

        /** An ideal hard-edge Solenoid magnet
         *
         * @param ds Segment length in m.
         * @param ks Solenoid strength in m^(-1) (MADX convention)
         *           = (magnetic field Bz in T) / (rigidity in T-m)
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param aperture_x horizontal half-aperture in m
         * @param aperture_y vertical half-aperture in m
         * @param nslice number of slices used for the application of space charge
         */
        SolGPU (
            amrex::ParticleReal ds,
            amrex::ParticleReal ks,
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
          m_ks(ks)
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

            // access reference particle values to find beta*gamma^2
            amrex::ParticleReal const betgam2 = powi<2>(refpart.pt) - 1.0_prt;
            m_ibeta = 1.0_prt / refpart.beta();

            // compute phase advance per unit length (in rad/m) and
            // rotation angle (in rad)
            amrex::ParticleReal const alpha = m_ks * 0.5_prt;
            amrex::ParticleReal const theta = alpha*slice_ds;
            auto const [sin_theta, cos_theta] = amrex::Math::sincos(theta);

            m_sin_theta = sin_theta;
            m_cos_theta = cos_theta;
            m_sin_theta_ialpha = sin_theta / alpha;
            m_neg_alpha_sin_theta = -alpha * sin_theta;
            m_ds_bg2 = slice_ds / betgam2;

            // access reference particle values required for spin push
            amrex::ParticleReal G = refpart.gyromagnetic_anomaly;
            amrex::ParticleReal const gamma = refpart.gamma();
            m_rel_const = 0.5_prt * (gamma - 1_prt);
            auto const [sin_G, cos_G] = amrex::Math::sincos(2_prt * G * theta);
            m_sin_G = sin_G;
            m_1cos_G = 1_prt - cos_G;
            m_neg2_theta_1pG = -2_prt * theta * (1_prt + G);

        }

        /** This is a sol functor, so that a variable of this type can be used like a sol function.
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
            T_Real pxout = px;
            T_Real pyout = py;
            // T_Real ptout = pt;

            // advance positions and momenta using map for focusing
            xout = m_cos_theta        * x
                 + m_sin_theta_ialpha * px;
            pxout = m_neg_alpha_sin_theta * x
                  + m_cos_theta           * px;

            yout = m_cos_theta        * y
                 + m_sin_theta_ialpha * py;
            pyout = m_neg_alpha_sin_theta * y
                  + m_cos_theta           * py;

            tout = t + m_ds_bg2 * pt;
            // ptout = pt;

            // assign intermediate momenta
            px = pxout;
            py = pyout;
            // pt = ptout;

            // advance positions and momenta using map for rotation
            x = m_cos_theta * xout
              + m_sin_theta * yout;
            pxout = m_cos_theta * px
                  + m_sin_theta * py;

            y = -m_sin_theta * xout
               + m_cos_theta * yout;
            pyout = -m_sin_theta * px
                   + m_cos_theta * py;

            t = tout;
            // ptout = pt;

            // assign updated values
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
            amrex::ParticleReal const step = slice_ds / std::sqrt(powi<2>(pt) - 1.0_prt);

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

            // initialize the three components of the axis-angle vector
            T_Real lambdax = 0_prt;
            T_Real lambday = 0_prt;
            T_Real lambdaz = 0_prt;

            // axis-angle vector components generating the remainder spin map
            lambdax = m_rel_const * ((2_prt*px + m_ks*y)*m_sin_G - (2_prt*py - m_ks*x)*m_1cos_G);
            lambday = m_rel_const * ((2_prt*py - m_ks*x)*m_sin_G + (2_prt*px + m_ks*y)*m_1cos_G);
            lambdaz = m_neg2_theta_1pG * pt * m_ibeta;

            // push the spin vector using the generator just determined
            rotate_spin(lambdax,lambday,lambdaz,sx,sy,sz);

            // axis-angle vector components generating the reference spin map
            lambdax = 0.0_prt;
            lambday = 0.0_prt;
            lambdaz = m_neg2_theta_1pG;

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

            // initialize linear map matrix elements
            Map6x6 R = Map6x6::Identity();
            Map6x6 Rrotation = Map6x6::Identity();
            Map6x6 Rfocusing = Map6x6::Identity();

            // access reference particle values to find beta*gamma^2
            amrex::ParticleReal const betgam2 = powi<2>(refpart.pt) - 1.0_prt;

            // compute phase advance per unit length (in rad/m) and
            // rotation angle (in rad)
            amrex::ParticleReal const alpha = m_ks * 0.5_prt;
            amrex::ParticleReal const theta = alpha*slice_ds;
            auto const [sin_theta, cos_theta] = amrex::Math::sincos(theta);

            amrex::ParticleReal const sin_theta_ialpha = sin_theta / alpha;
            amrex::ParticleReal const neg_alpha_sin_theta = -alpha * sin_theta;
            amrex::ParticleReal const ds_bg2 = slice_ds / betgam2;

            Rfocusing(1,1) = cos_theta;
            Rfocusing(1,2) = sin_theta_ialpha;
            Rfocusing(2,1) = neg_alpha_sin_theta;
            Rfocusing(2,2) = cos_theta;
            Rfocusing(3,3) = cos_theta;
            Rfocusing(3,4) = sin_theta_ialpha;
            Rfocusing(4,3) = neg_alpha_sin_theta;
            Rfocusing(4,4) = cos_theta;
            Rfocusing(5,6) = ds_bg2;

            Rrotation(1,1) = cos_theta;
            Rrotation(1,3) = sin_theta;
            Rrotation(2,2) = cos_theta;
            Rrotation(2,4) = sin_theta;
            Rrotation(3,1) = -sin_theta;
            Rrotation(3,3) = cos_theta;
            Rrotation(4,2) = -sin_theta;
            Rrotation(4,4) = cos_theta;

            R = Rrotation * Rfocusing;

            // apply the transverse rotation (roll) alignment error
            return rotate_aligned_map(R);
        }

        amrex::ParticleReal m_ks; //! solenoid strength in 1/m

      private:
        // constants that are independent of the individually tracked particle,
        // see: compute_constants() to refresh
        amrex::ParticleReal m_ds_bg2, m_cos_theta, m_sin_theta, m_sin_theta_ialpha, m_neg_alpha_sin_theta;
        amrex::ParticleReal m_rel_const, m_sin_G, m_1cos_G, m_neg2_theta_1pG, m_ibeta;
    };

    /** Sol: GPU state plus common host-only metadata */
    using Sol = Element<SolGPU>;

} // namespace impactx

IMPACTX_PUSH_EXTERN_TEMPLATE(impactx::elements::SolGPU)

#endif // IMPACTX_SOL_H
