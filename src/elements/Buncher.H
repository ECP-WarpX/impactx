/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_BUNCHER_H
#define IMPACTX_BUNCHER_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/gpuelement.H"
#include "mixin/alignment.H"
#include "mixin/beamoptic.H"
#include "mixin/thin.H"
#include "mixin/lineartransport.H"
#include "mixin/spintransport.H"
#include "mixin/nofinalize.H"
#include "mixin/pipeaperture.H"

#include <AMReX_Extension.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>
#include <AMReX_SIMD.H>

#include <cmath>
#include <stdexcept>

namespace impactx::elements
{
    struct BuncherGPU
    : public mixin::BeamOptic<BuncherGPU>,
      public mixin::LinearTransport<BuncherGPU>,
      public mixin::Thin,
      public mixin::Alignment,
      public mixin::PipeAperture,
      public mixin::SpinTransport,
      public mixin::NoFinalize
      // At least on Intel AVX512, there is a small overhead to vectorize this element, see
      // https://github.com/BLAST-ImpactX/impactx/pull/1002
      // public amrex::simd::Vectorized<amrex::simd::native_simd_size_particlereal>
    {
        static constexpr auto type = "Buncher";
        using PType = ImpactXParticleContainer::ParticleType;

        /** A short RF cavity element at zero crossing for bunching
         *
         * @param V Normalized RF voltage drop V = Emax*L/(c*Brho)
         * @param k Wavenumber of RF in 1/m
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param aperture_x horizontal half-aperture in m
         * @param aperture_y vertical half-aperture in m
         */
        BuncherGPU (
            amrex::ParticleReal V,
            amrex::ParticleReal k,
            amrex::ParticleReal dx = DEFAULT_dx,
            amrex::ParticleReal dy = DEFAULT_dy,
            amrex::ParticleReal rotation_degree = DEFAULT_rotation_degree,
            amrex::ParticleReal aperture_x = DEFAULT_aperture_x,
            amrex::ParticleReal aperture_y = DEFAULT_aperture_y
        )
        : Alignment(dx, dy, rotation_degree),
          PipeAperture(aperture_x, aperture_y),
          m_V(V), m_k(k)
        {
        }

        /** Reverse the element (negate V) */
        void reverse () { m_V = -m_V; }

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

            // find beta*gamma^2
            amrex::ParticleReal const betgam2 = powi<2>(refpart.pt) - 1.0_prt;

            m_neg_kV = -m_k * m_V;
            m_kV_r2bg2 = -m_neg_kV / (2.0_prt * betgam2);

            m_beta = refpart.beta();
        }

        /** This is a buncher functor, so that a variable of this type can be used like a
         *  buncher function.
         *
         * The @see compute_constants method must be called before pushing particles through this operator.
         *
         * @param x particle position in x
         * @param y particle position in y
         * @param t particle position in t
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t
         * @param idcpu particle global index (unused)
         * @param refpart reference particle (unused)
         */
        template<typename T_Real=amrex::ParticleReal, typename T_IdCpu=uint64_t>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
            T_Real const & AMREX_RESTRICT x,
            T_Real const & AMREX_RESTRICT y,
            T_Real const & AMREX_RESTRICT t,
            T_Real & AMREX_RESTRICT px,
            T_Real & AMREX_RESTRICT py,
            T_Real & AMREX_RESTRICT pt,
            [[maybe_unused]] T_IdCpu const & AMREX_RESTRICT idcpu,
            [[maybe_unused]] RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // initialize output values of momenta
            T_Real pxout = px;
            T_Real pyout = py;
            T_Real ptout = pt;

            // advance position and momentum
            // xout = x;
            pxout = px + m_kV_r2bg2 * x;
            // yout = y;
            pyout = py + m_kV_r2bg2 * y;
            // tout = t;
            ptout = pt + m_neg_kV   * t;

            // assign updated momenta
            px = pxout;
            py = pyout;
            pt = ptout;
        }

        /** This pushes the reference particle. */
        using Thin::operator();


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
            T_Real const & AMREX_RESTRICT x,
            T_Real const & AMREX_RESTRICT y,
            T_Real const & AMREX_RESTRICT t,
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
            using namespace std; // for cmath(float)

            // reference particle relativistic factors
            amrex::ParticleReal const gamma_ref = refpart.gamma();
            amrex::ParticleReal const gyro_anomaly = refpart.gyromagnetic_anomaly;

            // initialize the three components of the axis-angle vector
            T_Real lambdax = 0_prt;
            T_Real lambday = 0_prt;
            T_Real lambdaz = 0_prt;

            // Integrated magnetic field normalized by q/mc:
            T_Real const Bx = m_kV_r2bg2 * y * m_beta * gamma_ref;
            T_Real const By = -m_kV_r2bg2 * x * m_beta * gamma_ref;
            T_Real const Bz = 0_prt;

            // Integrated electric field normalized by q/mc^2:
            T_Real const Ex = 0.0_prt;
            T_Real const Ey = 0.0_prt;
            T_Real const Ez = -m_neg_kV * t * m_beta * gamma_ref;

            // The Thomas-BMT generator is evaluated at the momentum halfway
            // through the kick, so that reversing the element (negating V)
            // produces the exact inverse rotation. The buncher also changes the
            // energy deviation, so pt is averaged along with the transverse
            // momenta.
            T_Real const px_in = px;
            T_Real const py_in = py;
            T_Real const pt_in = pt;

            // phase space push
            (*this)(x, y, t, px, py, pt, idcpu, refpart);

            // momenta at the kick midpoint
            T_Real const px_mid = 0.5_prt * (px_in + px);
            T_Real const py_mid = 0.5_prt * (py_in + py);
            T_Real const pt_mid = 0.5_prt * (pt_in + pt);

            // Quantities required to evaluate the full Thomas-BMT precession
            // vector, at the midpoint momentum.
            T_Real const Pnorm = sqrt(1_prt - 2_prt * pt_mid / m_beta + pt_mid*pt_mid);
            T_Real const iPnorm = 1_prt/Pnorm;
            T_Real const gamma = gamma_ref * (1_prt - pt_mid*m_beta);
            T_Real const Ps = sqrt(Pnorm*Pnorm - px_mid*px_mid - py_mid*py_mid);
            T_Real const ux = px_mid * iPnorm;
            T_Real const uy = py_mid * iPnorm;
            T_Real const uz = Ps * iPnorm;

            // Zero curvature in a buncher
            amrex::ParticleReal const h = 0.0_prt;

            // Evaluation of the full Thomas-BMT precession vector.
            tbmt_precession_vector(x,ux,uy,uz,gamma,h,gyro_anomaly,Bx,By,Bz,Ex,Ey,Ez,lambdax,lambday,lambdaz);

            // push the spin vector using the generator just determined
            rotate_spin(lambdax, lambday, lambdaz, sx, sy, sz);
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

            // find beta*gamma^2
            amrex::ParticleReal const betgam2 = powi<2>(refpart.pt) - 1.0_prt;

            amrex::ParticleReal const neg_kV = -m_k * m_V;
            amrex::ParticleReal const kV_r2bg2 = -neg_kV / (2.0_prt * betgam2);

            // initialize linear map matrix elements
            Map6x6 R = Map6x6::Identity();

            // off-identity matrix elements
            R(2,1) = kV_r2bg2;
            R(4,3) = kV_r2bg2;
            R(6,5) = neg_kV;

            // apply the transverse rotation (roll) alignment error
            return rotate_aligned_map(R);
        }

        amrex::ParticleReal m_V; //! normalized (max) RF voltage drop.
        amrex::ParticleReal m_k; //! RF wavenumber in 1/m.

      private:
        // constants that are independent of the individually tracked particle,
        // see: compute_constants() to refresh
        amrex::ParticleReal m_neg_kV;    //! -m_k*m_V
        amrex::ParticleReal m_kV_r2bg2;  //! m_k*m_V/(2.0_prt*betgam2)
        amrex::ParticleReal m_beta;  // relativistic factor
    };

    /** Buncher: GPU state plus common host-only metadata */
    using Buncher = Element<BuncherGPU>;

} // namespace impactx

IMPACTX_PUSH_EXTERN_TEMPLATE(impactx::elements::BuncherGPU)

#endif // IMPACTX_BUNCHER_H
