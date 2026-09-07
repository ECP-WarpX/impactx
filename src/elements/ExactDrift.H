/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_EXACTDRIFT_H
#define IMPACTX_EXACTDRIFT_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/alignment.H"
#include "mixin/pipeaperture.H"
#include "mixin/beamoptic.H"
#include "mixin/thick.H"
#include "mixin/lineartransport.H"
#include "mixin/spintransport.H"
#include "mixin/named.H"
#include "mixin/nofinalize.H"

#include <AMReX_Extension.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>
#include <AMReX_SIMD.H>

#include <cmath>
#include <stdexcept>

namespace impactx::elements
{
    struct ExactDrift
    : public mixin::Named,
      public mixin::BeamOptic<ExactDrift>,
      public mixin::LinearTransport<ExactDrift>,
      public mixin::Thick,
      public mixin::Alignment,
      public mixin::PipeAperture,
      public mixin::SpinTransport,
      public mixin::NoFinalize,
      public amrex::simd::Vectorized<amrex::simd::native_simd_size_particlereal>
    {
        static constexpr auto type = "ExactDrift";
        using PType = ImpactXParticleContainer::ParticleType;

        /** A drift using the exact nonlinear transfer map
         *
         * @param ds Segment length in m
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param aperture_x horizontal half-aperture in m
         * @param aperture_y vertical half-aperture in m
         * @param nslice number of slices used for the application of space charge
         * @param name a user defined and not necessarily unique name of the element
         */
        ExactDrift (
            amrex::ParticleReal ds,
            amrex::ParticleReal dx = DEFAULT_dx,
            amrex::ParticleReal dy = DEFAULT_dy,
            amrex::ParticleReal rotation_degree = DEFAULT_rotation_degree,
            amrex::ParticleReal aperture_x = DEFAULT_aperture_x,
            amrex::ParticleReal aperture_y = DEFAULT_aperture_y,
            int nslice = DEFAULT_nslice,
            std::optional<std::string> name = DEFAULT_name
        )
        : Named(std::move(name)),
          Thick(ds, nslice),
          Alignment(dx, dy, rotation_degree),
          PipeAperture(aperture_x, aperture_y)
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

            // access reference particle values to find beta, 1/(beta*gamma)**2
            m_ibeta = 1_prt / refpart.beta();
            m_ibetagam2 = 1_prt / powi<2>(refpart.beta_gamma());
        }

        /** This is an exactdrift functor, so that a variable of this type can be used like
         *  an exactdrift function.
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
            T_Real const & AMREX_RESTRICT px,
            T_Real const & AMREX_RESTRICT py,
            T_Real const & AMREX_RESTRICT pt,
            T_IdCpu & AMREX_RESTRICT idcpu,
            [[maybe_unused]] RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using namespace std; // for cmath(float)
            using amrex::Math::powi;

            // initialize output values
            T_Real const xout = x;
            T_Real const yout = y;
            T_Real const tout = t;

            // initialize output values of momenta
            // T_Real const pxout = px;
            // T_Real const pyout = py;
            // T_Real const ptout = pt;

            // Longitudinal momentum squared, pz^2 = (pt - 1/beta)^2 - 1/(beta*gamma)^2 - px^2 - py^2.
            // A particle whose transverse momentum exceeds its total momentum has no real
            // longitudinal momentum (pz^2 <= 0): it lies outside the domain of the exact drift
            // map. Flag it invalid at the source so it is dropped at the next redistribute,
            // mirroring the exact bend (ExactSbend), instead of the otherwise unguarded
            // sqrt(<0) = NaN poisoning the whole-beam moments. The mask is a comparison, so it
            // survives -ffast-math; it is never set in double precision, leaving that path
            // numerically unchanged.
            T_Real const pzden2 =
                powi<2>(pt - m_ibeta) -
                m_ibetagam2 -
                powi<2>(px) -
                powi<2>(py);
            amrex::ParticleIDWrapper<T_IdCpu>{idcpu}.make_invalid(pzden2 <= 0_prt);
            T_Real const inv_pzden = 1_prt / sqrt(pzden2);

            // advance position and momentum (exact drift)
            x = xout + m_slice_ds * px * inv_pzden;
            // pxout = px;
            y = yout + m_slice_ds * py * inv_pzden;
            // pyout = py;
            t = tout - m_slice_ds * (m_ibeta + (pt - m_ibeta) * inv_pzden);
            // ptout = pt;

            // assign updated momenta
            // px = pxout;
            // py = pyout;
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
            amrex::ParticleReal const step = slice_ds /std::sqrt(powi<2>(pt)-1.0_prt);

            // advance position and momentum (drift)
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
         * @param sx particle spin x-component (unused)
         * @param sy particle spin y-component (unused)
         * @param sz particle spin z-component (unused)
         * @param idcpu particle global index
         * @param refpart reference particle (unused)
         */
        template<typename T_Real=amrex::ParticleReal, typename T_IdCpu=uint64_t>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void spin_and_phasespace_push (
            T_Real & AMREX_RESTRICT x,
            T_Real & AMREX_RESTRICT y,
            T_Real & AMREX_RESTRICT t,
            T_Real const & AMREX_RESTRICT px,
            T_Real const & AMREX_RESTRICT py,
            T_Real const & AMREX_RESTRICT pt,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT sx,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT sy,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT sz,
            T_IdCpu & AMREX_RESTRICT idcpu,
            RefPart const & AMREX_RESTRICT refpart
        ) const
        {
             // spin is unaffected

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

            // assign linear map matrix elements
            Map6x6 R = Map6x6::Identity();
            R(1,2) = slice_ds;
            R(3,4) = slice_ds;
            R(5,6) = slice_ds / betgam2;

            // apply the transverse rotation (roll) alignment error
            return rotate_aligned_map(R);
        }

      private:
        // constants that are independent of the individually tracked particle,
        // see: compute_constants() to refresh
        amrex::ParticleReal m_slice_ds;  //! m_ds / nslice();
        amrex::ParticleReal m_ibeta, m_ibetagam2;
    };

} // namespace impactx

IMPACTX_PUSH_EXTERN_TEMPLATE(impactx::elements::ExactDrift)

#endif // IMPACTX_EXACTDRIFT_H
