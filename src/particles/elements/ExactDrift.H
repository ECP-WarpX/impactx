/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
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
#include "mixin/beamoptic.H"
#include "mixin/thick.H"
#include "mixin/named.H"
#include "mixin/nofinalize.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>

#include <cmath>


namespace impactx
{
    struct ExactDrift
    : public elements::Named,
      public elements::BeamOptic<ExactDrift>,
      public elements::Thick,
      public elements::Alignment,
      public elements::NoFinalize
    {
        static constexpr auto type = "ExactDrift";
        using PType = ImpactXParticleContainer::ParticleType;

        /** A drift using the exact nonlinear transfer map
         *
         * @param ds Segment length in m
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param nslice number of slices used for the application of space charge
         * @param name a user defined and not necessarily unique name of the element
         */
        ExactDrift (
            amrex::ParticleReal ds,
            amrex::ParticleReal dx = 0,
            amrex::ParticleReal dy = 0,
            amrex::ParticleReal rotation_degree = 0,
            int nslice = 1,
            std::optional<std::string> name = std::nullopt
        )
        : Named(name),
          Thick(ds, nslice),
          Alignment(dx, dy, rotation_degree)
        {
        }

        /** Push all particles */
        using BeamOptic::operator();

        /** This is an exactdrift functor, so that a variable of this type can be used like
         *  an exactdrift function.
         *
         * @param x particle position in x
         * @param y particle position in y
         * @param t particle position in t
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t
         * @param idcpu particle global index (unused)
         * @param refpart reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
                amrex::ParticleReal & AMREX_RESTRICT x,
                amrex::ParticleReal & AMREX_RESTRICT y,
                amrex::ParticleReal & AMREX_RESTRICT t,
                amrex::ParticleReal & AMREX_RESTRICT px,
                amrex::ParticleReal & AMREX_RESTRICT py,
                amrex::ParticleReal & AMREX_RESTRICT pt,
                [[maybe_unused]] uint64_t & AMREX_RESTRICT idcpu,
                RefPart const & refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // shift due to alignment errors of the element
            shift_in(x, y, px, py);

            // initialize output values
            amrex::ParticleReal const xout = x;
            amrex::ParticleReal const yout = y;
            amrex::ParticleReal const tout = t;

            // initialize output values of momenta
            amrex::ParticleReal const pxout = px;
            amrex::ParticleReal const pyout = py;
            amrex::ParticleReal const ptout = pt;

            // length of the current slice
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // access reference particle values to find beta, beta*gamma
            amrex::ParticleReal const bet = refpart.beta();
            amrex::ParticleReal const betgam = refpart.beta_gamma();

            // compute the radical in the denominator (= pz):
            amrex::ParticleReal const pzden = std::sqrt(std::pow(pt-1_prt/bet,2) -
                                1_prt / std::pow(betgam,2) - std::pow(px,2) - std::pow(py,2));

            // advance position and momentum (exact drift)
            x = xout + slice_ds * px / pzden;
            // pxout = px;
            y = yout + slice_ds * py / pzden;
            // pyout = py;
            t = tout - slice_ds * (1_prt / bet + (pt-1_prt/bet)/pzden);
            // ptout = pt;

            // assign updated momenta
            px = pxout;
            py = pyout;
            pt = ptout;

            // undo shift due to alignment errors of the element
            shift_out(x, y, px, py);
        }

        /** This pushes the reference particle.
         *
         * @param[in,out] refpart reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (RefPart & AMREX_RESTRICT refpart) const {

            using namespace amrex::literals; // for _rt and _prt

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
            amrex::ParticleReal const step = slice_ds /std::sqrt(std::pow(pt,2)-1.0_prt);

            // advance position and momentum (drift)
            refpart.x = x + step*px;
            refpart.y = y + step*py;
            refpart.z = z + step*pz;
            refpart.t = t - step*pt;

            // advance integrated path length
            refpart.s = s + slice_ds;

        }
    };

} // namespace impactx

#endif // IMPACTX_EXACTDRIFT_H
