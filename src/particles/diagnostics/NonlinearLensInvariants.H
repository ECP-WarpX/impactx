/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_INVARIANTS_H
#define IMPACTX_INVARIANTS_H

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuComplex.H>

#include <cmath>


namespace impactx::diagnostics
{
    /** Compute invariants within the nonlinear magnetic insert element
     *
     *  Invariants associated with a single short segment of the
     *  nonlinear magnetic insert described by V. Danilov and
     *  S. Nagaitsev, PRSTAB 13, 084002 (2010), Sect. V.A.
     */
    struct NonlinearLensInvariants
    {
        using PType = ImpactXParticleContainer::ParticleType;

        /** Two independent phase space functions in involution
         *
         * Independent functions that are tracked as invariants of motion
         * in the IOTA nonlinear magnetic insert.
         */
        struct Data {
            amrex::ParticleReal H; ///< first phase space function (Hamiltonian)
            amrex::ParticleReal I; ///< second phase space function ("second invariant")
        };

        /** Initialize the parameters for the invariants based on the beam
         *  distribution and parameters of the nonlinear insert element
         *
         * @param alpha - Twiss alpha for the bare lattice
         * @param beta - Twiss beta for the bare lattice (m)
         * @param tn - dimensionless strength of the nonlinear insert
         * @param cn - scale parameter of the nonlinear insert (m^[1/2])
         *
         */
        NonlinearLensInvariants (
            amrex::ParticleReal const alpha,
            amrex::ParticleReal const beta,
            amrex::ParticleReal const tn,
            amrex::ParticleReal const cn )
        : m_alpha(alpha), m_beta(beta), m_tn(tn), m_cn(cn)
        {
        }

        /** This is a nonlinear invariant functor, so that a variable of this type
         *  can be used like a nonlinear invariant function.
         *
         * @param x particle position in x
         * @param y particle position in y
         * @param px particle momentum in x
         * @param py particle momentum in y
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        Data operator() (
            amrex::ParticleReal const x,
            amrex::ParticleReal const y,
            amrex::ParticleReal const px,
            amrex::ParticleReal const py
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // a complex type with two amrex::ParticleReal
            using Complex = amrex::GpuComplex<amrex::ParticleReal>;

            // convert transverse phase space coordinates to normalized units
            amrex::ParticleReal const xn = x/(m_cn*std::sqrt(m_beta));
            amrex::ParticleReal const yn = y/(m_cn*std::sqrt(m_beta));
            amrex::ParticleReal const pxn = px*std::sqrt(m_beta)/m_cn + m_alpha*xn;
            amrex::ParticleReal const pyn = py*std::sqrt(m_beta)/m_cn + m_alpha*yn;

            // assign complex position zeta = x + iy
            Complex const zeta(xn, yn);
            Complex const zetaconj(xn, -yn);
            Complex const re1(1.0_prt, 0.0_prt);
            Complex const im1(0.0_prt, 1.0_prt);

            // compute croot = sqrt(1-zeta**2)
            Complex croot = amrex::pow(zeta, 2);
            croot = re1 - croot;
            croot = amrex::sqrt(croot);

            // compute carcsin = arcsin(zeta)
            Complex carcsin = im1*zeta + croot;
            carcsin = -im1*amrex::log(carcsin);

            // compute complex potential appearing in H invariant
            Complex Hpotential = zeta/croot;
            Hpotential = Hpotential*carcsin;

            // compute complex potential appearing in I invariant
            Complex Ipotential = (zeta+zetaconj)/croot;
            Ipotential = Ipotential*carcsin;

            // evaluate real parts
            amrex::ParticleReal Hinv = Hpotential.m_real;
            amrex::ParticleReal Iinv = Ipotential.m_real;

            // compute invariants H and I
            amrex::ParticleReal const Jz = xn*pyn - yn*pxn;
            Hinv = (std::pow(xn,2) + std::pow(yn,2) + std::pow(pxn,2) + std::pow(pyn,2))/2
                 + m_tn*Hinv;
            Iinv = std::pow(Jz,2) + std::pow(pxn,2) + std::pow(xn,2) + m_tn*Iinv;

            return {Hinv, Iinv};

        }

    private:
        amrex::ParticleReal m_alpha; //! Twiss alpha
        amrex::ParticleReal m_beta; //! Twiss beta (m)
        amrex::ParticleReal m_tn; //! dimensionless strength of the nonlinear insert
        amrex::ParticleReal m_cn; //! scale parameter of the nonlinear insert (m^[1/2])
    };

} // namespace impactx

#endif // IMPACTX_INVARIANTS_H
