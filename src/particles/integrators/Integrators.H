/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_INTEGRATORS_H_
#define IMPACTX_INTEGRATORS_H_

#include <AMReX_Extension.H>  // for AMREX_RESTRICT
#include <AMReX_REAL.H>       // for ParticleReal


namespace impactx::integrators
{

   /** A second-order symplectic integrator based on a Hamiltonian
    *  splitting H = H_1 + H_2.  This is a generalization of the second-
    *  order leapfrog algorithm.  For a detailed overview:
    *
    *  E. Hairer et al, Geometric Numerical Integration:  Structure-
    *  Preserving Algorithms for Ordinary Differential Equations,
    *  2nd ed, Springer, Berlin, 2006.
    *
    * @param refpart  Reference particle data
    * @param zin  Initial value of independent variable (z-location)
    * @param zout  Final value of independent variable (z-location)
    * @param nsteps  Number of integration steps
    * @param element  Element defining the two maps associated with H_1 and H_2
    */
    template <typename T_Element>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void symp2_integrate (
        RefPart & refpart,
        amrex::ParticleReal const zin,
        amrex::ParticleReal const zout,
        int const nsteps,
        T_Element const & element
    )
    {
        using namespace amrex::literals; // for _rt and _prt

        // initialize numerical integration parameters
        amrex::ParticleReal const dz = (zout-zin)/nsteps;
        amrex::ParticleReal const tau1 = dz/2.0_prt;
        amrex::ParticleReal const tau2 = dz;

        // initialize the value of the independent variable
        amrex::ParticleReal zeval = zin;

        // loop over integration steps
        for(int j=0; j < nsteps; ++j)
        {
            element.map1(tau1,refpart,zeval);
            element.map2(tau2,refpart,zeval);
            element.map1(tau1,refpart,zeval);
        }
    }

   /** A second-order symplectic integrator based on a Hamiltonian
    *  splitting H = H_1 + H_2 + H_3.  This is the result of applying
    *  the two-terms splitting of symp2_integrate twice in succession.
    *  The form shown here appears, for example, in:
    *
    *  R. D. Ryne, Computational Methods in Accelerator Physics,
    *  USPAS lecture notes, 2009.
    *
    * @param refpart  Reference particle data
    * @param zin  Initial value of independent variable (z-location)
    * @param zout  Final value of independent variable (z-location)
    * @param nsteps  Number of integration steps
    * @param element  Element defining the two maps associated with H_1 and H_2
    */
    template <typename T_Element>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void symp2_integrate_split3 (
        RefPart & refpart,
        amrex::ParticleReal const zin,
        amrex::ParticleReal const zout,
        int const nsteps,
        T_Element const & element
    )
    {
        using namespace amrex::literals; // for _rt and _prt

        // initialize numerical integration parameters
        amrex::ParticleReal const dz = (zout-zin)/nsteps;
        amrex::ParticleReal const tau1 = dz/2.0_prt;
        amrex::ParticleReal const tau2 = dz/2.0_prt;
        amrex::ParticleReal const tau3 = dz;

        // initialize the value of the independent variable
        amrex::ParticleReal zeval = zin;

        // loop over integration steps
        for(int j=0; j < nsteps; ++j)
        {
            element.map1(tau1,refpart,zeval);
            element.map2(tau2,refpart,zeval);
            element.map3(tau3,refpart,zeval);
            element.map2(tau2,refpart,zeval);
            element.map1(tau1,refpart,zeval);
        }
    }

   /** A fourth-order symplectic integrator based on a Hamiltonian
    *  splitting H = H_1 + H_2.  This is the result of applying
    *  the two-terms splitting of symp2_integrate together with
    *  the method of Yoshida.
    *
    *  H. Yoshida, Phys. Lett. A 150, 292-268 (1990).
    *
    *  The form shown here appears, for example, in:
    *
    *  E. Forest and R. D. Ruth, Physica D 43, 105-117 (1990).
    *
    * @param refpart  Reference particle data
    * @param zin  Initial value of independent variable (z-location)
    * @param zout  Final value of independent variable (z-location)
    * @param nsteps  Number of integration steps
    * @param element  Element defining the two maps associated with H_1 and H_2
    */
    template <typename T_Element>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void symp4_integrate (
        RefPart & refpart,
        amrex::ParticleReal const zin,
        amrex::ParticleReal const zout,
        int const nsteps,
        T_Element const & element
    )
    {
        using namespace amrex::literals; // for _rt and _prt

        // initialize numerical integration parameters
        amrex::ParticleReal const dz = (zout-zin)/nsteps;
        amrex::ParticleReal const alpha = 1.0_prt - pow(2.0_prt,1.0/3.0);
        amrex::ParticleReal const tau2 = dz/(1.0_prt + alpha);
        amrex::ParticleReal const tau1 = tau2/2.0_prt;
        amrex::ParticleReal const tau3 = alpha*tau1;
        amrex::ParticleReal const tau4 = (alpha - 1.0_prt)*tau2;

        // initialize the value of the independent variable
        amrex::ParticleReal zeval = zin;

        // loop over integration steps
        for (int j=0; j < nsteps; ++j)
        {
            element.map1(tau1,refpart,zeval);
            element.map2(tau2,refpart,zeval);
            element.map1(tau3,refpart,zeval);
            element.map2(tau4,refpart,zeval);
            element.map1(tau3,refpart,zeval);
            element.map2(tau2,refpart,zeval);
            element.map1(tau1,refpart,zeval);
        }
    }
} // namespace impactx::integrators

#endif // IMPACTX_INTEGRATORS_H_
