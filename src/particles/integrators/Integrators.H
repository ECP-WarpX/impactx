/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
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

#include "particles/ReferenceParticle.H"

#include <AMReX_Extension.H>     // for AMREX_RESTRICT
#include <AMReX_REAL.H>          // for ParticleReal
#include <AMReX_SmallMatrix.H>   // for SmallVector


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
        amrex::ParticleReal const dz = (zout - zin) / static_cast<amrex::ParticleReal>(nsteps);
        amrex::ParticleReal const tau1 = dz * 0.5_prt;
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
        amrex::ParticleReal const dz = (zout - zin) / static_cast<amrex::ParticleReal>(nsteps);
        amrex::ParticleReal const tau1 = dz * 0.5_prt;
        amrex::ParticleReal const tau2 = dz * 0.5_prt;
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
        amrex::ParticleReal const dz = (zout - zin) / static_cast<amrex::ParticleReal>(nsteps);
        constexpr amrex::ParticleReal alpha = -0.2599210498948732_prt;  // C++26: 1_prt - std::pow(2_prt,1_prt/3_prt);
        amrex::ParticleReal const tau2 = dz / (1_prt + alpha);
        amrex::ParticleReal const tau1 = tau2*0.5_prt;
        amrex::ParticleReal const tau3 = alpha*tau1;
        amrex::ParticleReal const tau4 = (alpha - 1_prt)*tau2;

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

   /** A second-order symplectic integrator based on a Hamiltonian
    *  splitting H = H_1 + H_2.  This is a generalization of the second-
    *  order leapfrog algorithm.  For a detailed overview:
    *
    *  E. Hairer et al, Geometric Numerical Integration:  Structure-
    *  Preserving Algorithms for Ordinary Differential Equations,
    *  2nd ed, Springer, Berlin, 2006.
    *
    * @param particle  Particle phase space data in a 6-vector
    * @param zin  Initial value of independent variable (z-location)
    * @param zout  Final value of independent variable (z-location)
    * @param nsteps  Number of integration steps
    * @param refpart  Reference particle data
    * @param element  Element defining the two maps associated with H_1 and H_2
    */
    template <typename T_Real, typename T_Element>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void symp2_integrate_particle (
        amrex::SmallVector<T_Real, 6, 1> & particle,
        amrex::ParticleReal const zin,
        amrex::ParticleReal const zout,
        int const nsteps,
        RefPart const & refpart,
        T_Element const & element
    )
    {
        using namespace amrex::literals; // for _rt and _prt

        // initialize numerical integration parameters
        amrex::ParticleReal const dz = (zout - zin) / static_cast<amrex::ParticleReal>(nsteps);
        amrex::ParticleReal const tau1 = dz * 0.5_prt;
        amrex::ParticleReal const tau2 = dz;

        // initialize the value of the independent variable
        amrex::ParticleReal zeval = zin;

        // loop over integration steps
        for(int j=0; j < nsteps; ++j)
        {
            element.map1(tau1,particle,zeval,refpart);
            element.map2(tau2,particle,zeval,refpart);
            element.map1(tau1,particle,zeval,refpart);
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
    * @param particle  Particle phase space data in a 6-vector
    * @param zin  Initial value of independent variable (z-location)
    * @param zout  Final value of independent variable (z-location)
    * @param nsteps  Number of integration steps
    * @param refpart  Reference particle data
    * @param element  Element defining the maps associated with H_1, H_2, H_3
    */
    template <typename T_Real, typename T_Element>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void symp2_integrate_split3_particle (
        amrex::SmallVector<T_Real, 6, 1> & particle,
        amrex::ParticleReal const zin,
        amrex::ParticleReal const zout,
        int const nsteps,
        RefPart const & refpart,
        T_Element const & element
    )
    {
        using namespace amrex::literals; // for _rt and _prt

        // initialize numerical integration parameters
        amrex::ParticleReal const dz = (zout - zin) / static_cast<amrex::ParticleReal>(nsteps);
        amrex::ParticleReal const tau1 = dz * 0.5_prt;
        amrex::ParticleReal const tau2 = dz * 0.5_prt;
        amrex::ParticleReal const tau3 = dz;

        // initialize the value of the independent variable
        amrex::ParticleReal zeval = zin;

        // loop over integration steps
        for(int j=0; j < nsteps; ++j)
        {
            element.map1(tau1,particle,zeval,refpart);
            element.map2(tau2,particle,zeval,refpart);
            element.map3(tau3,particle,zeval,refpart);
            element.map2(tau2,particle,zeval,refpart);
            element.map1(tau1,particle,zeval,refpart);
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
    * @param particle  Particle phase space data in a 6-vector
    * @param zin  Initial value of independent variable (z-location)
    * @param zout  Final value of independent variable (z-location)
    * @param nsteps  Number of integration steps
    * @param refpart  Reference particle data
    * @param element  Element defining the two maps associated with H_1 and H_2
    */
    template <typename T_Real, typename T_Element>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void symp4_integrate_particle (
        amrex::SmallVector<T_Real, 6, 1> & particle,
        amrex::ParticleReal const zin,
        amrex::ParticleReal const zout,
        int const nsteps,
        RefPart const & refpart,
        T_Element const & element
    )
    {
        using namespace amrex::literals; // for _rt and _prt

        // initialize numerical integration parameters
        amrex::ParticleReal const dz = (zout - zin) / static_cast<amrex::ParticleReal>(nsteps);
        constexpr amrex::ParticleReal alpha = -0.2599210498948732_prt;  // C++26: 1_prt - std::pow(2_prt,1_prt/3_prt);
        amrex::ParticleReal const tau2 = dz/(1_prt + alpha);
        amrex::ParticleReal const tau1 = tau2*0.5_prt;
        amrex::ParticleReal const tau3 = alpha*tau1;
        amrex::ParticleReal const tau4 = (alpha - 1_prt)*tau2;

        // initialize the value of the independent variable
        amrex::ParticleReal zeval = zin;

        // loop over integration steps
        for (int j=0; j < nsteps; ++j)
        {
            element.map1(tau1,particle,zeval,refpart);
            element.map2(tau2,particle,zeval,refpart);
            element.map1(tau3,particle,zeval,refpart);
            element.map2(tau4,particle,zeval,refpart);
            element.map1(tau3,particle,zeval,refpart);
            element.map2(tau2,particle,zeval,refpart);
            element.map1(tau1,particle,zeval,refpart);
        }
    }


   /** A sixth-order symplectic integrator based on a Hamiltonian
    *  splitting H = H_1 + H_2.  The implementation below
    *  requires a reduced number of map evaluations (15 total)
    *  compared to direct use of the Yoshida construction, which
    *  requires three compositions of a fourth-order integrator step,
    *  each involving seven map evaluations (when expressed
    *  in the form obtained by Forest and Ruth, above).
    *  It is based on Section 5 of:
    *
    *  H. Yoshida, Phys. Lett. A 150, 292-268 (1990),
    *
    *  and appears also in Section V.3.2 (3.11) of:
    *
    *  E. Hairer et al, Geometric Numerical Integration, 2nd ed.,
    *  Springer-Verlag, Berlin, 2006.
    *
    * @param particle  Particle phase space data in a 6-vector
    * @param zin  Initial value of independent variable (z-location)
    * @param zout  Final value of independent variable (z-location)
    * @param nsteps  Number of integration steps
    * @param refpart  Reference particle data
    * @param element  Element defining the two maps associated with H_1 and H_2
    */
    template <typename T_Real, typename T_Element>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void symp6_integrate_particle (
        amrex::SmallVector<T_Real, 6, 1> & particle,
        amrex::ParticleReal const zin,
        amrex::ParticleReal const zout,
        int const nsteps,
        RefPart const & refpart,
        T_Element const & element
    )
    {
        using namespace amrex::literals; // for _rt and _prt

        // initialize numerical integration parameters
        amrex::ParticleReal const dz = (zout - zin) / static_cast<amrex::ParticleReal>(nsteps);
        amrex::ParticleReal const w0 = 1.31518632068391121888424973_prt;
        amrex::ParticleReal const w1 = -1.17767998417887100694641568_prt;
        amrex::ParticleReal const w2 = 0.23557321335935813368479318_prt;
        amrex::ParticleReal const w3 = 0.78451361047755726381949763_prt;

        amrex::ParticleReal const tau1 = w3 * dz / 2_prt;
        amrex::ParticleReal const tau2 = w3 * dz;
        amrex::ParticleReal const tau3 = (w3 + w2) * dz / 2_prt;
        amrex::ParticleReal const tau4 = w2 * dz;
        amrex::ParticleReal const tau5 = (w2 + w1) * dz / 2_prt;
        amrex::ParticleReal const tau6 = w1 * dz;
        amrex::ParticleReal const tau7 = (w1 + w0) * dz / 2_prt;
        amrex::ParticleReal const tau8 = w0 * dz;

        // initialize the value of the independent variable
        amrex::ParticleReal zeval = zin;

        // loop over integration steps
        for (int j=0; j < nsteps; ++j)
        {
            element.map1(tau1,particle,zeval,refpart);
            element.map2(tau2,particle,zeval,refpart);
            element.map1(tau3,particle,zeval,refpart);
            element.map2(tau4,particle,zeval,refpart);
            element.map1(tau5,particle,zeval,refpart);
            element.map2(tau6,particle,zeval,refpart);
            element.map1(tau7,particle,zeval,refpart);
            element.map2(tau8,particle,zeval,refpart);
            element.map1(tau7,particle,zeval,refpart);
            element.map2(tau6,particle,zeval,refpart);
            element.map1(tau5,particle,zeval,refpart);
            element.map2(tau4,particle,zeval,refpart);
            element.map1(tau3,particle,zeval,refpart);
            element.map2(tau2,particle,zeval,refpart);
            element.map1(tau1,particle,zeval,refpart);
        }
    }

   /** A second-order symplectic integration step based on a Hamiltonian
    *  splitting H = H_1 + H_2, with spin included.  The effect of spin
    *  is included by adding a third term in the splitting of the vector field,
    *  corresponding to the update from the Thomas-BMT equation.  This results
    *  in a net 3-term splitting, where map3 applies the spin map.
    *
    *  The algorithm is correct to second-order in the stepsize, with respect
    *  to both orbit and spin data.
    *
    * @param particle  Particle phase space data in a 6-vector
    * @param particle_spin  Particle spin data in a 3-vector
    * @param tau  Size of the step in s
    * @param zeval  Evaluation point in s
    * @param refpart  Reference particle data
    * @param element  Element defining the maps associated with H_1, H_2, and spin.
    */
    template <typename T_Real, typename T_Element>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void symp2_integrate_particle_spin_step (
        amrex::SmallVector<T_Real, 6, 1> & particle,
        amrex::SmallVector<T_Real, 3, 1> & particle_spin,
        amrex::ParticleReal const tau,
        amrex::ParticleReal & zeval,
        RefPart const & refpart,
        T_Element const & element
    )
    {
        using namespace amrex::literals; // for _rt and _prt

        // initialize numerical integration parameters
        amrex::ParticleReal const tau1 = tau*0.5_prt;
        amrex::ParticleReal const tau2 = tau*0.5_prt;
        amrex::ParticleReal const tau3 = tau;

        // a single second-order step
        element.map1(tau1,particle,zeval,refpart);
        element.map2(tau2,particle,zeval,refpart);
        element.map3(tau3,particle,particle_spin,zeval,refpart);
        element.map2(tau2,particle,zeval,refpart);
        element.map1(tau1,particle,zeval,refpart);
    }

   /** A fourth-order symplectic integration step based on a Hamiltonian
    *  splitting H = H_1 + H_2, with spin included.  The effect of spin
    *  is included by adding a third term in the splitting of the vector field,
    *  corresponding to the update from the Thomas-BMT equation.  This results
    *  in a net 3-term splitting, where map3 applies the spin map.
    *
    *  The algorithm is correct to fourth-order in the stepsize, with respect
    *  to both orbit and spin data.  Currently, this is obtained by using the
    *  Yoshida method.
    *
    * @param particle  Particle phase space data in a 6-vector
    * @param particle_spin  Particle spin data in a 3-vector
    * @param tau  Size of the step in s
    * @param zeval  Evaluation point in s
    * @param refpart  Reference particle data
    * @param element  Element defining the maps associated with H_1, H_2, and spin.
    */
    template <typename T_Real, typename T_Element>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void symp4_integrate_particle_spin_step (
        amrex::SmallVector<T_Real, 6, 1> & particle,
        amrex::SmallVector<T_Real, 3, 1> & particle_spin,
        amrex::ParticleReal const tau,
        amrex::ParticleReal & zeval,
        RefPart const & refpart,
        T_Element const & element
    )
    {
        using namespace amrex::literals; // for _rt and _prt

        // initialize numerical integration parameters (Yoshida)
        constexpr amrex::ParticleReal a = 1.2599210498948732_prt;  // C++26 std::pow(2_prt,1_prt/3_prt);
        constexpr amrex::ParticleReal z0 = 1_prt / (2_prt - a);
        constexpr amrex::ParticleReal z1 = -a / (2_prt - a);
        amrex::ParticleReal const tau1 = tau*z0;
        amrex::ParticleReal const tau2 = tau*z1;

        // a single fourth-order step
        symp2_integrate_particle_spin_step(particle,particle_spin,tau1,zeval,refpart,element);
        symp2_integrate_particle_spin_step(particle,particle_spin,tau2,zeval,refpart,element);
        symp2_integrate_particle_spin_step(particle,particle_spin,tau1,zeval,refpart,element);
    }


   /** A sixth-order symplectic integration step based on a Hamiltonian
    *  splitting H = H_1 + H_2, with spin included.  The effect of spin
    *  is included by adding a third term in the splitting of the vector field,
    *  corresponding to the update from the Thomas-BMT equation.  This results
    *  in a net 3-term splitting, where map3 applies the spin map.
    *
    *  The algorithm is correct to sixth-order in the stepsize, with respect
    *  to both orbit and spin data.  Currently, this is obtained by using the
    *  Yoshida method.
    *
    * @param particle  Particle phase space data in a 6-vector
    * @param particle_spin  Particle spin data in a 3-vector
    * @param tau  Size of the step in s
    * @param zeval  Evaluation point in s
    * @param refpart  Reference particle data
    * @param element  Element defining the maps associated with H_1, H_2, and spin.
    */
    template <typename T_Real, typename T_Element>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void symp6_integrate_particle_spin_step (
        amrex::SmallVector<T_Real, 6, 1> & particle,
        amrex::SmallVector<T_Real, 3, 1> & particle_spin,
        amrex::ParticleReal const tau,
        amrex::ParticleReal & zeval,
        RefPart const & refpart,
        T_Element const & element
    )
    {
        using namespace amrex::literals; // for _rt and _prt

        // initialize numerical integration parameters (Yoshida)
        constexpr amrex::ParticleReal a = 1.148698354997035_prt;  // C++26 std::pow(2_prt,1_prt/5_prt);
        constexpr amrex::ParticleReal z0 = 1_prt / (2_prt - a);
        constexpr amrex::ParticleReal z1 = -a / (2_prt - a);
        amrex::ParticleReal const tau1 = tau * z0;
        amrex::ParticleReal const tau2 = tau * z1;

        // a single sixth-order step
        symp4_integrate_particle_spin_step(particle,particle_spin,tau1,zeval,refpart,element);
        symp4_integrate_particle_spin_step(particle,particle_spin,tau2,zeval,refpart,element);
        symp4_integrate_particle_spin_step(particle,particle_spin,tau1,zeval,refpart,element);
    }


   /** A symplectic integrator based on a Hamiltonian
    *  splitting H = H_1 + H_2, with spin included.  The effect of spin
    *  is included by adding a third term in the splitting of the vector field,
    *  corresponding to the update from the Thomas-BMT equation.  This results
    *  in a net 3-term splitting, where map3 applies the spin map.
    *
    *  The algorithm is correct to second-order in the stepsize, with respect
    *  to both orbit and spin data.
    *
    * @param particle  Particle phase space data in a 6-vector
    * @param particle_spin  Particle spin data in a 3-vector
    * @param zin  Initial value of independent variable (z-location)
    * @param zout  Final value of independent variable (z-location)
    * @param nsteps  Number of integration steps
    * @param int_order Order of numerical integration (2, 4, or 6)
    * @param refpart  Reference particle data
    * @param element  Element defining the maps associated with H_1, H_2, and spin.
    */
    template <typename T_Real, typename T_Element>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void symp_integrate_particle_spin (
        amrex::SmallVector<T_Real, 6, 1> & particle,
        amrex::SmallVector<T_Real, 3, 1> & particle_spin,
        amrex::ParticleReal const zin,
        amrex::ParticleReal const zout,
        int const nsteps,
        int int_order,
        RefPart const & refpart,
        T_Element const & element
    )
    {
        using namespace amrex::literals; // for _rt and _prt

        // initialize numerical integration parameters
        amrex::ParticleReal const dz = (zout - zin) / static_cast<amrex::ParticleReal>(nsteps);

        // initialize the value of the independent variable
        amrex::ParticleReal zeval = zin;

        // loop over integration steps
        for(int j=0; j < nsteps; ++j)
        {
            if (int_order == 2) {
                symp2_integrate_particle_spin_step(particle,particle_spin,dz,zeval,refpart,element);
            } else if (int_order == 4) {
                symp4_integrate_particle_spin_step(particle,particle_spin,dz,zeval,refpart,element);
            } else if (int_order == 6) {
                symp6_integrate_particle_spin_step(particle,particle_spin,dz,zeval,refpart,element);
            } else {
                symp2_integrate_particle_spin_step(particle,particle_spin,dz,zeval,refpart,element);
            }
        }
    }


} // namespace impactx::integrators

#endif // IMPACTX_INTEGRATORS_H_
