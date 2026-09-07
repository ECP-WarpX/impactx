/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_EXACTQUAD_H
#define IMPACTX_EXACTQUAD_H

#include "particles/ImpactXParticleContainer.H"
#include "particles/integrators/Integrators.H"
#include "mixin/alignment.H"
#include "mixin/pipeaperture.H"
#include "mixin/beamoptic.H"
#include "mixin/thick.H"
#include "mixin/named.H"
#include "mixin/nofinalize.H"
#include "mixin/lineartransport.H"
#include "mixin/spintransport.H"

#include <AMReX_Extension.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>
#include <AMReX_SIMD.H>

#include <cmath>
#include <stdexcept>

namespace impactx::elements
{
    struct ExactQuad
    : public mixin::Named,
      public mixin::BeamOptic<ExactQuad>,
      public mixin::LinearTransport<ExactQuad>,
      public mixin::Thick,
      public mixin::Alignment,
      public mixin::PipeAperture,
      public mixin::SpinTransport,
      public mixin::NoFinalize,
      public amrex::simd::Vectorized<amrex::simd::native_simd_size_particlereal>
    {
        static constexpr auto type = "ExactQuad";
        using PType = ImpactXParticleContainer::ParticleType;

        static constexpr int DEFAULT_unit = 0;
        static constexpr int DEFAULT_int_order = 2;
        static constexpr int DEFAULT_mapsteps = 10;

        /** A Quadrupole magnet treated using the exact Hamiltonian, which
         *  includes all nonlinear kinematic effects.  Like the element Quad,
         *  a hard-edge model of the magnetic field is used.  Symplectic
         *  integration is used based on a splitting of the Hamiltonian
         *  described in [to be cited].
         *
         * @param ds Segment length in m.
         * @param k  Quadrupole strength in m^(-2) (MADX convention)
         *           = (gradient in T/m) / (rigidity in T-m)
         *        OR Quadrupole strength in T/m (MaryLie convention)
         *           k > 0 horizontal focusing
         *           k < 0 horizontal defocusing
         * @param unit  Unit specification
         *           unit = 0 MADX convention (default)
         *           unit = 1 MaryLie convention
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param aperture_x horizontal half-aperture in m
         * @param aperture_y vertical half-aperture in m
         * @param int_order order used for symplectic integration (2, 4, or 6)
         * @param mapsteps number of integration steps per slice used for
         *        particle push in the applied fields
         * @param nslice number of slices used for the application of space charge
         * @param name a user defined and not necessarily unique name of the element
         */
        ExactQuad (
            amrex::ParticleReal ds,
            amrex::ParticleReal k,
            int unit = DEFAULT_unit,
            amrex::ParticleReal dx = DEFAULT_dx,
            amrex::ParticleReal dy = DEFAULT_dy,
            amrex::ParticleReal rotation_degree = DEFAULT_rotation_degree,
            amrex::ParticleReal aperture_x = DEFAULT_aperture_x,
            amrex::ParticleReal aperture_y = DEFAULT_aperture_y,
            int int_order = DEFAULT_int_order,
            int mapsteps = DEFAULT_mapsteps,
            int nslice = DEFAULT_nslice,
            std::optional<std::string> name = DEFAULT_name
        )
        : Named(std::move(name)),
          Thick(ds, nslice),
          Alignment(dx, dy, rotation_degree),
          PipeAperture(aperture_x, aperture_y),
          m_k(k),m_unit(unit),m_int_order(int_order),m_mapsteps(mapsteps)
        {
            if (m_int_order != 2 && m_int_order != 4 && m_int_order != 6)
                throw std::invalid_argument("ExactQuad: The order used for symplectic integration must be 2, 4 or 6.");

            if (m_mapsteps <= 0)
                throw std::invalid_argument("ExactQuad: mapsteps must be > 0!");
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

            // find beta*gamma^2
            amrex::ParticleReal const pt_ref = refpart.pt;
            m_betgam2 = powi<2>(pt_ref) - 1.0_prt;
            m_ibetgam2 = 1_prt / m_betgam2;
            m_beta = refpart.beta();
            m_ibeta = 1_prt / m_beta;

            // normalize quad units to MAD-X convention if needed
            m_g = m_unit == 1 ? m_k / refpart.rigidity_Tm() : m_k;

            // compute phase advance per unit length in s (in rad/m)
            m_omega = std::sqrt(std::abs(m_g));
        }

        /** Wrapper that carries the particle global index through the (id-agnostic)
         *  symplectic integrator, so that the exact-drift substep (map2) can invalidate a
         *  particle that leaves the map domain at the source. It forwards map1 unchanged and
         *  adds the domain guard around map2; the integrator drives it exactly like the
         *  element itself.
         */
        template<typename T_IdCpu>
        struct MapWithDomainGuard
        {
            ExactQuad const & AMREX_RESTRICT m_elem;
            T_IdCpu & AMREX_RESTRICT m_idcpu;

            template<typename T_Real>
            AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
            void map1 (amrex::ParticleReal const tau,
                       amrex::SmallVector<T_Real, 6, 1> & particle,
                       amrex::ParticleReal & zeval,
                       RefPart const & AMREX_RESTRICT refpart) const
            {
                m_elem.map1(tau, particle, zeval, refpart);
            }

            template<typename T_Real>
            AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
            void map2 (amrex::ParticleReal const tau,
                       amrex::SmallVector<T_Real, 6, 1> & particle,
                       amrex::ParticleReal & zeval,
                       RefPart const & AMREX_RESTRICT refpart) const
            {
                using namespace amrex::literals; // for _prt
                using amrex::Math::powi;

                // Source domain guard for the exact-drift substep. The longitudinal momentum
                // squared is pz^2 = (pt - 1/beta)^2 - 1/(beta*gamma)^2 - px^2 - py^2. A deep-tail
                // particle amplified by a strong quadrupole can reach a transverse momentum that
                // exceeds its total momentum (pz^2 <= 0) -- a genuinely lost particle, outside
                // the domain of the map. Flag it invalid BEFORE the sqrt so it is dropped at the
                // next redistribute, mirroring the exact bend (ExactSbend), instead of the
                // otherwise unguarded sqrt(<0) = NaN poisoning the whole-beam moments. The check
                // is a comparison, so it survives -ffast-math (a post-hoc NaN test would not);
                // it is never set in double precision, leaving that path numerically unchanged.
                T_Real const px = particle(2);
                T_Real const py = particle(4);
                T_Real const pt = particle(6);
                T_Real const pzden2 = powi<2>(pt - m_elem.m_ibeta) - m_elem.m_ibetgam2
                                      - powi<2>(px) - powi<2>(py);
                amrex::ParticleIDWrapper<T_IdCpu>{m_idcpu}.make_invalid(pzden2 <= 0_prt);

                m_elem.map2(tau, particle, zeval, refpart);
            }
        };

        /** This is a quad functor, so that a variable of this type can be used like a quad function.
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
            T_IdCpu & AMREX_RESTRICT idcpu,
            RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // numerical integration parameters
            amrex::ParticleReal const zin = 0_prt;
            amrex::ParticleReal const zout = m_slice_ds;
            int const nsteps = m_mapsteps;

            // initialize phase space 6-vector
            amrex::SmallVector<T_Real, 6, 1> particle{
                x, px, y, py, t, pt
            };

            // Drive the symplectic integrator through a wrapper that carries the particle id,
            // so the exact-drift substep (map2) can reject a particle that leaves the map domain.
            MapWithDomainGuard<T_IdCpu> const guarded{*this, idcpu};

            // call integrator to advance through slice (int_order = 2 or 4, otherwise default 2)
            if (m_int_order == 2) {
                integrators::symp2_integrate_particle(particle,zin,zout,nsteps,refpart,guarded);
            } else if (m_int_order == 4) {
                integrators::symp4_integrate_particle(particle,zin,zout,nsteps,refpart,guarded);
            } else if (m_int_order == 6) {
                integrators::symp6_integrate_particle(particle,zin,zout,nsteps,refpart,guarded);
            } else {
                integrators::symp2_integrate_particle(particle,zin,zout,nsteps,refpart,guarded);
            }

            // assign updated values
            x = particle(1);
            px = particle(2);
            y = particle(3);
            py = particle(4);
            t = particle(5);
            pt = particle(6);
        }

        /** This pushes a particle in an ExactQuad element through
         *  the symplectic map associated with H_1 in the Hamiltonian
         *  splitting H = H_1 + H_2.  This is exactly the linear
         *  map of a quadrupole.
         *
         * @param tau Map step size in m
         * @param[in,out] particle particle phase space 6-vector
         * @param[in,out] zeval Longitudinal on-axis location in m
         * @param[in] refpart reference particle
         */
        template<typename T_Real>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void map1 (amrex::ParticleReal const tau,
                   amrex::SmallVector<T_Real, 6, 1> & particle,
                   amrex::ParticleReal & zeval,
                   [[maybe_unused]] RefPart const & AMREX_RESTRICT refpart) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using namespace std; // for cmath(float)

            T_Real const x = particle(1);
            T_Real const px = particle(2);
            T_Real const y = particle(3);
            T_Real const py = particle(4);
            T_Real const t = particle(5);
            T_Real const pt = particle(6);

            T_Real xout = x;
            T_Real pxout = px;
            T_Real yout = y;
            T_Real pyout = py;
            T_Real tout = t;
            T_Real ptout = pt;

            amrex::ParticleReal const sin_omega_ds = sin(m_omega*tau);
            amrex::ParticleReal const cos_omega_ds = cos(m_omega*tau);
            amrex::ParticleReal const sinh_omega_ds = sinh(m_omega*tau);
            amrex::ParticleReal const cosh_omega_ds = cosh(m_omega*tau);
            amrex::ParticleReal const slice_bg = tau / m_betgam2;

            if (m_g > 0.0_prt)
            {
                // advance position and momentum (focusing quad)
                xout = cos_omega_ds*x + sin_omega_ds/m_omega*px;
                pxout = -m_omega*sin_omega_ds*x + cos_omega_ds*px;

                yout = cosh_omega_ds*y + sinh_omega_ds/m_omega*py;
                pyout = m_omega*sinh_omega_ds*y + cosh_omega_ds*py;

                tout = t + slice_bg*pt;
                // ptout = pt;
            } else if (m_g < 0.0_prt)
            {
                // advance position and momentum (defocusing quad)
                xout = cosh_omega_ds*x + sinh_omega_ds/m_omega*px;
                pxout = m_omega*sinh_omega_ds*x + cosh_omega_ds*px;

                yout = cos_omega_ds*y + sin_omega_ds/m_omega*py;
                pyout = -m_omega*sin_omega_ds*y + cos_omega_ds*py;

                tout = t + slice_bg*pt;
                // ptout = pt;
            } else {
                // advance position and momentum (zero strength = drift)
                xout = x + tau * px;
                // pxout = px;
                yout = y + tau * py;
                // pyout = py;
                tout = t + slice_bg * pt;
                // ptout = pt;
            }

            // push particle coordinates
            particle(1) = xout;
            particle(2) = pxout;
            particle(3) = yout;
            particle(4) = pyout;
            particle(5) = tout;
            particle(6) = ptout;

            zeval += 0_prt;
        }

        /** This pushes a particle in an ExactQuad element through
         *  the symplectic map associated with H_2 in the Hamiltonian
         *  splitting H = H_1 + H_2.
         *
         * @param tau Map step size in m
         * @param[in,out] particle particle phase space 6-vector
         * @param[in,out] zeval Longitudinal on-axis location in m
         * @param[in] refpart reference particle
         */
        template<typename T_Real>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void map2 (amrex::ParticleReal const tau,
                   amrex::SmallVector<T_Real, 6, 1> & particle,
                   amrex::ParticleReal & zeval,
                   [[maybe_unused]] RefPart const & AMREX_RESTRICT refpart) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using namespace std; // for cmath(float)
            using amrex::Math::powi;

            T_Real const x = particle(1);
            T_Real const px = particle(2);
            T_Real const y = particle(3);
            T_Real const py = particle(4);
            T_Real const t = particle(5);
            T_Real const pt = particle(6);

            // compute the radical in the denominator (= pz):
            T_Real const inv_pzden = 1_prt / sqrt(
                powi<2>(pt - m_ibeta) -
                m_ibetgam2 -
                powi<2>(px) -
                powi<2>(py)
            );

            // The momenta are not affected.
            particle(1) = x + tau * px * (-1_prt + inv_pzden);
            particle(2) = px;
            particle(3) = y + tau * py * (-1_prt + inv_pzden);
            particle(4) = py;
            particle(5) = t + tau * (-m_ibeta - m_ibetgam2*pt + inv_pzden*(1_prt - pt*m_beta)*m_ibeta);
            particle(6) = pt;

            zeval += tau;
        }

        /** This pushes particle spin in an ExactQuad element over a
         *  single step of size tau, based on the Thomas-BMT equation.
         *
         * @param tau Map step size in m
         * @param[in,out] particle particle phase space 6-vector
         * @param[in,out] particle_spin particle spin 3-vector
         * @param[in,out] zeval Longitudinal on-axis location in m
         * @param[in] refpart reference particle
         */
        template<typename T_Real>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void map3 (amrex::ParticleReal const tau,
                   amrex::SmallVector<T_Real, 6, 1> & particle,
                   amrex::SmallVector<T_Real, 3, 1> & particle_spin,
                   amrex::ParticleReal & zeval,
                   RefPart const & AMREX_RESTRICT refpart) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using namespace std; // for cmath(float)
            using amrex::Math::powi;

            // reference particle relativistic factors
            amrex::ParticleReal const gamma_ref = refpart.gamma();
            amrex::ParticleReal const gyro_anomaly = refpart.gyromagnetic_anomaly;

            // initialize the three components of the axis-angle vector
            T_Real lambdax = 0_prt;
            T_Real lambday = 0_prt;
            T_Real lambdaz = 0_prt;

            T_Real const x = particle(1);
            T_Real const px = particle(2);
            T_Real const y = particle(3);
            T_Real const py = particle(4);
            //T_Real const t = particle(5); // not used
            T_Real const pt = particle(6);
            T_Real sx = particle_spin(1);
            T_Real sy = particle_spin(2);
            T_Real sz = particle_spin(3);

            // Magnetic field normalized by q/mc:
            T_Real const Bx = m_g * y * m_beta * gamma_ref;
            T_Real const By = m_g * x * m_beta * gamma_ref;
            T_Real const Bz = 0_prt;

            // Electric field normalized by q/mc^2:
            T_Real const Ex = 0.0_prt;
            T_Real const Ey = 0.0_prt;
            T_Real const Ez = 0.0_prt;

            // Quantities required to evaluate the full Thomas-BMT precession vector.
            T_Real const Pnorm = sqrt(1_prt - 2_prt * pt * m_ibeta + pt*pt);
            T_Real const iPnorm = 1_prt/Pnorm;
            T_Real const Ps = sqrt(Pnorm*Pnorm - px*px - py*py);
            T_Real const gamma = gamma_ref * (1_prt - pt*m_beta);
            T_Real const ux = px * iPnorm;
            T_Real const uy = py * iPnorm;
            T_Real const uz = Ps * iPnorm;

            // Zero curvature in a quadrupole
            amrex::ParticleReal const h = 0.0_prt;

            // Evaluation of the full Thomas-BMT precession vector.
            tbmt_precession_vector(x,ux,uy,uz,gamma,h,gyro_anomaly,Bx,By,Bz,Ex,Ey,Ez,lambdax,lambday,lambdaz);

            // Generator of the spin rotation for a single step
            lambdax *= tau;
            lambday *= tau;
            lambdaz *= tau;

            // push the spin vector using the generator just determined
            rotate_spin(lambdax,lambday,lambdaz,sx,sy,sz);

            // update the spin variables
            particle_spin(1) = sx;
            particle_spin(2) = sy;
            particle_spin(3) = sz;

            // The coordinates and momenta are not affected.

            zeval += 0_prt;
        }

        /** This pushes the reference particle.
         *
         * @param[in,out] refpart reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (RefPart & AMREX_RESTRICT refpart) const {  // TODO: update as well, but needs more careful placement of calc_constants

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

        /** This function returns the linear transport map.
         *
         * @returns 6x6 transport matrix
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        Map6x6
        transport_map (RefPart const & AMREX_RESTRICT refpart) const  // TODO: update as well, but needs more careful placement of calc_constants
        {
            using namespace amrex::literals; // for _rt and _prt
            using amrex::Math::powi;

            // length of the current slice
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // access reference particle values to find beta*gamma^2
            amrex::ParticleReal const pt_ref = refpart.pt;
            amrex::ParticleReal const betgam2 = powi<2>(pt_ref) - 1.0_prt;

            // normalize quad units to MAD-X convention if needed
            amrex::ParticleReal const g = m_unit == 1 ? m_k / refpart.rigidity_Tm() : m_k;

            // compute phase advance per unit length in s (in rad/m)
            amrex::ParticleReal const omega = std::sqrt(std::abs(g));

            // initialize linear map matrix elements
            Map6x6 R = Map6x6::Identity();

            if (g > 0.0_prt) {
                R(1,1) = std::cos(omega*slice_ds);
                R(1,2) = std::sin(omega*slice_ds)/omega;
                R(2,1) = -omega*std::sin(omega*slice_ds);
                R(2,2) = std::cos(omega*slice_ds);
                R(3,3) = std::cosh(omega*slice_ds);
                R(3,4) = std::sinh(omega*slice_ds)/omega;
                R(4,3) = omega*std::sinh(omega*slice_ds);
                R(4,4) = std::cosh(omega*slice_ds);
                R(5,6) = slice_ds/betgam2;
            } else if (g < 0.0_prt) {
                R(1,1) = std::cosh(omega*slice_ds);
                R(1,2) = std::sinh(omega*slice_ds)/omega;
                R(2,1) = omega*std::sinh(omega*slice_ds);
                R(2,2) = std::cosh(omega*slice_ds);
                R(3,3) = std::cos(omega*slice_ds);
                R(3,4) = std::sin(omega*slice_ds)/omega;
                R(4,3) = -omega*std::sin(omega*slice_ds);
                R(4,4) = std::cos(omega*slice_ds);
                R(5,6) = slice_ds/betgam2;
            } else {
                R(1,2) = slice_ds;
                R(3,4) = slice_ds;
                R(5,6) = slice_ds / betgam2;
            }

            // apply the transverse rotation (roll) alignment error
            return rotate_aligned_map(R);
        }

        /** This operator pushes the particle spin.
         *
         * @param[in,out] x particle position in x
         * @param[in,out] y particle position in y
         * @param[in,out] t particle position in t
         * @param[in,out] px particle momentum in x
         * @param[in,out] py particle momentum in y
         * @param[in,out] pt particle momentum in t
         * @param[in,out] sx particle spin x-component
         * @param[in,out] sy particle spin y-component
         * @param[in,out] sz particle spin z-component
         * @param[in,out] idcpu particle global index (unused)
         * @param[in] refpart reference particle
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
            [[maybe_unused]] T_IdCpu const & AMREX_RESTRICT idcpu,
            RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // numerical integration parameters
            amrex::ParticleReal const zin = 0_prt;
            amrex::ParticleReal const zout = m_slice_ds;
            int const nsteps = m_mapsteps;

            // initialize phase space 6-vector
            amrex::SmallVector<T_Real, 6, 1> particle{
                x, px, y, py, t, pt
            };

            // initialize spin 3-vector
            amrex::SmallVector<T_Real, 3, 1> particle_spin{
                sx, sy, sz
            };

            // call integrator to advance through slice
            integrators::symp_integrate_particle_spin(particle,particle_spin,zin,zout,nsteps,m_int_order,refpart,*this);

            // assign updated values
            x = particle(1);
            px = particle(2);
            y = particle(3);
            py = particle(4);
            t = particle(5);
            pt = particle(6);
            sx = particle_spin(1);
            sy = particle_spin(2);
            sz = particle_spin(3);
        }

        /** This pushes the covariance matrix. */
        using LinearTransport::operator();

        amrex::ParticleReal m_k; //! quadrupole strength in 1/m^2 (or T/m)
        int m_unit; //! unit specification for quad strength
        int m_int_order; //! order used for the symplectic integration (2 or 4)
        int m_mapsteps; //! number of integration steps per slice

    private:
        // constants that are independent of the individually tracked particle,
        // see: compute_constants() to refresh
        amrex::ParticleReal m_slice_ds;       //! m_ds / nslice();
        amrex::ParticleReal m_betgam2;        //! beta*gamma^2
        amrex::ParticleReal m_ibetgam2;       //! 1 / m_betgam2
        amrex::ParticleReal m_beta;           //! beta
        amrex::ParticleReal m_ibeta;          //! 1 / m_beta
        amrex::ParticleReal m_g;              //! quadrupole strength in 1/m^2
        amrex::ParticleReal m_omega;
    };

} // namespace impactx

IMPACTX_PUSH_EXTERN_TEMPLATE(impactx::elements::ExactQuad)

#endif // IMPACTX_EXACTQUAD_H
