/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_EXACTMULTIPOLE_H
#define IMPACTX_EXACTMULTIPOLE_H

#include "particles/ImpactXParticleContainer.H"
#include "Element.H"
#include "particles/integrators/Integrators.H"
#include "mixin/alignment.H"
#include "mixin/beamoptic.H"
#include "mixin/devicebackedarrays.H"
#include "mixin/lineartransport.H"
#include "mixin/spintransport.H"
#include "mixin/pipeaperture.H"
#include "mixin/thick.H"

#include <AMReX_Extension.H>
#include <AMReX_GpuComplex.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>
#include <AMReX_SIMD.H>
#include <AMReX_TrackedVector.H>

#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

namespace impactx::elements
{

    /** Multipole coefficients for magnetic field.
     *
     * Host data persists across simulations; device data is uploaded
     * on-demand "lazily" on first use and released when AMReX shuts down.
     */
    struct MultipoleCoefficients
    {
        amrex::Gpu::TrackedVector<amrex::ParticleReal> k_normal;
        amrex::Gpu::TrackedVector<amrex::ParticleReal> k_skew;
    };

    /** The prepared execution snapshot of an exact multipole.
     *
     * Holds the particle-facing configuration, the constants derived from the reference
     * particle, and a bounded, non-owning view of the multipole coefficients.
     */
    struct ExactMultipoleCore
    : public mixin::ParticlePush<ExactMultipoleCore>,
      public mixin::Thick,
      public mixin::Alignment,
      public mixin::PipeAperture,
      public mixin::SpinTransport
#ifndef _WIN32
      // https://github.com/AMReX-Codes/amrex/issues/4609
      , public amrex::simd::Vectorized<amrex::simd::native_simd_size_particlereal>
#endif
    {
        static constexpr auto type = "ExactMultipole";
        using PType = ImpactXParticleContainer::ParticleType;

        static constexpr int DEFAULT_unit = 0;
        static constexpr int DEFAULT_int_order = 2;
        static constexpr int DEFAULT_mapsteps = 10;

        /** The particle-facing configuration of an exact multipole
         *
         * @param ds Segment length in m.
         * @param unit unit specification for the multipole strengths
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param aperture_x horizontal half-aperture in m
         * @param aperture_y vertical half-aperture in m
         * @param int_order order used for symplectic integration (2, 4, or 6)
         * @param mapsteps number of integration steps per slice
         * @param nslice number of slices used for the application of space charge
         */
        ExactMultipoleCore (
            amrex::ParticleReal ds,
            int unit,
            amrex::ParticleReal dx,
            amrex::ParticleReal dy,
            amrex::ParticleReal rotation_degree,
            amrex::ParticleReal aperture_x,
            amrex::ParticleReal aperture_y,
            int int_order,
            int mapsteps,
            int nslice
        )
          : Thick(ds, nslice),
            Alignment(dx, dy, rotation_degree),
            PipeAperture(aperture_x, aperture_y),
            m_unit(unit), m_int_order(int_order), m_mapsteps(mapsteps)
        {
        }

        /** Reverse the element (negate ds) */
        void reverse () { Thick::reverse(); }

        /** Push all particles on a tile */
        using mixin::ParticlePush<ExactMultipoleCore>::operator();

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
            m_brho = refpart.rigidity_Tm();
        }

        /** This is a Multipole functor, so that a variable of this type can be used like a Multipole function.
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
            [[maybe_unused]] T_IdCpu const & AMREX_RESTRICT idcpu,
            RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // coordinates within the multipole (not needed)
            //amrex::ParticleReal const s = refpart.s;
            //amrex::ParticleReal const sedge = refpart.sedge;

            // numerical integration parameters
            amrex::ParticleReal const zin = 0_prt;
            amrex::ParticleReal const zout = m_slice_ds;
            int const nsteps = m_mapsteps;

            // initialize phase space 6-vector
            amrex::SmallVector<T_Real, 6, 1> particle{
                x, px, y, py, t, pt
            };

            // call integrator to advance through slice (int_order = 2 or 4, otherwise default 2)
            if (m_int_order == 2) {
                integrators::symp2_integrate_particle(particle,zin,zout,nsteps,refpart,*this);
            } else if (m_int_order == 4) {
                integrators::symp4_integrate_particle(particle,zin,zout,nsteps,refpart,*this);
            } else if (m_int_order == 6) {
                integrators::symp6_integrate_particle(particle,zin,zout,nsteps,refpart,*this);
            } else {
                integrators::symp2_integrate_particle(particle,zin,zout,nsteps,refpart,*this);
            }

            // assign updated values
            x = particle(1);
            px = particle(2);
            y = particle(3);
            py = particle(4);
            t = particle(5);
            pt = particle(6);
        }

        /** This pushes a particle in an ExactMultipole element through
         *  the symplectic map associated with H_1 in the Hamiltonian
         *  splitting H = H_1 + H_2.  This is the exact map for a drift.
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
            using amrex::Math::powi;
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

            // compute the radical in the denominator (= pz):
            T_Real const inv_pzden = 1_prt / sqrt(
                powi<2>(pt - m_ibeta) -
                m_ibetgam2 -
                powi<2>(px) -
                powi<2>(py)
            );

            // advance position and momentum (exact drift)
            xout = x + tau * px * inv_pzden;
            // pxout = px;
            yout = y + tau * py * inv_pzden;
            // pyout = py;
            tout = t - tau * (m_ibeta + (pt - m_ibeta) * inv_pzden);
            // ptout = pt;

            // push particle coordinates
            particle(1) = xout;
            particle(2) = pxout;
            particle(3) = yout;
            particle(4) = pyout;
            particle(5) = tout;
            particle(6) = ptout;

            zeval += 0_prt;
        }

        /** This pushes a particle in an ExactMultipole element through
         *  the symplectic map associated with H_2 in the Hamiltonian
         *  splitting H = H_1 + H_2.  This is the exact kick obtained
         *  in the ideal multipole field.
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

            amrex::ParticleReal const * const k_normal = m_coefficients.first;
            amrex::ParticleReal const * const k_skew = m_coefficients.second;

            // a complex type with two T_Real
            using Complex = amrex::GpuComplex<T_Real>;

            // assign complex position and complex multipole strength
            Complex const zeta(x, y);

            // initialize the momentum kick
            T_Real dpx = 0_prt;
            T_Real dpy = 0_prt;

            int factorial = 1;

            // loop over array of multipole coefficients
            for (int m = 0; m < m_coefficients.size; m++) {

               // normalize units to MAD-X convention if needed
               amrex::ParticleReal knormal_value = m_unit == 1 ? k_normal[m] / m_brho : k_normal[m];
               amrex::ParticleReal kskew_value = m_unit == 1 ? k_skew[m] / m_brho : k_skew[m];

               // multipole contribution to complex momentum kick
               Complex alpha(knormal_value, kskew_value);
               Complex kick = amrex::pow(zeta, m) * alpha;
               dpx += -1.0_prt * kick.m_real / factorial;
               dpy += kick.m_imag / factorial;

               factorial *= m+1;

            }

            // advance momentum by a step of length tau
            pxout = pxout + tau * dpx;
            pyout = pyout + tau * dpy;

            //update the particle coordinates
            particle(1) = xout;
            particle(2) = pxout;
            particle(3) = yout;
            particle(4) = pyout;
            particle(5) = tout;
            particle(6) = ptout;

            zeval += tau;
        }

        /** This pushes particle spin in an ExactMultipole element over a
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

            amrex::ParticleReal const * const k_normal = m_coefficients.first;
            amrex::ParticleReal const * const k_skew = m_coefficients.second;

            // a complex type with two T_Real
            using Complex = amrex::GpuComplex<T_Real>;

            // assign complex position and complex multipole strength
            Complex const zeta(x, y);

            // Initialization of magnetic field variables:
            T_Real Bx = 0_prt;
            T_Real By = 0_prt;
            T_Real Bz = 0_prt;

            int factorial = 1;

            // loop over array of multipole coefficients
            for (int m = 0; m < m_coefficients.size; m++) {

               // normalize units to MAD-X convention if needed
               amrex::ParticleReal knormal_value = m_unit == 1 ? k_normal[m] / m_brho : k_normal[m];
               amrex::ParticleReal kskew_value = m_unit == 1 ? k_skew[m] / m_brho : k_skew[m];

               // multipole contribution to magnetic field normalized by q/mc:
               Complex alpha(knormal_value, kskew_value);
               Complex kick = amrex::pow(zeta, m) * alpha;
               Bx += kick.m_imag / factorial * m_beta * gamma_ref;
               By += kick.m_real / factorial * m_beta * gamma_ref;

               factorial *= m+1;

            }

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

            // Zero curvature in a straight-axis multipole
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

        int m_unit; //! unit specification for Multipole strength
        int m_int_order; //! order used for the symplectic integration (2 or 4)
        int m_mapsteps; //! number of integration steps per slice
        /** Coefficients resolved to this snapshot's execution space; non-owning */
        mixin::DualArrayView m_coefficients{};


    private:
        // constants that are independent of the individually tracked particle,
        // see: compute_constants() to refresh
        amrex::ParticleReal m_slice_ds;       //! m_ds / nslice();
        amrex::ParticleReal m_betgam2;        //! beta*gamma^2
        amrex::ParticleReal m_ibetgam2;       //! 1 / m_betgam2
        amrex::ParticleReal m_beta;           //! beta
        amrex::ParticleReal m_ibeta;          //! 1 / m_beta
        amrex::ParticleReal m_brho;           //! magnetic ridigity in T-m
    };

    /** An exact multipole: the canonical element that owns its coefficients.
     *
     * The reference-particle push, the transfer map, serialization and finalization run
     * here on the host. Particle pushes run against @see ExactMultipoleCore.
     */
    struct ExactMultipolePhysics
    : public mixin::BeamOptic<ExactMultipolePhysics>,
      public mixin::LinearTransport<ExactMultipolePhysics>,
      public ExactMultipoleCore
    {
        /** A general (thick) straight-axis magnetic element whose field consists
         *  of a superposition of ideal hard-edge multipoles.  Like the element
         *  ExactQuad, the magnet is treated using the exact Hamiltonian, which
         *  includes all nonlinear kinematic effects.  This is a generalization
         *  of the (thin) element Multipole, which applies the kick associated
         *  with a single ideal multipole.  Here, symplectic integration
         *  is used based on a splitting of the underlying Hamiltonian.
         *
         * @param ds Segment length in m.
         * @param k_normal Array of normal multipole coefficients (units below)
         * @param k_skew   Array of skew multipole coefficients (units below)
         * @param unit  Unit specification
         *           unit = 0 multipole strengths in units of 1/m^n (default)
         *           unit = 1 multipole strengths in units of T/m^(n-1)
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param aperture_x horizontal half-aperture in m
         * @param aperture_y vertical half-aperture in m
         * @param int_order order used for symplectic integration (2, 4, or 6)
         * @param mapsteps number of integration steps per slice used for
         *        particle push in the applied fields
         * @param nslice number of slices used for the application of space charge
         */
        ExactMultipolePhysics (
            amrex::ParticleReal ds,
            std::vector<amrex::ParticleReal> k_normal,
            std::vector<amrex::ParticleReal> k_skew,
            int unit = DEFAULT_unit,
            amrex::ParticleReal dx = DEFAULT_dx,
            amrex::ParticleReal dy = DEFAULT_dy,
            amrex::ParticleReal rotation_degree = DEFAULT_rotation_degree,
            amrex::ParticleReal aperture_x = DEFAULT_aperture_x,
            amrex::ParticleReal aperture_y = DEFAULT_aperture_y,
            int int_order = DEFAULT_int_order,
            int mapsteps = DEFAULT_mapsteps,
            int nslice = DEFAULT_nslice
        )
          : ExactMultipoleCore(ds, unit, dx, dy, rotation_degree, aperture_x, aperture_y, int_order, mapsteps, nslice),
            m_coefficients_data(std::move(k_normal), std::move(k_skew), "ExactMultipole: normal and skew coefficient")
        {
            if (m_int_order != 2 && m_int_order != 4 && m_int_order != 6)
                throw std::invalid_argument("ExactMultipole: The order used for symplectic integration must be 2, 4 or 6.");

            if (m_mapsteps <= 0)
                throw std::invalid_argument("ExactMultipole: mapsteps must be > 0!");
        }

        /** Push all particles */
        using mixin::BeamOptic<ExactMultipolePhysics>::operator();

        /** This pushes the covariance matrix. */
        using mixin::LinearTransport<ExactMultipolePhysics>::operator();

        /** This pushes the reference particle.
         *
         * @param[in,out] refpart reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (RefPart & AMREX_RESTRICT refpart) const  // TODO: update as well, but needs more careful placement of calc_constants
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

            // pick the right data depending if we are on the host side
            // (reference particle push) or device side (particles):
            // host view, obtained at the point of use
            auto const coefficients = m_coefficients_data.host_view();
            amrex::ParticleReal const * k_normal = coefficients.first;
            amrex::ParticleReal const * k_skew = coefficients.second;

            // normalize quad units to MAD-X convention if needed
            // a pure dipole (or an element with no coefficients at all) has no quadrupole
            // component, and then supplies no index 1
            amrex::ParticleReal const kn_raw = coefficients.size > 1 ? k_normal[1] : 0.0_prt;
            amrex::ParticleReal const ks_raw = coefficients.size > 1 ? k_skew[1] : 0.0_prt;
            amrex::ParticleReal const kn = m_unit == 1 ? kn_raw / refpart.rigidity_Tm() : kn_raw;
            amrex::ParticleReal const ks = m_unit == 1 ? ks_raw / refpart.rigidity_Tm() : ks_raw;
            amrex::ParticleReal const kabs = std::sqrt(powi<2>(kn) + powi<2>(ks));

            // initialize linear map matrix elements
            Map6x6 R = Map6x6::Identity();

            // No quadrupole component (e.g. a pure sextupole or higher
            // multipole, or a pure dipole): the first-order map about the
            // reference orbit reduces to a drift.
            if (kabs == 0.0_prt) {
                R(1,2) = slice_ds;
                R(3,4) = slice_ds;
                R(5,6) = slice_ds / betgam2;
                // apply the transverse rotation (roll) alignment error
                return rotate_aligned_map(R);
            }

            amrex::ParticleReal const kplus = kn + kabs;
            amrex::ParticleReal const kminus = kn - kabs;

            // compute phase advance per unit length in s (in rad/m)
            amrex::ParticleReal const omega = std::sqrt(kabs);
            amrex::ParticleReal const sin_omega_ds = std::sin(omega*slice_ds);
            amrex::ParticleReal const cos_omega_ds = std::cos(omega*slice_ds);
            amrex::ParticleReal const sinh_omega_ds = std::sinh(omega*slice_ds);
            amrex::ParticleReal const cosh_omega_ds = std::cosh(omega*slice_ds);

            // linear matrix for a quadrupole with both normal and skew coefficients
            R(1,1) = (kplus * cos_omega_ds - kminus * cosh_omega_ds) / (2_prt * kabs);
            R(1,2) = (kplus * sin_omega_ds - kminus * sinh_omega_ds) / (2_prt * kabs * omega);
            R(1,3) = ks * (-cos_omega_ds + cosh_omega_ds) / (2_prt * kabs);
            R(1,4) = ks * (-sin_omega_ds + sinh_omega_ds) / (2_prt * kabs * omega);
            R(2,1) = (-kplus * sin_omega_ds - kminus * sinh_omega_ds) / (2_prt * omega);
            R(2,2) = (kplus * cos_omega_ds - kminus * cosh_omega_ds) / (2_prt * kabs);
            R(2,3) = ks * (sin_omega_ds + sinh_omega_ds) / (2_prt * omega);
            R(2,4) = ks * (-cos_omega_ds + cosh_omega_ds) / (2_prt * kabs);
            R(3,1) = ks * (-cos_omega_ds + cosh_omega_ds) / (2_prt * kabs);
            R(3,2) = ks * (-sin_omega_ds + sinh_omega_ds) / (2_prt * kabs * omega);
            R(3,3) = (-kminus * cos_omega_ds + kplus * cosh_omega_ds) / (2_prt * kabs);
            R(3,4) = (-kminus * sin_omega_ds + kplus * sinh_omega_ds) / (2_prt * kabs * omega);
            R(4,1) = ks * (sin_omega_ds + sinh_omega_ds) / (2_prt * omega);
            R(4,2) = ks * (-cos_omega_ds + cosh_omega_ds) / (2_prt * kabs);
            R(4,3) = (kminus * sin_omega_ds + kplus * sinh_omega_ds) / (2_prt * omega);
            R(4,4) = (-kminus * cos_omega_ds + kplus * cosh_omega_ds) / (2_prt * kabs);
            R(5,6) = slice_ds/betgam2;

            /* Debugging:
            for (int j=1; j<=6; j++) {
               for (int k=1; k<=6; k++) {
                  amrex::Print() << "R(" << j << "," << k << ") = " << R(j,k) << "\n";
               }
            }
            */

            // apply the transverse rotation (roll) alignment error
            return rotate_aligned_map(R);
        }

        /** Prepare the execution snapshot for one element slice
         *
         * @param refpart reference particle, already advanced for this slice
         * @return a core holding an execution-space view of the coefficients
         */
        [[nodiscard]] ExactMultipoleCore make_particle_pusher (RefPart const & refpart) const
        {
            ExactMultipoleCore core = *this;
            core.compute_constants(refpart);
            core.m_coefficients = m_coefficients_data.device_view();
            return core;
        }

        /** Normal multipole coefficients */
        [[nodiscard]] std::vector<amrex::ParticleReal> k_normal_coefficients () const
        {
            return m_coefficients_data.first_host();
        }

        /** Skew multipole coefficients */
        [[nodiscard]] std::vector<amrex::ParticleReal> k_skew_coefficients () const
        {
            return m_coefficients_data.second_host();
        }

        /** Replace both arrays at once
         *
         * The two are validated together and committed as one step, so a rejected update
         * leaves the element as it was. They are always the same length, which is why
         * this takes both rather than offering to change one of them alone.
         *
         * @param first the new normal values
         * @param second the new skew values
         */
        void set_coefficients (
            std::vector<amrex::ParticleReal> first,
            std::vector<amrex::ParticleReal> second
        )
        {
            m_coefficients_data.commit(std::move(first), std::move(second), "ExactMultipole: normal and skew coefficient");
        }

        /** Release device storage; host configuration survives for reuse */
        void finalize ()
        {
            m_coefficients_data.release_gpu();
        }

        //! the coefficient arrays, owned by this element
        mixin::DeviceBackedArrays m_coefficients_data;
    };

    /** ExactMultipole: physics plus common host-only metadata */
    using ExactMultipole = Element<ExactMultipolePhysics>;

} // namespace impactx

IMPACTX_PUSH_EXTERN_TEMPLATE(impactx::elements::ExactMultipolePhysics)

#endif // IMPACTX_EXACTMULTIPOLE_H
