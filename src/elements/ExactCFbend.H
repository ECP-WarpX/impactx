/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_EXACTCFBEND_H
#define IMPACTX_EXACTCFBEND_H

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
#include <AMReX_REAL.H>
#include <AMReX_Print.H>
#include <AMReX_GpuComplex.H>
#include <AMReX_Math.H>
#include <AMReX_TrackedVector.H>

#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

namespace impactx::elements
{

    /** CFbend multipole coefficients for magnetic field.
     *
     * Host data persists across simulations; device data is uploaded
     * on-demand "lazily" on first use and released when AMReX shuts down.
     */
    struct CFbendCoefficients
    {
        amrex::Gpu::TrackedVector<amrex::ParticleReal> k_normal;
        amrex::Gpu::TrackedVector<amrex::ParticleReal> k_skew;
    };

    /** The prepared execution snapshot of an exact combined-function bend.
     *
     * Holds the particle-facing configuration, the constants derived from the reference
     * particle, and a bounded, non-owning view of the multipole coefficients.
     */
    struct ExactCFbendCore
    : public mixin::ParticlePush<ExactCFbendCore>,
      public mixin::Thick,
      public mixin::Alignment,
      public mixin::PipeAperture,
      public mixin::SpinTransport
    {
        static constexpr auto type = "ExactCFbend";
        using PType = ImpactXParticleContainer::ParticleType;

        static constexpr int DEFAULT_unit = 0;
        static constexpr int DEFAULT_int_order = 2;
        static constexpr int DEFAULT_mapsteps = 10;

        /** The particle-facing configuration of an exact combined-function bend
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
        ExactCFbendCore (
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
        using mixin::ParticlePush<ExactCFbendCore>::operator();

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

            Alignment::compute_constants(refpart);

            // access reference particle values to find relativistic factors
            m_beta = refpart.beta();
            m_ibeta = 1_prt / m_beta;
            m_ibetgam2 = 1_prt / amrex::Math::powi<2>(refpart.beta_gamma());
            m_brho = refpart.rigidity_Tm();

            amrex::ParticleReal const * k_normal = m_coefficients.first;

            // find magnetic field and curvature; normalize units to MAD-X convention if needed
            amrex::ParticleReal curvature = m_unit == 1 ? k_normal[0] / m_brho : k_normal[0];
            m_rc = curvature == 0.0 ? 0.0_prt : 1.0_prt/curvature;

            // arc length and angle of the current slice
            m_slice_ds = m_ds / nslice();
        }

        /** This is an ExactCFbend functor, so that a variable of this type can be used like an ExactCFbend function.
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
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
            amrex::ParticleReal & AMREX_RESTRICT x,
            amrex::ParticleReal & AMREX_RESTRICT y,
            amrex::ParticleReal & AMREX_RESTRICT t,
            amrex::ParticleReal & AMREX_RESTRICT px,
            amrex::ParticleReal & AMREX_RESTRICT py,
            amrex::ParticleReal & AMREX_RESTRICT pt,
            [[maybe_unused]] uint64_t const & AMREX_RESTRICT idcpu,
            RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // numerical integration parameters
            amrex::ParticleReal const zin = 0_prt;
            amrex::ParticleReal const zout = m_slice_ds;
            int const nsteps = m_mapsteps;

            // initialize phase space 6-vector
            amrex::SmallVector<amrex::ParticleReal, 6, 1> particle{
                x, px, y, py, t, pt
            };

            // call integrator to advance through slice (int_order = 2 or 4, otherwise default 2)
            if (m_int_order == 2) {
                integrators::symp2_integrate_particle(particle,zin,zout,nsteps,refpart,*this);
            } else if (m_int_order == 4) {
                integrators::symp4_integrate_particle(particle,zin,zout,nsteps,refpart,*this);
            } else if (m_int_order == 6) {
                integrators::symp6_integrate_particle(particle,zin,zout,nsteps,refpart,*this);
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
         *  splitting H = H_1 + H_2.  This is the exact map for a sector bend.
         *
         * @param tau Map step size in m
         * @param[in,out] particle particle phase space 6-vector
         * @param[in,out] zeval Longitudinal on-axis location in m
         * @param[in] refpart reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void map1 (amrex::ParticleReal const tau,
                   amrex::SmallVector<amrex::ParticleReal, 6, 1> & particle,
                   amrex::ParticleReal & zeval,
                   [[maybe_unused]] RefPart const & AMREX_RESTRICT refpart) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using namespace amrex::Math;  // for powi

            amrex::ParticleReal const x = particle(1);
            amrex::ParticleReal const px = particle(2);
            amrex::ParticleReal const y = particle(3);
            amrex::ParticleReal const py = particle(4);
            amrex::ParticleReal const t = particle(5);
            amrex::ParticleReal const pt = particle(6);

            amrex::ParticleReal xout = x;
            amrex::ParticleReal pxout = px;
            amrex::ParticleReal yout = y;
            amrex::ParticleReal pyout = py;
            amrex::ParticleReal tout = t;
            amrex::ParticleReal ptout = pt;

            // treat the special case of zero field
            if (m_rc==0_prt) {
               // compute the radical in the denominator (= pz):
               amrex::ParticleReal const inv_pzden = 1_prt / std::sqrt(
                   powi<2>(pt - m_ibeta) -
                   m_ibetgam2 -
                   powi<2>(px) -
                   powi<2>(py)
               );

               // advance position and momentum (exact drift)
               xout = x + tau * px * inv_pzden;
               yout = y + tau * py * inv_pzden;
               tout = t - tau * (m_ibeta + (pt - m_ibeta) * inv_pzden);

               // assign updated momenta
               pxout = px;
               pyout = py;
               ptout = pt;

            } else {

               // assign intermediate quantities
               amrex::ParticleReal const pperp2 = powi<2>(pt)-2.0_prt * m_ibeta * pt - powi<2>(py)+1.0_prt;
               amrex::ParticleReal const px2 = powi<2>(px);

               // trigonometric evaluations
               amrex::ParticleReal const slice_phi = tau / m_rc;
               auto const [sin_phi, cos_phi] = amrex::Math::sincos(slice_phi);

               // determine if particle lies within the domain of map definition
               if (pperp2 > px2)
               {
                   amrex::ParticleReal const pperp = std::sqrt(pperp2);
                   amrex::ParticleReal const pzi = std::sqrt(powi<2>(pperp) - powi<2>(px));
                   amrex::ParticleReal const rho = m_rc + x;

                   // update momenta
                   pxout = px * cos_phi + (pzi - rho / m_rc) * sin_phi;
                   pyout = py;
                   ptout = pt;

                   // angle of momentum rotation
                   amrex::ParticleReal const px2f = powi<2>(pxout);
                   // determine if particle lies within domain of map definition
                   if (pperp2 > px2f)
                   {
                       amrex::ParticleReal const pzf = std::sqrt(powi<2>(pperp)-powi<2>(pxout));
                       amrex::ParticleReal const theta = slice_phi + std::asin(px/pperp) - std::asin(pxout/pperp);

                       // update position coordinates
                       xout = -m_rc + rho*cos_phi + m_rc*(pzf + px*sin_phi - pzi*cos_phi);
                       yout = y + theta * m_rc * py;
                       tout = t - theta * m_rc * (pt - 1.0_prt * m_ibeta) - slice_phi * m_rc * m_ibeta;

                   } // else { amrex::Print() << "Warning: outside map domain." << "\n"; }
               } // else { amrex::Print() << "Warning: outside map domain." << "\n"; }
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
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void map2 (amrex::ParticleReal const tau,
                   amrex::SmallVector<amrex::ParticleReal, 6, 1> & particle,
                   amrex::ParticleReal & zeval,
                   [[maybe_unused]] RefPart const & AMREX_RESTRICT refpart) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using namespace amrex::Math;  // for powi

            amrex::ParticleReal const x = particle(1);
            amrex::ParticleReal const px = particle(2);
            amrex::ParticleReal const y = particle(3);
            amrex::ParticleReal const py = particle(4);
            amrex::ParticleReal const t = particle(5);
            amrex::ParticleReal const pt = particle(6);

            amrex::ParticleReal xout = x;
            amrex::ParticleReal pxout = px;
            amrex::ParticleReal yout = y;
            amrex::ParticleReal pyout = py;
            amrex::ParticleReal tout = t;
            amrex::ParticleReal ptout = pt;

            amrex::ParticleReal const * k_normal = m_coefficients.first;
            amrex::ParticleReal const * k_skew = m_coefficients.second;

            // defined dimensionless variables scaled by radius of curvature
            amrex::ParticleReal rho = 1_prt + x / m_rc;
            amrex::ParticleReal y_scl = y / m_rc;
            amrex::ParticleReal ln_rho = std::log(rho);

            // initialize the momentum kick
            amrex::ParticleReal dpx = 0_prt;
            amrex::ParticleReal dpy = 0_prt;

            // loop over array of multipole coefficients
            for (int m = 1; m < m_coefficients.size; m++) {

               // normalize units to MAD-X convention if needed
               amrex::ParticleReal kn = m_unit == 1 ? k_normal[m] / m_brho : k_normal[m];
               amrex::ParticleReal ks = m_unit == 1 ? k_skew[m] / m_brho : k_skew[m];

               int n = m+1;  // multipole order

               // define dimensionless multipole coefficients (scaled by radius of curvature)
               kn = std::pow(m_rc, static_cast<amrex::ParticleReal>(m)) * kn;
               ks = std::pow(m_rc, static_cast<amrex::ParticleReal>(m)) * ks;

               switch (n) {

                  case 1:   // dipole-like contribution (not used - contributing to map 1)
                     dpx += -rho * kn;
                     dpy += ks;
                     break;
                  case 2:   // quadrupole-like contribution
                     dpx += -rho*ln_rho*kn + y_scl*rho*ks;
                     dpy += y_scl*kn + 0.5_prt*(rho*rho - 1_prt)*ks;
                     break;
                  case 3:  // sextupole-like contribution
                     dpx += 0.25_prt*rho*(1_prt + 2_prt*y_scl*y_scl - rho*rho + 2_prt*ln_rho)*kn + y_scl*rho*ln_rho*ks;
                     dpy += 0.5_prt*y_scl*(rho*rho - 1_prt)*kn + 0.25_prt*(1_prt - 2_prt*y_scl*y_scl - rho*rho + 2_prt*rho*rho*ln_rho)*ks;
                     break;
                  case 4:  // octupole-like contribution
                     dpx += -0.25_prt*(rho - powi<3>(rho) + (rho - 2_prt*y_scl*y_scl*rho + powi<3>(rho))*ln_rho)*kn
                            -1_prt/12_prt*y_scl*rho*(3_prt + 2_prt*y_scl*y_scl - 3_prt*rho*rho + 6_prt*ln_rho)*ks;
                     dpy += -1_prt/12_prt*y_scl*(-3_prt + 2_prt*y_scl*y_scl + 3_prt*rho*rho - 6_prt*rho*rho*ln_rho)*kn
                            -1_prt/16_prt*((-1_prt + 4_prt*y_scl*y_scl - rho*rho)*(-1_prt + rho*rho) + 4_prt*rho*rho*ln_rho)*ks;
                     break;
                  case 5:  // decapole-like contribution
                     dpx += -1_prt/192_prt*rho*(8_prt*powi<4>(y_scl) - 24_prt*powi<2>(y_scl)*(-1_prt + rho*rho)
                            + 3_prt*(-5_prt + 4_prt*powi<2>(rho) + powi<4>(rho))
                            +12_prt*(-1_prt + 4_prt*y_scl*y_scl - 2_prt*rho*rho)*ln_rho)*kn
                            +1_prt/12_prt*y_scl*rho*(3_prt*(-1_prt + rho*rho) + (2_prt*y_scl*y_scl - 3_prt*(1_prt + rho*rho))*ln_rho)*ks;
                     dpy += -1_prt/48_prt*y_scl*((-1_prt + rho*rho)*(4_prt*y_scl*y_scl - 3_prt*(1_prt + rho*rho)) + 12_prt*rho*rho*ln_rho)*kn
                            + 1_prt/192_prt*(3_prt + 8_prt*powi<4>(y_scl) + 12_prt*rho*rho
                            - 15_prt*powi<4>(rho) + 24_prt*y_scl*y_scl*(-1_prt + rho*rho)
                            +12_prt*rho*rho*(2_prt - 4_prt*y_scl*y_scl + rho*rho)*ln_rho)*ks;
                     break;
                // default:
                //   TODO: print a warning here that order > 5 is not currently supported
               }

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

        /** This pushes particle spin in an ExactCFbend element over a
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

            // defined dimensionless variables scaled by radius of curvature
            amrex::ParticleReal rho = 1_prt + x / m_rc;
            amrex::ParticleReal y_scl = y / m_rc;
            amrex::ParticleReal ln_rho = std::log(rho);

            // Initialization of magnetic field variables:
            T_Real Bx = 0_prt;
            T_Real By = 0_prt;
            T_Real Bz = 0_prt;

            // loop over array of multipole coefficients
            // NOTE:  Here we start at m = 0 in order to include the dipole contribution to the spin
            for (int m = 0; m < m_coefficients.size; m++) {

               // normalize units to MAD-X convention if needed
               amrex::ParticleReal kn = m_unit == 1 ? k_normal[m] / m_brho : k_normal[m];
               amrex::ParticleReal ks = m_unit == 1 ? k_skew[m] / m_brho : k_skew[m];

               int n = m+1;  // multipole order

               // define dimensionless multipole coefficients (scaled by radius of curvature)
               kn = std::pow(m_rc, static_cast<amrex::ParticleReal>(m)) * kn;
               ks = std::pow(m_rc, static_cast<amrex::ParticleReal>(m)) * ks;

               // Below, we use the consistency relationships Bx = dpy/rho and By = -dpx/rho.  See T. Zolkin, PRAB 20, 043501 (2017).

               switch (n) {

                  case 1:   // dipole-like contribution (not used - contributing to map 1)
                     Bx += ks;
                     By += rho * kn;
                     break;
                  case 2:   // quadrupole-like contribution
                     Bx += y_scl*kn + 0.5_prt*(rho*rho - 1_prt)*ks;
                     By += -(-rho*ln_rho*kn + y_scl*rho*ks);
                     break;
                  case 3:  // sextupole-like contribution
                     Bx += 0.5_prt*y_scl*(rho*rho - 1_prt)*kn + 0.25_prt*(1_prt - 2_prt*y_scl*y_scl - rho*rho + 2_prt*rho*rho*ln_rho)*ks;
                     By += -(0.25_prt*rho*(1_prt + 2_prt*y_scl*y_scl - rho*rho + 2_prt*ln_rho)*kn + y_scl*rho*ln_rho*ks);
                     break;
                  case 4:  // octupole-like contribution
                     Bx += -1_prt/12_prt*y_scl*(-3_prt + 2_prt*y_scl*y_scl + 3_prt*rho*rho - 6_prt*rho*rho*ln_rho)*kn
                            -1_prt/16_prt*((-1_prt + 4_prt*y_scl*y_scl - rho*rho)*(-1_prt + rho*rho) + 4_prt*rho*rho*ln_rho)*ks;
                     By += -(-0.25_prt*(rho - powi<3>(rho) + (rho - 2_prt*y_scl*y_scl*rho + powi<3>(rho))*ln_rho)*kn
                            -1_prt/12_prt*y_scl*rho*(3_prt + 2_prt*y_scl*y_scl - 3_prt*rho*rho + 6_prt*ln_rho)*ks);
                     break;
                  case 5:  // decapole-like contribution
                     Bx += -1_prt/48_prt*y_scl*((-1_prt + rho*rho)*(4_prt*y_scl*y_scl - 3_prt*(1_prt + rho*rho)) + 12_prt*rho*rho*ln_rho)*kn
                            + 1_prt/192_prt*(3_prt + 8_prt*powi<4>(y_scl) + 12_prt*rho*rho
                            - 15_prt*powi<4>(rho) + 24_prt*y_scl*y_scl*(-1_prt + rho*rho)
                            +12_prt*rho*rho*(2_prt - 4_prt*y_scl*y_scl + rho*rho)*ln_rho)*ks;
                     By += -(-1_prt/192_prt*rho*(8_prt*powi<4>(y_scl) - 24_prt*powi<2>(y_scl)*(-1_prt + rho*rho)
                            + 3_prt*(-5_prt + 4_prt*powi<2>(rho) + powi<4>(rho))
                            +12_prt*(-1_prt + 4_prt*y_scl*y_scl - 2_prt*rho*rho)*ln_rho)*kn
                            +1_prt/12_prt*y_scl*rho*(3_prt*(-1_prt + rho*rho) + (2_prt*y_scl*y_scl - 3_prt*(1_prt + rho*rho))*ln_rho)*ks);
                     break;
               }

            }

            // Scaling magnetic field to units used in the T-BMT representation (e.g. q*B/mc).

            Bx = Bx * m_beta * gamma_ref / rho;
            By = By * m_beta * gamma_ref / rho;

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

            // Design-orbit curvature
            amrex::ParticleReal const h = 1_prt/m_rc;

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


    protected:
        // constants that are independent of the individually tracked particle,
        // see: compute_constants() to refresh.
        // protected: the owning element reads the rigidity in its transfer map.
        amrex::ParticleReal m_slice_ds;       //! m_ds / nslice();
        amrex::ParticleReal m_ibetgam2;       //! 1 / m_betgam2
        amrex::ParticleReal m_beta;           //! beta
        amrex::ParticleReal m_ibeta;          //! 1 / m_beta
        amrex::ParticleReal m_brho;           //! magnetic ridigity in T-m
        amrex::ParticleReal m_rc;             //! magnet radius of curvature
    };

    /** An exact combined-function bend: the canonical element that owns its coefficients.
     *
     * The reference-particle push, the transfer map, serialization and finalization run
     * here on the host. Particle pushes run against @see ExactCFbendCore.
     */
    struct ExactCFbendPhysics
    : public mixin::BeamOptic<ExactCFbendPhysics>,
      public mixin::LinearTransport<ExactCFbendPhysics>,
      public ExactCFbendCore
    {
        /** A general (thick) combined-function sector dipole whose field consists
         *  of a superposition of ideal curvilinear hard-edge multipoles.
         *  The magnet is treated using the exact Hamiltonian, which
         *  includes all nonlinear kinematic effects.  This is a generalization
         *  of the element CFbend, which is a linear sector bend that
         *  supports quadrupole focusing.  Here, symplectic integration
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
        ExactCFbendPhysics (
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
          : ExactCFbendCore(ds, unit, dx, dy, rotation_degree, aperture_x, aperture_y, int_order, mapsteps, nslice),
            m_coefficients_data(std::move(k_normal), std::move(k_skew), "ExactCFbend: normal and skew coefficient")
        {
            if (m_int_order != 2 && m_int_order != 4 && m_int_order != 6)
                throw std::invalid_argument("ExactCFbend: The order used for symplectic integration must be 2, 4 or 6.");

            if (m_mapsteps <= 0)
                throw std::invalid_argument("ExactCFbend: mapsteps must be > 0!");
        }

        /** Push all particles */
        using mixin::BeamOptic<ExactCFbendPhysics>::operator();

        /** This pushes the covariance matrix. */
        using mixin::LinearTransport<ExactCFbendPhysics>::operator();

        /** This pushes the reference particle.
         *
         * @param[in,out] refpart reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (RefPart & AMREX_RESTRICT refpart) const
        {
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
            amrex::ParticleReal const brho = refpart.rigidity_Tm();

            // host view, obtained at the point of use
            amrex::ParticleReal const * k_normal = m_coefficients_data.host_view().first;
            // find magnetic field and curvature; normalize units to MAD-X convention if needed
            amrex::ParticleReal const curvature = m_unit == 1 ? k_normal[0] / brho : k_normal[0];
            amrex::ParticleReal const rc = curvature == 0.0 ? 0.0_prt : 1.0_prt/curvature;

            // length of the current slice
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // special case of zero field = an exact drift
            if (rc==0_prt) {
               // advance position and momentum (drift)
               amrex::ParticleReal const step = slice_ds / std::sqrt(amrex::Math::powi<2>(pt) - 1.0_prt);
               refpart.x = x + step*px;
               refpart.y = y + step*py;
               refpart.z = z + step*pz;
               refpart.t = t - step*pt;

            } else {

               // assign intermediate parameters
               amrex::ParticleReal const B = refpart.beta_gamma() /rc;
               amrex::ParticleReal const theta = slice_ds / rc;

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

            }

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

            // length of the current slice
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // find beta*gamma^2, beta
            amrex::ParticleReal const betgam2 = amrex::Math::powi<2>(refpart.pt) - 1_prt;
            amrex::ParticleReal const bet = refpart.beta();
            amrex::ParticleReal const ibetgam2 = 1_prt / betgam2;

            // host view, obtained at the point of use
            amrex::ParticleReal const * k_normal = m_coefficients_data.host_view().first;

            // find magnetic field and curvature; normalize units to MAD-X convention if needed
            // magnetic rigidity of the reference particle
            // note: read from the reference particle and not from the cached m_brho, because
            //       compute_constants() does not run on the transport-map paths (map_trace,
            //       envelope tracking, lattice.transfer_map)
            amrex::ParticleReal const brho = refpart.rigidity_Tm();

            amrex::ParticleReal const curvature = m_unit == 1 ? k_normal[0] / brho : k_normal[0];
            amrex::ParticleReal const rc = curvature == 0.0 ? 0.0_prt : 1.0_prt/curvature;
            amrex::ParticleReal const kn = m_unit == 1 ? k_normal[1] / brho : k_normal[1];
            amrex::ParticleReal const b2rc2 = amrex::Math::powi<2>(bet) * amrex::Math::powi<2>(rc);

            // update horizontal and longitudinal phase space variables
            amrex::ParticleReal const gx = kn + amrex::Math::powi<-2>(rc);
            amrex::ParticleReal const omega_x = std::sqrt(std::abs(gx));

            // update vertical phase space variables
            amrex::ParticleReal const gy = -kn;
            amrex::ParticleReal const omega_y = std::sqrt(std::abs(gy));

            // trigonometry
            auto const [sinx, cosx] = amrex::Math::sincos(omega_x * slice_ds);
            amrex::ParticleReal const sinhx = std::sinh(omega_x * slice_ds);
            amrex::ParticleReal const coshx = std::cosh(omega_x * slice_ds);
            auto const [siny, cosy] = amrex::Math::sincos(omega_y * slice_ds);
            amrex::ParticleReal const sinhy = std::sinh(omega_y * slice_ds);
            amrex::ParticleReal const coshy = std::cosh(omega_y * slice_ds);

            amrex::ParticleReal const igbrc = 1_prt / (gx * bet * rc);
            amrex::ParticleReal const iobrc = 1_prt / (omega_x * bet * rc);
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

        /** Prepare the execution snapshot for one element slice
         *
         * @param refpart reference particle, already advanced for this slice
         * @return a core holding an execution-space view of the coefficients
         */
        [[nodiscard]] ExactCFbendCore make_particle_pusher (RefPart const & refpart) const
        {
            ExactCFbendCore core = *this;
            // compute_constants() derives the curvature from the first coefficient, so it
            // needs host storage; the device view is installed afterwards.
            core.m_coefficients = m_coefficients_data.host_view();
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
        {            m_coefficients_data.commit(std::move(first), std::move(second), "ExactCFbend: normal and skew coefficient");
        }

        /** Release device storage; host configuration survives for reuse */
        void finalize ()
        {
            m_coefficients_data.release_gpu();
        }

        //! the coefficient arrays, owned by this element
        mixin::DeviceBackedArrays m_coefficients_data;
    };

    /** ExactCFbend: physics plus common host-only metadata */
    using ExactCFbend = Element<ExactCFbendPhysics>;

} // namespace impactx

IMPACTX_PUSH_EXTERN_TEMPLATE(impactx::elements::ExactCFbendPhysics)

#endif // IMPACTX_EXACTCFBEND_H
