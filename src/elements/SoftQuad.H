/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_SOFTQUAD_H
#define IMPACTX_SOFTQUAD_H

#include "particles/ImpactXParticleContainer.H"
#include "Element.H"
#include "particles/integrators/Integrators.H"
#include "mixin/alignment.H"
#include "mixin/beamoptic.H"
#include "mixin/devicebackedarrays.H"
#include "mixin/lineartransport.H"
#include "mixin/pipeaperture.H"
#include "mixin/thick.H"
#include "mixin/spintransport.H"

#include <ablastr/constant.H>

#include <AMReX.H>
#include <AMReX_Extension.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>
#include <AMReX_SIMD.H>
#include <AMReX_SmallMatrix.H>
#include <AMReX_TrackedVector.H>

#include <cmath>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace impactx::elements
{
   /** Default Fourier coefficients
    *
    * Fourier coefficients for the on-axis quadrupole field gradient
    * k = 1/(Brho)*dBy/dx as a function of longitudinal position z.
    * The default on-axis gradient is described by the function:
    * k(z) = gscale/2*(tanh((z-zi)/g) - tanh((z-zf)/g)),
    * where zf = -zi = L/2, L is the magnetic length of the
    * quadrupole, and g is the full gap parameter.  Here the
    * parameters are chosen so that g/L = 1/10.  This is
    * based on a similar model used for soft-edge elements
    * in MaryLie and ML/I, eg.
    * http://www.physics.umd.edu/dsat/docs/MaryLieMan.pdf
    *
    * This is equivalent to a single-parameter Enge model for
    * the fields in the fringe region.
    *
    * Usually superseded by user-provided input values.
    */
    struct Quad_field_data
    {
       amrex::Vector<amrex::ParticleReal> default_cos_coef = {
             0.834166514794446,
             0.598104328994702,
             0.141852844428785,
            -0.118211272182381,
            -9.056149864743113E-002,
             1.803476331179615E-002,
             4.464887700797893E-002,
             7.364410636252136E-003,
            -1.697541023436736E-002,
            -9.012679515542771E-003,
             4.367667630047725E-003,
             5.444030542119803E-003,
            -5.889959910931886E-005,
            -2.409098101409192E-003,
            -7.962712154165590E-004,
             7.855814707106538E-004,
             6.174930463182168E-004,
            -1.340154094301854E-004,
            -3.167213724698439E-004,
            -4.925292460592617E-005,
             1.221580597451921E-004,
             6.331025910961789E-005,
            -3.202416719002774E-005,
            -3.872103053895529E-005,
             8.212882937116278E-007
            };

       amrex::Vector<amrex::ParticleReal> default_sin_coef = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0
            };
    };

    /** Fourier coefficients for on-axis quadrupole field gradient.
     *
     * Host data persists across simulations; device data is uploaded
     * on-demand "lazily" on first use and released when AMReX shuts down.
     */
    struct QuadrupoleFourierCoefficients
    {
        amrex::Gpu::TrackedVector<amrex::ParticleReal> cos;
        amrex::Gpu::TrackedVector<amrex::ParticleReal> sin;
    };

    /** The prepared, resource-free execution snapshot of a soft-edge quadrupole.
     *
     * This is what a particle-push kernel is handed. It owns nothing: the coefficient
     * arrays reach it as a bounded, non-owning @see mixin::DualArrayView, resolved to
     * host or device storage by @see SoftQuadrupolePhysics::make_particle_pusher for the
     * lifetime of one element slice.
     *
     * It deliberately derives @see mixin::ParticlePush rather than @see mixin::BeamOptic,
     * so it can be pushed over a tile but cannot start a lattice push or advance the
     * reference particle. Those need the canonical element, which owns the arrays.
     */
    struct SoftQuadrupoleCore
    : public mixin::ParticlePush<SoftQuadrupoleCore>,
      public mixin::Thick,
      public mixin::Alignment,
      public mixin::PipeAperture,
      public mixin::SpinTransport,
      public amrex::simd::Vectorized<amrex::simd::native_simd_size_particlereal>
    {
        static constexpr auto type = "SoftQuadrupole";
        using PType = ImpactXParticleContainer::ParticleType;

        static constexpr int DEFAULT_mapsteps = 10;

        /** The particle-facing configuration of a soft-edge quadrupole
         *
         * @param ds Segment length in m
         * @param gscale Scaling factor for on-axis field gradient Bz in 1/m^2
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param aperture_x horizontal half-aperture in m
         * @param aperture_y vertical half-aperture in m
         * @param mapsteps number of integration steps per slice
         * @param nslice number of slices used for the application of space charge
         */
        SoftQuadrupoleCore (
            amrex::ParticleReal ds,
            amrex::ParticleReal gscale,
            amrex::ParticleReal dx,
            amrex::ParticleReal dy,
            amrex::ParticleReal rotation_degree,
            amrex::ParticleReal aperture_x,
            amrex::ParticleReal aperture_y,
            int mapsteps,
            int nslice
        )
          : Thick(ds, nslice),
            Alignment(dx, dy, rotation_degree),
            PipeAperture(aperture_x, aperture_y),
            m_gscale(gscale), m_mapsteps(mapsteps)
        {
        }

        /** Reverse the element (negate ds) */
        void reverse () {
            // Reversing ds traverses the Fourier profile in the opposite z direction,
            // so the odd sine terms change sign implicitly via sin(-kz) = -sin(kz).
            Thick::reverse();
        }

        /** Push all particles on a tile */
        using mixin::ParticlePush<SoftQuadrupoleCore>::operator();

        /** This is a soft-edge quadrupole functor, so that a variable of this type can be used
         *  like a soft-edge quadrupole function.
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
            [[maybe_unused]] RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // get the linear map
            amrex::SmallMatrix<amrex::ParticleReal, 6, 6, amrex::Order::F, 1> const R = m_map;

            // symplectic linear map for a quadrupole is computed using the
            // Hamiltonian formalism as described in:
            // https://uspas.fnal.gov/materials/09UNM/ComputationalMethods.pdf .
            // R denotes the transfer matrix in the basis (x,px,y,py,t,pt),
            // so that, e.g., R(3,4) = dyf/dpyi.
            amrex::SmallVector<T_Real, 6, 1> const v{x, px, y, py, t, pt};

            // push particles using the linear map
            auto const out = R * v;

            // assign updated values
            x  = out[1];
            px = out[2];
            y  = out[3];
            py = out[4];
            t  = out[5];
            pt = out[6];
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
         * @param refpart reference particle
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
            T_IdCpu const & AMREX_RESTRICT idcpu,
            RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // initialize the three components of the axis-angle vector
            T_Real lambdax = 0_prt;
            T_Real lambday = 0_prt;
            T_Real lambdaz = 0_prt;

            // store the phase space variables in vector form
            amrex::SmallVector<T_Real, 6, 1> const v{x, px, y, py, t, pt};

            // get the spin-orbit coupling matrix
            amrex::SmallMatrix<amrex::ParticleReal, 3, 6, amrex::Order::F, 1> const A = m_spin_coupling;

            // use phase space variables to obtain the angle-axis generator of spin rotation
            auto const out = A * v;

            // update the angle-axis generator
            lambdax = out[1];
            lambday = out[2];
            lambdaz = out[3];

            // push the spin vector using the generator just determined
            rotate_spin(lambdax,lambday,lambdaz,sx,sy,sz);

            // phase space push
            (*this)(x, y, t, px, py, pt, idcpu, refpart);
        }

        /** This evaluates the on-axis magnetic field Bz at a fixed location
         *  z, together with certain required integrals and derivatives.
         *  The field returned is normalized to a peak value of 1.
         *
         * @param zeval Longitudinal on-axis location in m
         */
        std::tuple<amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        Quad_Bfield (amrex::ParticleReal const zeval) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // The coefficients were already resolved to the right execution space when
            // this snapshot was prepared, so there is no host/device branch here and only
            // one copy of the field expansion below.
            amrex::ParticleReal const * cos_data = m_coefficients.first;
            amrex::ParticleReal const * sin_data = m_coefficients.second;
            int const ncoef = m_coefficients.size;

            // specify constants
            using ablastr::constant::math::pi;
            amrex::ParticleReal const zlen = std::abs(m_ds);
            amrex::ParticleReal const zmid = zlen * 0.5_prt;

            // compute on-axis magnetic field (z is relative to quadrupole midpoint)
            amrex::ParticleReal bfield = 0.0;
            amrex::ParticleReal bfieldp = 0.0;
            amrex::ParticleReal bfieldint = 0.0;
            amrex::ParticleReal const z = zeval - zmid;

            if (std::abs(z) <= zmid)
            {
               bfield = 0.5_prt*cos_data[0];
               bfieldint = z*bfield;
               for (int j=1; j < ncoef; ++j)
               {
                 bfield = bfield + cos_data[j] *std::cos(j * 2 * pi * z / zlen) +
                         sin_data[j] *std::sin(j * 2 * pi * z / zlen);
                 bfieldp = bfieldp - j * 2 * pi * cos_data[j] *std::sin(j * 2 * pi * z / zlen) / zlen +
                           j * 2 * pi * sin_data[j] *std::cos(j * 2 * pi * z / zlen) / zlen;
                 bfieldint = bfieldint + zlen * cos_data[j] *std::sin(j * 2 * pi * z / zlen) / (j * 2 * pi) -
                             zlen * sin_data[j] *std::cos(j * 2 * pi * z / zlen) / (j * 2 * pi);
               }
            }
            return std::make_tuple(bfield, bfieldp, bfieldint);
        }

        /** This pushes the reference particle and the linear map matrix
         *  elements for a quadrupole through the symplectic map associated
         *  with H_1 in the Hamiltonian splitting H = H_1 + H_2.
         *
         * @param tau Map step size in m
         * @param[in,out] refpart reference particle
         * @param[in,out] zeval Longitudinal on-axis location in m
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void map1 (amrex::ParticleReal const tau,
                   RefPart & refpart,
                   [[maybe_unused]] amrex::ParticleReal & zeval) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using amrex::Math::powi;

            // push the reference particle
            amrex::ParticleReal const t = refpart.t;
            amrex::ParticleReal const pt = refpart.pt;
            amrex::ParticleReal const z = zeval;

            if (pt < -1.0_prt) {
                refpart.t = t + tau/std::sqrt(1.0_prt - powi<-2>(pt));
                refpart.pt = pt;
            }
            else {
                refpart.t = t;
                refpart.pt = pt;
            }

            zeval = z + tau;

            // push the linear map equations
            amrex::SmallMatrix<amrex::ParticleReal, 6, 6, amrex::Order::F, 1> const R = m_map;
            amrex::ParticleReal const betgam = refpart.beta_gamma();

            m_map(1,1) = R(1,1) + tau*R(2,1);
            m_map(1,2) = R(1,2) + tau*R(2,2);
            m_map(1,3) = R(1,3) + tau*R(2,3);
            m_map(1,4) = R(1,4) + tau*R(2,4);

            m_map(3,1) = R(3,1) + tau*R(4,1);
            m_map(3,2) = R(3,2) + tau*R(4,2);
            m_map(3,3) = R(3,3) + tau*R(4,3);
            m_map(3,4) = R(3,4) + tau*R(4,4);

            m_map(5,5) = R(5,5) + tau*R(6,5) / powi<2>(betgam);
            m_map(5,6) = R(5,6) + tau*R(6,6) / powi<2>(betgam);
        }

        /** This pushes the reference particle and the linear map matrix
         *  elements for a quadrupole through the symplectic map associated
         *  with H_2 in the Hamiltonian splitting H = H_1 + H_2.
         *
         * @param tau Map step size in m
         * @param[in,out] refpart reference particle
         * @param[in,out] zeval Longitudinal on-axis location in m
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void map2 (amrex::ParticleReal const tau,
                   RefPart & refpart,
                   amrex::ParticleReal & zeval) const
        {
            using namespace amrex::literals; // for _rt and _prt

            amrex::ParticleReal const t = refpart.t;
            amrex::ParticleReal const pt = refpart.pt;

            // Define parameters and intermediate constants
            amrex::ParticleReal const G0 = m_gscale;

            // push the reference particle
            auto [bz, bzp, bzint] = Quad_Bfield(zeval);
            amrex::ignore_unused(bzp, bzint);

            refpart.t = t;
            refpart.pt = pt;

            // push the linear map equations
            amrex::SmallMatrix<amrex::ParticleReal, 6, 6, amrex::Order::F, 1> const R = m_map;
            amrex::ParticleReal const alpha = G0*bz;

            m_map(2,1) = R(2,1) - tau*alpha*R(1,1);
            m_map(2,2) = R(2,2) - tau*alpha*R(1,2);
            m_map(2,3) = R(2,3) - tau*alpha*R(1,3);
            m_map(2,4) = R(2,4) - tau*alpha*R(1,4);

            m_map(4,1) = R(4,1) + tau*alpha*R(3,1);
            m_map(4,2) = R(4,2) + tau*alpha*R(3,2);
            m_map(4,3) = R(4,3) + tau*alpha*R(3,3);
            m_map(4,4) = R(4,4) + tau*alpha*R(3,4);

            // BELOW: if spin is needed only:
            amrex::SmallMatrix<amrex::ParticleReal, 3, 6, amrex::Order::F, 1> const A = m_spin_coupling;
            amrex::SmallMatrix<amrex::ParticleReal, 3, 6, amrex::Order::F, 1> dA = {};
            amrex::ParticleReal const gamma = refpart.gamma();
            amrex::ParticleReal G = refpart.gyromagnetic_anomaly;

            // Update spin-orbit coupling matrix here
            dA(1,3) = -(1_prt + G*gamma) * tau * alpha;
            dA(2,1) = -(1_prt + G*gamma) * tau * alpha;

            // update the spin-orbit coupling matrix here
            m_spin_coupling = A + dA*R;

        }

        amrex::ParticleReal m_gscale; //! scaling factor for quad field gradient
        int m_mapsteps; //! number of map integration steps per slice
        /** Coefficients resolved to this snapshot's execution space.
         *
         * Non-owning and valid only for the launch this snapshot was prepared for.
         * @see SoftQuadrupolePhysics, which owns the arrays.
         */
        mixin::DualArrayView m_coefficients{};


        // Reference-trajectory linearization around the reference particle.
        // Computed during the reference-particle push and consumed during the
        // particle push.
        // mutable: written by the const reference push before the element is copied
        // into the particle-push functor.
        mutable amrex::SmallMatrix<amrex::ParticleReal, 6, 6, amrex::Order::F, 1> m_map{}; //! linearized map
        mutable amrex::SmallMatrix<amrex::ParticleReal, 3, 6, amrex::Order::F, 1> m_spin_coupling{}; //! linearized spin-orbit coupling matrix
    };

    /** A soft-edge quadrupole: the canonical element that owns its coefficients.
     *
     * Everything that runs on the host and needs the coefficient arrays lives here: the
     * reference-particle push, the transfer map, serialization and finalization. Particle
     * pushes run against @see SoftQuadrupoleCore, which this hands out already resolved to
     * the right execution space.
     */
    struct SoftQuadrupolePhysics
    : public mixin::BeamOptic<SoftQuadrupolePhysics>,
      public mixin::LinearTransport<SoftQuadrupolePhysics>,
      public SoftQuadrupoleCore
    {
        /** A soft-edge quadrupole
         *
         * @param ds Segment length in m
         * @param gscale Scaling factor for on-axis field gradient Bz in 1/m^2
         * @param cos_coef cosine coefficients in Fourier expansion of on-axis magnetic field Bz
         * @param sin_coef sine coefficients in Fourier expansion of on-axis magnetic field Bz
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param aperture_x horizontal half-aperture in m
         * @param aperture_y vertical half-aperture in m
         * @param mapsteps number of integration steps per slice used for
         *        map and reference particle push in applied fields
         * @param nslice number of slices used for the application of space charge
         */
        SoftQuadrupolePhysics (
            amrex::ParticleReal ds,
            amrex::ParticleReal gscale,
            std::vector<amrex::ParticleReal> cos_coef,
            std::vector<amrex::ParticleReal> sin_coef,
            amrex::ParticleReal dx = DEFAULT_dx,
            amrex::ParticleReal dy = DEFAULT_dy,
            amrex::ParticleReal rotation_degree = DEFAULT_rotation_degree,
            amrex::ParticleReal aperture_x = DEFAULT_aperture_x,
            amrex::ParticleReal aperture_y = DEFAULT_aperture_y,
            int mapsteps = DEFAULT_mapsteps,
            int nslice = DEFAULT_nslice
        )
          : SoftQuadrupoleCore(ds, gscale, dx, dy, rotation_degree, aperture_x, aperture_y, mapsteps, nslice),
            m_coefficients_data(std::move(cos_coef), std::move(sin_coef), "SoftQuadrupole: cos and sin coefficient", 1)
        {
            if (m_mapsteps <= 0)
                throw std::invalid_argument("SoftQuadrupole: mapsteps must be > 0!");
        }

        /** Push all particles */
        using mixin::BeamOptic<SoftQuadrupolePhysics>::operator();

        /** This pushes the reference particle.
         *
         * @param[in,out] refpart reference particle
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        void operator() (RefPart & AMREX_RESTRICT refpart) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using amrex::Math::powi;

            // Run the shared field expansion against a fresh *host* view. The view is
            // local to this call, so the canonical element never caches a pointer into
            // storage that a later mutation could replace.
            SoftQuadrupoleCore core = *this;
            core.m_coefficients = m_coefficients_data.host_view();

            // assign input reference particle values
            amrex::ParticleReal const x = refpart.x;
            amrex::ParticleReal const px = refpart.px;
            amrex::ParticleReal const y = refpart.y;
            amrex::ParticleReal const py = refpart.py;
            amrex::ParticleReal const z = refpart.z;
            amrex::ParticleReal const pz = refpart.pz;
            amrex::ParticleReal const pt = refpart.pt;
            amrex::ParticleReal const s = refpart.s;
            amrex::ParticleReal const sedge = refpart.sedge;

            // initialize linear map (deviation) values
            core.m_map = decltype(m_map)::Identity();

            // initialize the spin-orbit coupling matrix
            core.m_spin_coupling = {};

            // length of the current slice
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // compute initial value of beta*gamma
            amrex::ParticleReal const bgi = std::sqrt(powi<2>(pt) - 1.0_prt);

            // call integrator to advance (t,pt)
            amrex::ParticleReal const zin = s - sedge;
            amrex::ParticleReal const zout = zin + slice_ds;
            int const nsteps = m_mapsteps;

            integrators::symp2_integrate(refpart,zin,zout,nsteps,core);

            // keep the linearization the integrator just produced
            m_map = core.m_map;
            m_spin_coupling = core.m_spin_coupling;
            amrex::ParticleReal const ptf = refpart.pt;

            /*
            // print computed linear map:
               for(int i=1; i<7; ++i){
                 for(int j=1; j<7; ++j){
                    amrex::PrintToFile("QuadMap.txt") << i << " " <<
                    j << " " << m_map(i,j) << "\n";
                 }
               }
            //
            */

            // advance position (x,y,z)
            refpart.x = x + slice_ds*px/bgi;
            refpart.y = y + slice_ds*py/bgi;
            refpart.z = z + slice_ds*pz/bgi;

            // compute final value of beta*gamma
            amrex::ParticleReal const bgf = std::sqrt(powi<2>(ptf) - 1.0_prt);

            // advance momentum (px,py,pz)
            refpart.px = px*bgf/bgi;
            refpart.py = py*bgf/bgi;
            refpart.pz = pz*bgf/bgi;

            // advance integrated path length
            refpart.s = s + slice_ds;
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

            Map6x6 R = Map6x6::Identity();
            R = m_map;

            // apply the transverse rotation (roll) alignment error
            return rotate_aligned_map(R);
        }
        /** Prepare the execution snapshot for one element slice
         *
         * @param refpart reference particle, already advanced for this slice
         * @return a resource-free core holding an execution-space view of the coefficients
         */
        [[nodiscard]] SoftQuadrupoleCore make_particle_pusher (RefPart const & refpart) const
        {
            SoftQuadrupoleCore core = *this;
            core.Alignment::compute_constants(refpart);
            core.m_coefficients = m_coefficients_data.device_view();
            return core;
        }

        /** Cosine coefficients, e.g. for serialization */
        [[nodiscard]] std::vector<amrex::ParticleReal> cos_coefficients () const
        {
            return m_coefficients_data.first_host();
        }

        /** Sine coefficients, e.g. for serialization */
        [[nodiscard]] std::vector<amrex::ParticleReal> sin_coefficients () const
        {
            return m_coefficients_data.second_host();
        }

        /** Replace both arrays at once
         *
         * The two are validated together and committed as one step, so a rejected update
         * leaves the element as it was. They are always the same length, which is why
         * this takes both rather than offering to change one of them alone.
         *
         * @param first the new cosine values
         * @param second the new sine values
         */
        void set_coefficients (
            std::vector<amrex::ParticleReal> first,
            std::vector<amrex::ParticleReal> second
        )
        {
            m_coefficients_data.commit(std::move(first), std::move(second), "SoftQuadrupole: cos and sin coefficient", 1);
        }

        /** Release device storage; host configuration survives for reuse */
        void finalize ()
        {
            m_coefficients_data.release_gpu();
        }

        //! the coefficient arrays, owned by this element
        mixin::DeviceBackedArrays m_coefficients_data;
    };

    /** SoftQuadrupole: physics plus common host-only metadata */
    using SoftQuadrupole = Element<SoftQuadrupolePhysics>;

} // namespace impactx

IMPACTX_PUSH_EXTERN_TEMPLATE(impactx::elements::SoftQuadrupolePhysics)

#endif // IMPACTX_SOFTQUAD_H
