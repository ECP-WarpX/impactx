/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Eric G. Stern, Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_POLYGONAPERTURE_H
#define IMPACTX_POLYGONAPERTURE_H

#include "particles/ImpactXParticleContainer.H"
#include "Element.H"
#include "mixin/alignment.H"
#include "mixin/beamoptic.H"
#include "mixin/devicebackedarrays.H"
#include "mixin/lineartransport.H"
#include "mixin/spintransport.H"
#include "mixin/thin.H"

#include <ablastr/constant.H>

#include <AMReX_Enum.H>
#include <AMReX_Extension.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>
#include <AMReX_SIMD.H>
#include <AMReX_TrackedVector.H>

#include <cmath>
#include <map>
#include <memory>
#include <stdexcept>
#include <vector>

namespace impactx::elements
{
    /** Polygon vertices for aperture boundary definition.
     *
     * Host data persists across simulations; device data is uploaded
     * on-demand "lazily" on first use and released when AMReX shuts down.
     */
    struct PolygonVertices
    {
        amrex::Gpu::TrackedVector<amrex::ParticleReal> x;
        amrex::Gpu::TrackedVector<amrex::ParticleReal> y;
    };

    /** The prepared execution snapshot of a polygon aperture.
     *
     * Holds the particle-facing configuration and a bounded, non-owning view of the
     * polygon vertices, resolved to the execution space of one element slice.
     */
    struct PolygonApertureCore
    : public mixin::ParticlePush<PolygonApertureCore>,
      public mixin::Thin,
      public mixin::Alignment,
      public mixin::SpinTransport
      // public amrex::simd::Vectorized<amrex::simd::native_simd_size_particlereal>
    {
        static constexpr auto type = "PolygonAperture";
        using PType = ImpactXParticleContainer::ParticleType;

        AMREX_ENUM_IN_CLASS(Action,
            transmit,  //! transmit particles inside the aperture
            absorb,    //! absorb particles inside the aperture
        );

        static std::string
        action_name (Action const & action);

        static constexpr amrex::ParticleReal DEFAULT_min_radius2 = 0;
        static constexpr amrex::ParticleReal DEFAULT_repeat_x = 0;
        static constexpr amrex::ParticleReal DEFAULT_repeat_y = 0;
        static constexpr bool DEFAULT_shift_odd_x = false;
        static constexpr Action DEFAULT_action = Action::transmit;

        /** The particle-facing configuration of a polygon aperture
         *
         * @param min_radius2 a minimum radius**2 below which all particles are transmitted
         * @param repeat_x horizontal period for repeated masking, optional (m)
         * @param repeat_y vertical period for repeated masking, optional (m)
         * @param shift_odd_x horizontal shift of every 2nd (odd) vertical period by repeat_x / 2
         * @param action specify action of domain (transmit/absorb)
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         */
        PolygonApertureCore (
            amrex::ParticleReal min_radius2,
            amrex::ParticleReal repeat_x,
            amrex::ParticleReal repeat_y,
            bool shift_odd_x,
            Action action,
            amrex::ParticleReal dx,
            amrex::ParticleReal dy,
            amrex::ParticleReal rotation_degree
        )
          : Alignment(dx, dy, rotation_degree),
            m_action(action), m_repeat_x(repeat_x), m_repeat_y(repeat_y),
            m_shift_odd_x(shift_odd_x), m_min_radius2(min_radius2)
        {
            using namespace amrex::literals; // for _rt and _prt

            if (m_repeat_y == 0_prt && m_shift_odd_x) {
                throw std::invalid_argument("PolygonAperture: shift_odd_x can only be used if repeat_y is set, too!");
            }
        }

        /** Reverse the element (no-op: aperture is not invertible for absorbed particles) */
        void reverse () { /* no-op */ }

        /** Push all particles on a tile */
        using mixin::ParticlePush<PolygonApertureCore>::operator();


        /** This is an aperture functor, so that a variable of this type can be used like an
         *  aperture function.
         *
         * @param x particle position in x
         * @param y particle position in y
         * @param t particle position in t (unused)
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t (unused)
         * @param idcpu particle global index
         * @param refpart reference particle (unused)
         */
        template<typename T_Real=amrex::ParticleReal, typename T_IdCpu=uint64_t>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
            T_Real const & AMREX_RESTRICT x,
            T_Real const & AMREX_RESTRICT y,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT t,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT px,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT py,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT pt,
            T_IdCpu & AMREX_RESTRICT idcpu,
            [[maybe_unused]] RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            // a complex type with two T_Real
            using Complex = amrex::GpuComplex<T_Real>;

            using namespace amrex::literals; // for _rt and _prt
            namespace stdx = amrex::simd::stdx;
            using amrex::Math::powi;
            using std::fmod;
            using std::abs;
            using amrex::arg;
            using VBool = decltype(x>0_prt);

            // define intermediate variables
            amrex::ParticleReal const dx = m_repeat_x * 0.5_prt;
            amrex::ParticleReal const dy = m_repeat_y * 0.5_prt;

            // shift every 2nd row by half of m_repeat_x
            T_Real sx = x;
            T_Real const sy = y;
            VBool const shift_x = m_shift_odd_x == false ? VBool{false} : fmod(floor((y+dy)/m_repeat_y), 2) != T_Real{0};
#ifdef AMREX_USE_SIMD
            amrex::simd::stdx::where(shift_x, sx) += dx;
#else
            sx += shift_x ? dx : 0_prt;
#endif

            // if the aperture is periodic, shift sx,sy coordinates to the fundamental domain
            T_Real const u = (m_repeat_x == 0_prt) ? sx : (fmod(abs(sx)+dx, m_repeat_x)-dx);
            T_Real const v = (m_repeat_y == 0_prt) ? sy : (fmod(abs(sy)+dy, m_repeat_y)-dy);

            // compare against the aperture boundary
            VBool inside_aperture;

            // calculate whether point is inside polygon or not
            T_Real const r2 = powi<2>(u) + powi<2>(v);
            // are we inside the minimum radius?
            if (r2 < m_min_radius2) {
                inside_aperture = VBool(true); // yes!
            } else{

                amrex::ParticleReal const * vertices_x = m_vertices.first;
                amrex::ParticleReal const * vertices_y = m_vertices.second;

                // no, we have to do the polygon calculation

                /*
                 * The algorithms works as follows:
                 *
                 * You have a polygon (in 2D) described by points V_1, ..., V_n arranged
                 * counter-clockwise. You have a particle with coordinates w = (x,y).
                 * Get the difference vector from the particle to each of
                 * the vertices in turn expressed as a complex number z_k, k=1,n. Now the angle
                 * between the particle and vertex i and vertex j is
                 *
                 * theta_{ij} = arg(\bar{z_i} z_j)
                 * This is easily demonstrated by expressing the complex numbers in exponential
                 * polar notation and also the angle is positive in the counter-clockwise direction.
                 *
                 * Take the sum of all the angles. The angle sum for particle inside the polygon will
                 * be 2 pi. The sum of angles of particle outside of the polygon will be 0.
                */
                T_Real thetasum = 0.0;
                Complex v1conj, v2;

                for (int i = 0; i < m_vertices.size - 1; ++i) {
                    v1conj = Complex(u - vertices_x[i], -(v-vertices_y[i]));
                    v2 = Complex(u - vertices_x[i+1], v-vertices_y[i+1]);
                    thetasum += arg(v2 * v1conj);
                }
                inside_aperture = abs(thetasum/(2*ablastr::constant::math::pi)) > 1.0e-12_prt;
            }

            // For transmit (default): invalidate if outside aperture
            // For absorb: invalidate if inside aperture
            auto const invalid_mask = m_action == Action::transmit ? !inside_aperture : inside_aperture;
            amrex::ParticleIDWrapper<T_IdCpu>{idcpu}.make_invalid(invalid_mask);
        }

        /** This pushes the reference particle. */
        using Thin::operator();

        /** This operator pushes the particle spin.
         *
         * An aperture is field-free, so the spin is unaffected.
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
            T_Real const & AMREX_RESTRICT x,
            T_Real const & AMREX_RESTRICT y,
            T_Real const & AMREX_RESTRICT t,
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

        /** This function returns the linear transport map.
         *
         * A polygon aperture (potentially) removes particles; it does not
         * deflect or drift the beam, so the linear transport map is the
         * identity.
         *
         * @param[in] refpart reference particle (unused)
         * @returns 6x6 transport matrix
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        Map6x6
        transport_map ([[maybe_unused]] RefPart const & AMREX_RESTRICT refpart) const
        {
            // apply the transverse rotation (roll) alignment error
            return rotate_aligned_map(Map6x6::Identity());
        }

        Action m_action; //! action type (transmit, absorb)
        amrex::ParticleReal m_repeat_x; //! horizontal period for repeated masking
        amrex::ParticleReal m_repeat_y; //! vertical period for repeated masking
        bool m_shift_odd_x; //! for hexagonal/triangular mask patterns: horizontal shift of every 2nd (odd) vertical period by m_repeat_x / 2. Use Alignment offsets to move whole mask as needed.

        /** Vertices resolved to this snapshot's execution space; non-owning */
        mixin::DualArrayView m_vertices{};

        amrex::ParticleReal m_min_radius2; // minimum radius in which all pass squared


    };

    /** A polygon aperture: the canonical element that owns its vertices.
     *
     * Particle pushes run against @see PolygonApertureCore.
     */
    struct PolygonAperturePhysics
    : public mixin::BeamOptic<PolygonAperturePhysics>,
      public mixin::LinearTransport<PolygonAperturePhysics>,
      public PolygonApertureCore
    {
        /** A thin collimator element that applies a transverse aperture boundary defined
         *  by points of a polygon vertices.
         *  Particles outside the boundary are considered lost.
         *
         * @param action specify action of domain (transmit/absorb)
         * @param vertices_x array of coordinates of the horizontal vertex positions (m)
         * @param vertices_y array of coordinates of the vertical vertex positions (m)
         * @param min_radius2 a minimum radius**2 below which all particles are transmitted
         * @param repeat_x horizontal period for repeated masking, optional (m)
         * @param repeat_y vertical period for repeated masking, optional (m)
         * @param shift_odd_x for hexagonal/triangular mask patterns: horizontal shift of every 2nd (odd) vertical period by repeat_x / 2. Use alignment offsets dx,dy to move whole mask as needed.
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         */
        PolygonAperturePhysics (
            std::vector< amrex::ParticleReal> vertices_x,
            std::vector< amrex::ParticleReal> vertices_y,
            amrex::ParticleReal min_radius2 = DEFAULT_min_radius2,
            amrex::ParticleReal repeat_x = DEFAULT_repeat_x,
            amrex::ParticleReal repeat_y = DEFAULT_repeat_y,
            bool shift_odd_x = DEFAULT_shift_odd_x,
            Action action = DEFAULT_action,
            amrex::ParticleReal dx = DEFAULT_dx,
            amrex::ParticleReal dy = DEFAULT_dy,
            amrex::ParticleReal rotation_degree = DEFAULT_rotation_degree
        )
          : PolygonApertureCore(min_radius2, repeat_x, repeat_y, shift_odd_x, action, dx, dy, rotation_degree),
            m_vertices_data(std::move(vertices_x), std::move(vertices_y), "PolygonAperture: horizontal and vertical vertices")
        {
            require_closed_outline(m_vertices_data.first_host(), m_vertices_data.second_host());
        }

        /** Push all particles */
        using mixin::BeamOptic<PolygonAperturePhysics>::operator();

        /** This pushes the reference particle. */
        using mixin::Thin::operator();

        /** This pushes the covariance matrix. */
        using mixin::LinearTransport<PolygonAperturePhysics>::operator();

        /** Prepare the execution snapshot for one element slice
         *
         * @param refpart reference particle, already advanced for this slice
         * @return a core holding an execution-space view of the vertices
         */
        [[nodiscard]] PolygonApertureCore make_particle_pusher (RefPart const & refpart) const
        {
            PolygonApertureCore core = *this;
            core.Alignment::compute_constants(refpart);
            core.m_vertices = m_vertices_data.device_view();
            return core;
        }

        /** Horizontal vertex coordinates */
        [[nodiscard]] std::vector<amrex::ParticleReal> vertices_x () const
        {
            return m_vertices_data.first_host();
        }

        /** Vertical vertex coordinates */
        [[nodiscard]] std::vector<amrex::ParticleReal> vertices_y () const
        {
            return m_vertices_data.second_host();
        }

        /** Replace both arrays at once
         *
         * The two are validated together and committed as one step, so a rejected update
         * leaves the element as it was. They are always the same length, which is why
         * this takes both rather than offering to change one of them alone.
         *
         * @param first the new horizontal values
         * @param second the new vertical values
         */
        void set_vertices (
            std::vector<amrex::ParticleReal> first,
            std::vector<amrex::ParticleReal> second
        )
        {
            require_closed_outline(first, second);
            m_vertices_data.commit(std::move(first), std::move(second), "PolygonAperture: horizontal and vertical vertices");
        }

        /** Reject vertices that do not describe a closed outline
         *
         * Checked before anything is read positionally, so a mismatched or empty pair is
         * reported rather than indexed into.
         *
         * @param vx the horizontal values
         * @param vy the vertical values
         */
        template<typename T_Values>
        static void require_closed_outline (T_Values const & vx, T_Values const & vy)
        {
            if (vx.size() != vy.size()) {
                throw std::invalid_argument(
                    "PolygonAperture: horizontal and vertical vertices must have the same "
                    "length, got " + std::to_string(vx.size()) + " and " +
                    std::to_string(vy.size()) + "!");
            }
            if (vx.empty()) {
                throw std::invalid_argument("PolygonAperture: at least one vertex is required!");
            }
            if ((vx.front() != vx.back()) || (vy.front() != vy.back())) {
                throw std::invalid_argument("PolygonAperture: first and last vertex must be exactly equal!");
            }
        }

        /** Release device storage; host configuration survives for reuse */
        void finalize ()
        {
            m_vertices_data.release_gpu();
        }

        //! the vertex arrays, owned by this element
        mixin::DeviceBackedArrays m_vertices_data;
    };

    /** PolygonAperture: physics plus common host-only metadata */
    using PolygonAperture = Element<PolygonAperturePhysics>;

} // namespace impactx

IMPACTX_PUSH_EXTERN_TEMPLATE(impactx::elements::PolygonAperturePhysics)

#endif // IMPACTX_POLYGONAPERTURE_H
