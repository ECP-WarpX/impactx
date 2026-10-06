/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_PARTICLE_CONTAINER_H
#define IMPACTX_PARTICLE_CONTAINER_H

#include "ReferenceParticle.H"
#include "initialization/AmrCoreData_fwd.H"

#include <AMReX_AmrCoreFwd.H>
#include <AMReX_BaseFwd.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParIter.H>
#include <AMReX_Particles.H>
#include <AMReX_ParticleTile.H>
#include <AMReX_REAL.H>

#include <AMReX_IntVect.H>
#include <AMReX_Vector.H>

#include <list>
#include <optional>
#include <tuple>
#include <unordered_map>


namespace impactx
{
    enum class CoordSystem
    {
        s,  ///< fixed s as the independent variable
        t   ///< fixed t as the independent variable
    };

    /** This struct indexes the Real attributes
     *  stored in an SoA in ImpactXParticleContainer
     */
    struct RealSoA
    {
        enum
        {
            x,   ///< position in x [m] (at fixed s or t)
            y,   ///< position in y [m] (at fixed s or t)
            t,   ///< time-of-flight ct [m] (at fixed s)
            px,  ///< momentum in x, scaled by the magnitude of the reference momentum [unitless] (at fixed s or t)
            py,  ///< momentum in y, scaled by the magnitude of the reference momentum [unitless] (at fixed s or t)
            pt,  ///< energy deviation, scaled by speed of light * the magnitude of the reference momentum [unitless] (at fixed s)
            sx,  ///< spin vector x-component [unitless] (at fixed s or t)
            sy,  ///< spin vector y-component [unitless] (at fixed s or t)
            sz,  ///< spin vector z-component [unitless] (at fixed s or t)
            qm,  ///< charge to mass ratio, in q_e/m_e [q_e/eV]
            w,   ///< particle weight, number of real particles represented by this macroparticle [unitless]
            nattribs ///< the number of attributes above (always last)
        };

        // at fixed t, the third component represents the position z, the 6th component represents the momentum in z
        enum {
            z = t,   ///< position in z [m] (at fixed t)
            pz = pt  ///< momentum in z, scaled by the magnitude of the reference momentum [unitless] (at fixed t)
        };

        //! named labels for fixed s
        static constexpr auto names_s = { "position_x", "position_y", "position_t", "momentum_x", "momentum_y", "momentum_t", "spin_x", "spin_y", "spin_z", "qm", "weighting" };
        //! named labels for fixed t
        static constexpr auto names_t = { "position_x", "position_y", "position_z", "momentum_x", "momentum_y", "momentum_z", "spin_x", "spin_y", "spin_z", "qm", "weighting" };
        static_assert(names_s.size() == nattribs);
        static_assert(names_t.size() == nattribs);
    };

    /** This struct indexes the Integer attributes
     *  stored in an SoA in ImpactXParticleContainer
     */
    struct IntSoA
    {
        enum
        {
            nattribs ///< the number of attributes above (always last)
        };

        //! named labels for fixed s
        static constexpr std::initializer_list<char const *> names_s = {};
        //! named labels for fixed t
        static constexpr std::initializer_list<char const *> names_t = {};
        static_assert(names_s.size() == nattribs);
        static_assert(names_t.size() == nattribs);
    };

    /** AMReX iterator for particle boxes
     *
     * We subclass here to change the default threading strategy, which is
     * `static` in AMReX, to `dynamic` in ImpactX.
     */
    class ParIterSoA
        : public amrex::ParIterSoA<RealSoA::nattribs, IntSoA::nattribs, amrex::PolymorphicArenaAllocator>
    {
    public:
        using amrex::ParIterSoA<RealSoA::nattribs, IntSoA::nattribs, amrex::PolymorphicArenaAllocator>::ParIterSoA;

        ParIterSoA (ContainerType& pc, int level);

        ParIterSoA (ContainerType& pc, int level, amrex::MFItInfo& info);

        ContainerType * pc () { return m_pc; }
    };

    /** Const AMReX iterator for particle boxes - data is read only.
     *
     * We subclass here to change the default threading strategy, which is
     * `static` in AMReX, to `dynamic` in ImpactX.
     */
    class ParConstIterSoA
        : public amrex::ParConstIterSoA<RealSoA::nattribs, IntSoA::nattribs, amrex::PolymorphicArenaAllocator>
    {
    public:
        using amrex::ParConstIterSoA<RealSoA::nattribs, IntSoA::nattribs, amrex::PolymorphicArenaAllocator>::ParConstIterSoA;

        ParConstIterSoA (ContainerType& pc, int level);

        ParConstIterSoA (ContainerType& pc, int level, amrex::MFItInfo& info);

        ContainerType const * pc () const { return m_pc; }
    };

    /** Beam Particles in ImpactX
     *
     * This class stores particles, distributed over MPI ranks.
     */
    class ImpactXParticleContainer
        : public amrex::ParticleContainerPureSoA<RealSoA::nattribs, IntSoA::nattribs, amrex::PolymorphicArenaAllocator>
    {
    public:
        //! amrex iterator for particle boxes
        using iterator = impactx::ParIterSoA;

        //! amrex constant iterator for particle boxes (read-only)
        using const_iterator = impactx::ParConstIterSoA;

        //! Construct a new particle container
        ImpactXParticleContainer (initialization::AmrCoreData* amr_core);

        //! Destruct a particle container
        virtual ~ImpactXParticleContainer() = default;

        /** Prepares the container for use.
         *
         * This sets the tile size to a sensible default and resizes the
         * underlying data structures.
         */
        void prepare ();

        /** Empty the container and reset the reference particle
         *
         * Also clears the history of the beam moments, see reset_beam_moments_history()
         *
         * @param keep_mass do not reset the reference particle mass
         * @param keep_charge do not reset the reference particle charge
         */
        void clear (bool keep_mass=false, bool keep_charge=false);

        /** Add new particles to the container for fixed s.
         *
         * Either the total charge (bunch_charge) or the weight of each
         * particle (w) must be provided.
         *
         * Note: This can only be used *after* the initialization (grids) have
         *       been created, meaning after the call to AmrCore::InitFromScratch
         *       or AmrCore::InitFromCheckpoint has been made in the ImpactX
         *       class.
         *
         * @param x positions in x
         * @param y positions in y
         * @param t positions as time-of-flight in c*t
         * @param px momentum in x
         * @param py momentum in y
         * @param pt momentum in t
         * @param qm charge over mass in 1/eV
         * @param bunch_charge total charge within a bunch in C
         * @param w weight of each particle: how many real particles to represent
         * @param sx spin component in x
         * @param sy spin component in y
         * @param sz spin component in z
         */
        void
        AddNParticles (
            amrex::Gpu::DeviceVector<amrex::ParticleReal> const & x,
            amrex::Gpu::DeviceVector<amrex::ParticleReal> const & y,
            amrex::Gpu::DeviceVector<amrex::ParticleReal> const & t,
            amrex::Gpu::DeviceVector<amrex::ParticleReal> const & px,
            amrex::Gpu::DeviceVector<amrex::ParticleReal> const & py,
            amrex::Gpu::DeviceVector<amrex::ParticleReal> const & pt,
            amrex::ParticleReal qm,
            std::optional<amrex::ParticleReal> bunch_charge = std::nullopt,
            std::optional<amrex::Gpu::DeviceVector<amrex::ParticleReal>> w = std::nullopt,
            std::optional<amrex::Gpu::DeviceVector<amrex::ParticleReal>> sx = std::nullopt,
            std::optional<amrex::Gpu::DeviceVector<amrex::ParticleReal>> sy = std::nullopt,
            std::optional<amrex::Gpu::DeviceVector<amrex::ParticleReal>> sz = std::nullopt
        );

        /** Register storage for lost particles
         *
         * @param lost_pc particle container for lost particles
         */
        void
        SetLostParticleContainer (ImpactXParticleContainer * lost_pc);

        ImpactXParticleContainer *
        GetLostParticleContainer ();

        /** Set reference particle attributes
         *
         * @param refpart reference particle
         */
        void
        SetRefParticle (RefPart const & refpart);

        /** Get reference particle attributes
         *
         * @returns refpart
         */
        RefPart &
        GetRefParticle ();

        /** Get a copy to the reference particle attributes
         *
         * @returns refpart
         */
        RefPart const &
        GetRefParticle () const;

        /** Get particle shape
         */
        int
        GetParticleShape () const { return m_particle_shape.value(); }

        /** Set Particle Shape from amrex::ParmParse inputs
         *
         * Note: this can only be called once. All later calls are a logic error.
         * The reason for that is that subsequent calls would need to change
         * the guard size of all our MultiFabs, which is not implemented.
         */
        void SetParticleShape ();

        /** Set Particle Shape order
         *
         * Note: this can only be called once. All later calls are a logic error.
         * The reason for that is that subsequent calls would need to change
         * the guard size of all our MultiFabs, which is not implemented.
         *
         * @param order the order of the particle shape
         */
        void SetParticleShape (int order);

        /** Compute the min and max of the particle position in each dimension
         *
         * @returns x_min, y_min, z_min, x_max, y_max, z_max
         */
        std::tuple<
            amrex::ParticleReal, amrex::ParticleReal,
            amrex::ParticleReal, amrex::ParticleReal,
            amrex::ParticleReal, amrex::ParticleReal>
        MinAndMaxPositions ();

        /** Compute the mean and std of the particle position in each dimension
         *
         * @returns x_mean, x_std, y_mean, y_std, z_mean, z_std
         */
        std::tuple<
                amrex::ParticleReal, amrex::ParticleReal,
                amrex::ParticleReal, amrex::ParticleReal,
                amrex::ParticleReal, amrex::ParticleReal>
        MeanAndStdPositions ();

        /** Deposit the charge of the particles onto a grid
         *
         * This resets the values in rho to zero and then deposits the particle
         * charge. In MPI-parallel contexts, this also performs a communication
         * of boundary regions to sum neighboring contributions.
         *
         * @param[out] rho charge grid per level to deposit on
         * @param[in] ref_ratio mesh refinement ratios between levels
         */
        void
        DepositCharge (std::unordered_map<int, amrex::MultiFab> & rho,
                       amrex::Vector<amrex::IntVect> const & ref_ratio);

        /** Get the current coordinate system of particles in this container */
        CoordSystem
        GetCoordSystem () const;

        /** Set the current coordinate system of particles in this container
         *
         * This also names the SoA components after the coordinates they hold:
         * RealSoA::names_s at fixed s and RealSoA::names_t at fixed t.
         *
         * @param coord_system the new coordinate system the beam particles are in
         */
        void
        SetCoordSystem (CoordSystem coord_system);

        /** Set the bucket length used for particle boundaries in this container
         *
         * Note that the bucket length is only read and applied if ParmParse
         * "algo.particle_bc" is not set to "open" (default).
         *
         * @param bucket_length the length of the RF bucket in m
         */
        void SetBucketLength (amrex::ParticleReal bucket_length);

        /** Get the bucket length used for periodic particle BCs in this container
         *
         * @returns bucket length in m
         */
        amrex::ParticleReal
        GetBucketLength ();

        /** In situ calculate and store the beam moments for every simulation step */
        bool store_beam_moments = false;

        /** Calculate & record the beam moments at current s */
        void
        record_beam_moments ();

        /** Calculate beam moments at current s */
        std::unordered_map<std::string, amrex::ParticleReal>
        beam_moments ();

        /** Get the history of beam moments */
        std::list<std::unordered_map<std::string, amrex::ParticleReal> >
        beam_moments_history () { return m_beam_moments; }

        /** Reset the history of the beam moments */
        void
        reset_beam_moments_history () { m_beam_moments.clear(); }

      private:

        //! the reference particle for the beam in the particle container
        RefPart m_refpart;

        //! the particle shape
        std::optional<int> m_particle_shape;

        //! a non-owning reference to lost particles, i.e., due to apertures
        ImpactXParticleContainer* m_particles_lost = nullptr;

        //! the current coordinate system of particles in this container
        CoordSystem m_coordsystem = CoordSystem::s;

        //! history of the beam moments over s
        std::list<std::unordered_map<std::string, amrex::ParticleReal> > m_beam_moments;

        //! the RF bucket length for particles in the particle container
        //!   only active if algo.particle_bc is not set to "open" (default)
        amrex::ParticleReal m_bucket_length = 0.0;

    }; // ImpactXParticleContainer

} // namespace impactx

#endif // IMPACTX_PARTICLE_CONTAINER_H
