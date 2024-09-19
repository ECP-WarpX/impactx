/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
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
        static constexpr auto names_s = { "position_x", "position_y", "position_t", "momentum_x", "momentum_y", "momentum_t", "qm", "weighting" };
        //! named labels for fixed t
        static constexpr auto names_t = { "position_x", "position_y", "position_z", "momentum_x", "momentum_y", "momentum_z", "qm", "weighting" };
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
        : public amrex::ParIterSoA<RealSoA::nattribs, IntSoA::nattribs>
    {
    public:
        using amrex::ParIterSoA<RealSoA::nattribs, IntSoA::nattribs>::ParIterSoA;

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
        : public amrex::ParConstIterSoA<RealSoA::nattribs, IntSoA::nattribs>
    {
    public:
        using amrex::ParConstIterSoA<RealSoA::nattribs, IntSoA::nattribs>::ParConstIterSoA;

        ParConstIterSoA (ContainerType& pc, int level);

        ParConstIterSoA (ContainerType& pc, int level, amrex::MFItInfo& info);

        ContainerType const * pc () const { return m_pc; }
    };

    /** Beam Particles in ImpactX
     *
     * This class stores particles, distributed over MPI ranks.
     */
    class ImpactXParticleContainer
        : public amrex::ParticleContainerPureSoA<RealSoA::nattribs, IntSoA::nattribs>
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

        /** Check if a container has a ParticleReal component
         *
         * @param name component name to check
         * @return true if found, else false
         */
        bool HasRealComp (std::string const & name);

        /** Check if a container has an Integer component
         *
         * @param name component name to check
         * @return true if found, else false
         */
        bool HasIntComp (std::string const & name);

        /** Get the ParticleReal SoA index of a component
         *
         * This throws a runtime exception if the component does not exist.
         *
         * @param name component name to query index for
         * @return zero-based index
         */
        int GetRealCompIndex (std::string const & name);

        /** Get the Integer SoA index of a component
         *
         * This throws a runtime exception if the component does not exist.
         *
         * @param name component name to query index for
         * @return zero-based index
         */
        int GetIntCompIndex (std::string const & name);

        /** Add a ParticleReal component
         *
         * @param name a unique name of the component, not colliding with RealSoA
         * @param communicate participate in MPI communication when particles move
         */
        void AddRealComp (std::string const & name, bool communicate=true);

        /** Add an Integer component
         *
         * @param name a unique name of the component, not colliding with IntSoA
         * @param communicate participate in MPI communication when particles move
         */
        void AddIntComp (std::string const & name, bool communicate=true);

        /** Add new particles to the container for fixed s.
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
         * @param bchchg total charge within a bunch in C
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
            amrex::ParticleReal bchchg
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

        /** Update reference particle element edge
         *
         */
        void SetRefParticleEdge ();

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
         * @param rho charge grid per level to deposit on
         * @param ref_ratio mesh refinement ratios between levels
         */
        void
        DepositCharge (std::unordered_map<int, amrex::MultiFab> & rho,
                       amrex::Vector<amrex::IntVect> const & ref_ratio);

        /** Get the name of each ParticleReal SoA component */
        std::vector<std::string>
        RealSoA_names () const;

        /** Get the name of each int SoA component */
        std::vector<std::string>
        intSoA_names () const;

        /** Get the current coordinate system of particles in this container */
        CoordSystem
        GetCoordSystem () const;

        /** Set the current coordinate system of particles in this container
         *
         * @param coord_system the new coordinate system the beam particles are in
         */
        void
        SetCoordSystem (CoordSystem coord_system);

      private:

        //! the reference particle for the beam in the particle container
        RefPart m_refpart;

        //! the particle shape
        std::optional<int> m_particle_shape;

        //! a non-owning reference to lost particles, i.e., due to apertures
        ImpactXParticleContainer* m_particles_lost = nullptr;

        //! the current coordinate system of particles in this container
        CoordSystem m_coordsystem = CoordSystem::s;

        //! ParticleReal component names
        std::vector<std::string> m_real_soa_names;

        //! Int component names
        std::vector<std::string> m_int_soa_names;

    }; // ImpactXParticleContainer

} // namespace impactx

#endif // IMPACTX_PARTICLE_CONTAINER_H
