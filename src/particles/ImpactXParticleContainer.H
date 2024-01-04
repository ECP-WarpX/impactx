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

#include <AMReX_AmrCoreFwd.H>
#include <AMReX_BaseFwd.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParIter.H>
#include <AMReX_Particles.H>

#include <AMReX_IntVect.H>
#include <AMReX_Vector.H>

#include <optional>
#include <tuple>
#include <unordered_map>


namespace impactx
{
    /** AMReX pre-defined Real attributes
     *
     * These are the AMReX pre-defined struct indexes for the Real attributes
     * stored in an AoS in ImpactXParticleContainer. We document this here,
     * because we change the meaning of these "positions" depending on the
     * coordinate system we are currently in.
     */
    struct RealAoS
    {
        enum
        {
            x,  ///< position in x [m] (at fixed s OR fixed t)
            y,  ///< position in y [m] (at fixed s OR fixed t)
            t,  ///< c * time-of-flight [m] (at fixed s)
            nattribs ///< the number of attributes above (always last)
        };

        // at fixed t, the third component represents the position z
        enum {
            z = t  ///< position in z [m] (at fixed t)
        };

        //! named labels for fixed s
        static constexpr auto names_s = { "position_x", "position_y", "position_t" };
        //! named labels for fixed t
        static constexpr auto names_t = { "position_x", "position_y", "position_z" };
        static_assert(names_s.size() == nattribs);
        static_assert(names_t.size() == nattribs);
    };

    /** This struct indexes the additional Real attributes
     *  stored in an SoA in ImpactXParticleContainer
     */
    struct RealSoA
    {
        enum
        {
            px,  ///< momentum in x, scaled by the magnitude of the reference momentum [unitless] (at fixed s or t)
            py,  ///< momentum in y, scaled by the magnitude of the reference momentum [unitless] (at fixed s or t)
            pt,  ///< energy deviation, scaled by speed of light * the magnitude of the reference momentum [unitless] (at fixed s)
            qm,  ///< charge to mass ratio, in q_e/m_e [q_e/eV]
            w,   ///< particle weight, number of real particles represented by this macroparticle [unitless]
            nattribs ///< the number of attributes above (always last)
        };

        // at fixed t, the third component represents the momentum in z
        enum {
            pz = pt  ///< momentum in z, scaled by the magnitude of the reference momentum [unitless] (at fixed t)
        };

        //! named labels for fixed s
        static constexpr auto names_s = { "momentum_x", "momentum_y", "momentum_t", "qm", "weighting" };
        //! named labels for fixed t
        static constexpr auto names_t = { "momentum_x", "momentum_y", "momentum_z", "qm", "weighting" };
        static_assert(names_s.size() == nattribs);
        static_assert(names_t.size() == nattribs);
    };

    /** This struct indexes the additional Integer attributes
     *  stored in an SoA in ImpactXParticleContainer
     */
    struct IntSoA
    {
        enum
        {
            nattribs ///< the number of particles above (always last)
        };
    };

    /** AMReX iterator for particle boxes
     *
     * We subclass here to change the default threading strategy, which is
     * `static` in AMReX, to `dynamic` in ImpactX.
     */
    class ParIter
        : public amrex::ParIter<0, 0, RealSoA::nattribs, IntSoA::nattribs>
    {
    public:
        using amrex::ParIter<0, 0, RealSoA::nattribs, IntSoA::nattribs>::ParIter;

        ParIter (ContainerType& pc, int level);

        ParIter (ContainerType& pc, int level, amrex::MFItInfo& info);
    };

    /** Const AMReX iterator for particle boxes - data is read only.
     *
     * We subclass here to change the default threading strategy, which is
     * `static` in AMReX, to `dynamic` in ImpactX.
     */
    class ParConstIter
        : public amrex::ParConstIter<0, 0, RealSoA::nattribs, IntSoA::nattribs>
    {
    public:
        using amrex::ParConstIter<0, 0, RealSoA::nattribs, IntSoA::nattribs>::ParConstIter;

        ParConstIter (ContainerType& pc, int level);

        ParConstIter (ContainerType& pc, int level, amrex::MFItInfo& info);
    };

    /** Beam Particles in ImpactX
     *
     * This class stores particles, distributed over MPI ranks.
     */
    class ImpactXParticleContainer
        : public amrex::ParticleContainer<0, 0, RealSoA::nattribs, IntSoA::nattribs>
    {
    public:
        //! amrex iterator for particle boxes
        using iterator = impactx::ParIter;

        //! amrex constant iterator for particle boxes (read-only)
        using const_iterator = impactx::ParConstIter;

        //! Construct a new particle container
        ImpactXParticleContainer (amrex::AmrCore* amr_core);

        //! Destruct a particle container
        virtual ~ImpactXParticleContainer() = default;

        /** Add new particles to the container for fixed s.
         *
         * In this version of AddNParticles, when running on GPUs, the input
         * vectors of particle data live on the host, and will be copied to
         * the device as part of the method call.
         *
         * Note: This can only be used *after* the initialization (grids) have
         *       been created, meaning after the call to AmrCore::InitFromScratch
         *       or AmrCore::InitFromCheckpoint has been made in the ImpactX
         *       class.
         *
         * @param lev mesh-refinement level
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
        AddNParticles (int lev,
                       amrex::Vector<amrex::ParticleReal> const & x,
                       amrex::Vector<amrex::ParticleReal> const & y,
                       amrex::Vector<amrex::ParticleReal> const & t,
                       amrex::Vector<amrex::ParticleReal> const & px,
                       amrex::Vector<amrex::ParticleReal> const & py,
                       amrex::Vector<amrex::ParticleReal> const & pt,
                       amrex::ParticleReal qm,
                       amrex::ParticleReal bchchg);

        /** Add new particles to the container for fixed s.
         *
         * In this version of AddNParticles, when running on GPUs, the input
         * vectors of particle data have already been copied to the device.
         *
         * Note: This can only be used *after* the initialization (grids) have
         *       been created, meaning after the call to AmrCore::InitFromScratch
         *       or AmrCore::InitFromCheckpoint has been made in the ImpactX
         *       class.
         *
         * @param lev mesh-refinement level
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
        AddNParticles (int lev,
                       amrex::Gpu::DeviceVector<amrex::ParticleReal> const & x,
                       amrex::Gpu::DeviceVector<amrex::ParticleReal> const & y,
                       amrex::Gpu::DeviceVector<amrex::ParticleReal> const & t,
                       amrex::Gpu::DeviceVector<amrex::ParticleReal> const & px,
                       amrex::Gpu::DeviceVector<amrex::ParticleReal> const & py,
                       amrex::Gpu::DeviceVector<amrex::ParticleReal> const & pt,
                       amrex::ParticleReal qm,
                       amrex::ParticleReal bchchg);

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

        /** Get the name of each Real AoS component */
        std::vector<std::string>
        RealAoS_names () const;

        /** Get the name of each Real SoA component */
        std::vector<std::string>
        RealSoA_names () const;

      private:

        //! the reference particle for the beam in the particle container
        RefPart m_refpart;

        //! the particle shape
        std::optional<int> m_particle_shape;

        //! a non-owning reference to lost particles, i.e., due to apertures
        ImpactXParticleContainer* m_particles_lost = nullptr;

    }; // ImpactXParticleContainer

    /** Get the name of each Real AoS component */
    std::vector<std::string>
    get_RealAoS_names ();

    /** Get the name of each Real SoA component
     *
     * @param num_real_comps number of compile-time + runtime arrays
     * @return names
     */
    std::vector<std::string>
    get_RealSoA_names (int num_real_comps);

} // namespace impactx

#endif // IMPACTX_PARTICLE_CONTAINER_H
