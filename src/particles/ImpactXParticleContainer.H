/* Copyright 2022 The Regents of the University of California, through Lawrence
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
            x,  ///< position in x [m] (at fixed t OR fixed z)
            y,  ///< position in y [m] (at fixed t OR fixed z)
            z,  ///< position in z [m] (at fixed t) OR time-of-flight ct [m] (at fixed z)
        };
    };

    /** This struct indexes the additional Real attributes
     *  stored in an SoA in ImpactXParticleContainer
     */
    struct RealSoA
    {
        enum
        {
            ux,  ///< momentum in x, normalized to proper velocity
            uy,  ///< momentum in y, normalized to proper velocity
            pt,  ///< kinetic energy deviation, (E0-E)/mc^2
            m_qm, ///< charge to mass ratio, in q_e/m_e (q_e/eV)
            w,   ///< particle weight, unitless
            nattribs ///< the number of attributes above (always last)
        };
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

    /** This struct stores the reference particle attributes
     *  stored in ImpactXParticleContainer
     */
    struct RefPart
    {
        amrex::ParticleReal s = 0.0;  ///< integrated orbit path length, in meters
        amrex::ParticleReal x = 0.0;  ///< horizontal position x, in meters
        amrex::ParticleReal y = 0.0;  ///< vertical position y, in meters
        amrex::ParticleReal z = 0.0;  ///< longitudinal position y, in meters
        amrex::ParticleReal t = 0.0;  ///< clock time * c in meters
        amrex::ParticleReal px = 0.0; ///< momentum in x, normalized to proper velocity
        amrex::ParticleReal py = 0.0; ///< momentum in y, normalized to proper velocity
        amrex::ParticleReal pz = 0.0; ///< momentum in z, normalized to proper velocity
        amrex::ParticleReal pt = 0.0; ///< energy deviation, normalized by rest energy

        amrex::Array2D<amrex::ParticleReal,1,6,1,6> map; ///< linearized map

        /** Get reference particle relativistic gamma
         *
         * @returns relativistic gamma
         */
        amrex::ParticleReal
        gamma () const;

        /** Get reference particle relativistic beta
         *
         * @returns relativistic beta
         */
        amrex::ParticleReal
        beta () const;

        /** Get reference particle beta*gamma
         *
         * @returns relativistic beta*gamma
         */
        amrex::ParticleReal
        beta_gamma () const;

        /** Get reference particle energy
         *
         * @returns kinetic energy in MeV
         */
        amrex::ParticleReal
        energy_MeV (amrex::ParticleReal massE) const;

        /** Set reference particle energy
         *
         * @param energy initial kinetic energy (MeV)
         * @param massE particle rest energy (MeV)
         */
        void
        set_energy_MeV (
            amrex::ParticleReal const energy,
            amrex::ParticleReal const massE);

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
        using iterator = amrex::ParIter<0, 0, RealSoA::nattribs, IntSoA::nattribs>;

        //! amrex constant iterator for particle boxes (read-only)
        using const_iterator = amrex::ParConstIter<0, 0, RealSoA::nattribs, IntSoA::nattribs>;

        //! Construct a new particle container
        ImpactXParticleContainer (amrex::AmrCore* amr_core);

        //! Destruct a particle container
        virtual ~ImpactXParticleContainer() = default;

        /** Add new particles to the container
         *
         * Note: This can only be used *after* the initialization (grids) have
         *       been created, meaning after the call to AmrCore::InitFromScratch
         *       or AmrCore::InitFromCheckpoint has been made in the ImpactX
         *       class.
         *
         * @param lev mesh-refinement level
         * @param x positions in x
         * @param y positions in y
         * @param z positions in z
         * @param px momentum in x
         * @param py momentum in y
         * @param pz momentum in z
         * @param qm charge over mass in 1/eV
         * @param bchchg total charge within a bunch in C
         */
        void
        AddNParticles (int lev,
                       amrex::Vector<amrex::ParticleReal> const & x,
                       amrex::Vector<amrex::ParticleReal> const & y,
                       amrex::Vector<amrex::ParticleReal> const & z,
                       amrex::Vector<amrex::ParticleReal> const & px,
                       amrex::Vector<amrex::ParticleReal> const & py,
                       amrex::Vector<amrex::ParticleReal> const & pz,
                       amrex::ParticleReal const & qm,
                       amrex::ParticleReal const & bchchg);

        /** Set reference particle attributes
         *
         * @param refpart reference particle
         */
        void
        SetRefParticle (RefPart const refpart);

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
        RefPart const
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
        void SetParticleShape (int const order);

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

      private:

        //! the reference particle for the beam in the particle container
        RefPart m_refpart;

        //! the particle shape
        std::optional<int> m_particle_shape;

    }; // ImpactXParticleContainer

} // namespace impactx

#endif // IMPACTX_PARTICLE_CONTAINER_H
