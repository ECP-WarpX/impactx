/* Copyright 2021 Axel Huebl
 *
 * This file is part of ImpactX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_PARTICLE_CONTAINER_H
#define IMPACTX_PARTICLE_CONTAINER_H

#include <AMReX_BaseFwd.H>
#include <AMReX_AmrCoreFwd.H>
#include <AMReX_ParIter.H>
#include <AMReX_Particles.H>

#include <AMReX_Vector.H>

#include <tuple>


namespace impactx
{
    /** This struct indexes the additional Real attributes
     *  stored in an SoA in ImpactXParticleContainer
     */
    struct RealSoA
    {
        enum
        {
            ux,  ///< momentum in x, normalized to proper velocity
            uy,  ///< momentum in y, normalized to proper velocity
            t,   ///< phase in arrival time * 2pi * reference frequency, radians
            pt,  ///< kinetic energy deviation, (E0-E)/mc^2
            q_m, ///< charge to mass ratio, in q_e/m_e (?)
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

    /** This struct indexes the reference particle attributes
     *  stored in ImpactXParticleContainer
     */
    struct RefPart
    {
        amrex::ParticleReal x = 0.0;  ///< horizontal position x, in meters
        amrex::ParticleReal y = 0.0;  ///< vertical position y, in meters
        amrex::ParticleReal t = 0.0;  ///< arrival time * c in meters
        amrex::ParticleReal px = 0.0; ///< momentum in x, normalized to proper velocity
        amrex::ParticleReal py = 0.0; ///< momentum in y, normalized to proper velocity
        amrex::ParticleReal pt = 0.0; ///< energy deviation, normalized by rest energy
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
         */
        void
        AddNParticles (int lev,
                       amrex::Vector<amrex::ParticleReal> const & x,
                       amrex::Vector<amrex::ParticleReal> const & y,
                       amrex::Vector<amrex::ParticleReal> const & z,
                       amrex::Vector<amrex::ParticleReal> const & px,
                       amrex::Vector<amrex::ParticleReal> const & py,
                       amrex::Vector<amrex::ParticleReal> const & pz);

        void
        SetRefParticle (RefPart const refpart);


        /** Compute the min and max of the particle position in each dimension
         *
         * @returns x_min, x_max, y_min, y_max, z_min, z_max
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

      private:

        //! the reference particle for the beam in the particle container
        RefPart m_refpart;

    }; // ImpactXParticleContainer

} // namespace impactx

#endif // IMPACTX_PARTICLE_CONTAINER_H
