/* Copyright 2022 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell, Ji Qiang
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACT_X_H
#define IMPACT_X_H

#include "particles/distribution/All.H"
#include "particles/elements/All.H"
#include "particles/ImpactXParticleContainer.H"

#include <AMReX_AmrCore.H>

#include <list>
#include <memory>
#include <unordered_map>


namespace impactx
{
    /** An ImpactX simulation
     *
     * This is the central ImpactX simulation class
     */
    class ImpactX final
        : public amrex::AmrCore
    {
      public:
        /** Construct an ImpactX simulation object
         */
        ImpactX ();

        ~ImpactX () = default;

        // removed constructors/assignments
        ImpactX (ImpactX const&) = delete;
        ImpactX (ImpactX &&) = delete;
        void operator= (ImpactX const&) = delete;
        void operator= (ImpactX &&) = delete;

        /** Initialize AMReX blocks/grids
         *
         * This must come first, before particle beams and lattice elements are
         * initialized.
         */
        void initGrids ();

        /** Initialize the particle beam distribution
         *
         * This clears and initialized the particle beam from input file and
         * command-line options, as parsed by amrex::ParmParse.
         */
        void initBeamDistributionFromInputs ();

        /** Initialize the list of lattice elements
         *
         * This clears and initialized the accelerator lattice elements from
         * input file and command-line options, as parsed by amrex::ParmParse.
         */
        void initLatticeElementsFromInputs ();

        /** Generate and add n particles to the particle container
         *
         * Will also resize the geometry based on the updated particle
         * distribution's extent and then redistribute particles in according
         * AMReX grid boxes.
         *
         * @param qm charge/mass ratio (q_e/eV)
         * @param bunch_charge bunch charge (C)
         * @param distr distribution function to draw from (object)
         * @param npart number of particles to draw
         */
        void
        add_particles (
            amrex::ParticleReal qm,
            amrex::ParticleReal bunch_charge,
            distribution::KnownDistributions distr,
            int npart
        );

        /** Run the main simulation loop for a number of steps
         */
        void evolve ();

      private:
        //! Tag cells for refinement.  TagBoxArray tags is built on level lev grids.
        void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time,
                               int ngrow) override;

        //! Make a new level from scratch using provided BoxArray and DistributionMapping.
        //! Only used during initialization.
        void MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
                                              const amrex::DistributionMapping& dm) override;

        //! Make a new level using provided BoxArray and DistributionMapping and fill
        //  with interpolated coarse level data.
        void MakeNewLevelFromCoarse (int lev, amrex::Real time, const amrex::BoxArray& ba,
                                             const amrex::DistributionMapping& dm) override;

        //! Remake an existing level using provided BoxArray and DistributionMapping
        //  and fill with existing fine and coarse data.
        void RemakeLevel (int lev, amrex::Real time, const amrex::BoxArray& ba,
                                  const amrex::DistributionMapping& dm) override;

        //! Delete level data
        void ClearLevel (int lev) override;

      public:
        /** Resize the mesh, based on the extent of the bunch of particle
         *
         * This only changes the physical extent of the mesh, but not the
         * number of grid cells.
         */
        void ResizeMesh ();

        /** these are the physical/beam particles of the simulation */
        std::unique_ptr<ImpactXParticleContainer> m_particle_container;

        /** charge per level */
        std::unordered_map<int, amrex::MultiFab> m_rho;

        /** these are elements defining the lattice */
        std::list<KnownElements> m_lattice;

        std::list<std::variant<int>> m_lattice_test;
    };

} // namespace impactx

#endif // IMPACT_X_H
