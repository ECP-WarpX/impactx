/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
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

#include "initialization/AmrCoreData.H"

#include <AMReX_REAL.H>

#include <list>
#include <memory>


namespace impactx
{
    /** An ImpactX simulation
     *
     * This is the central ImpactX simulation class
     */
    class ImpactX
    {
      public:
        /** Construct an ImpactX simulation object
         */
        ImpactX ();

        // removed constructors/assignments
        ImpactX (ImpactX const&) = delete;
        ImpactX (ImpactX &&) = delete;
        void operator= (ImpactX const&) = delete;
        void operator= (ImpactX &&) = delete;

        /** Destruct the ImpactX simulation object
         */
        ~ImpactX ();

        /** Deallocate all contexts and data
         */
        void finalize ();

        /** Initialize AMReX blocks/grids for domain decomposition & space charge mesh.
         *
         * This must come first, before particle beams and lattice elements are
         * initialized.
         */
        void init_grids ();

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
         * @param bunch_charge bunch charge (C)
         * @param distr distribution function to draw from (object)
         * @param npart number of particles to draw
         */
        void
        add_particles (
            amrex::ParticleReal bunch_charge,
            distribution::KnownDistributions distr,
            int npart
        );

        /** Validate the simulation is ready to run via @see evolve
         */
        void validate ();

        /** Check inputs for unused parameters (e.g. typos) after step 1
         *
         * @return true if successful
         */
        bool early_param_check ();

        /** Run the main simulation loop for a number of steps
         */
        void evolve ();

        /** Query input for warning logger variables and set up warning logger accordingly
         *
         * Input variables are: ``always_warn_immediately`` and ``abort_on_warning_threshold``.
         */
        void init_warning_logger ();

        /** The AMReX core object that contains AMR hierarchy, \see amrex::AmrCore
         *
         * Inside here, we store our particle containers and space charge fields.
         */
        std::unique_ptr<initialization::AmrCoreData> amr_data;

        /** Resize the mesh, based on the extent of the bunch of particle
         *
         * This only changes the physical extent of the mesh, but not the
         * number of grid cells.
         */
        void ResizeMesh ();

        /** these are elements defining the accelerator lattice */
        std::list<KnownElements> m_lattice;

        /** Was init_grids already called?
         *
         * Some operations, like resizing a simulation in terms of cells and changing blocking
         * factors are not possible after they were initialized in AMReX. This keeps track
         * of this.
         */
        bool initialized ()
        {
            return m_grids_initialized;
        }

      private:
        /** Keeps track if init_grids was called.
         *
         * Some operations, like resizing a simulation in terms of cells and changing blocking
         * factors are not possible after they were initialized in AMReX. This keeps track
         * of this.
         */
        bool m_grids_initialized = false;
    };

} // namespace impactx

#endif // IMPACT_X_H
