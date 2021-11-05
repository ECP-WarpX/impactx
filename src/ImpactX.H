/* Copyright 2021 Axel Huebl, Chad Mitchell, Ji Qiang
 *
 * This file is part of ImpactX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACT_X_H
#define IMPACT_X_H

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_AmrCore.H>

#include <memory>

namespace impactx
{
    /** Overwrite defaults in AMReX Inputs
     *
     * This overwrites defaults in amrex::ParamParse for inputs.
     */
    void
    overwrite_amrex_parser_defaults();

    /** An ImpactX simulation
     *
     * This is the central ImpactX simulation class
     */
    class ImpactX final
        : public amrex::AmrCore
    {
      public:
        ImpactX () = delete;
        ImpactX (amrex::Geometry const& geom, amrex::AmrInfo const& amr_info);
        ~ImpactX () = default;
        ImpactX (ImpactX const&) = delete;
        ImpactX (ImpactX &&) = delete;
        void operator= (ImpactX const&) = delete;
        void operator= (ImpactX &&) = delete;

        /** Initialize particle containers and MultiFabs (fields) */
        void initData ();

        /** Run the main simulation loop */
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

        /** these are the physical/beam particles of the simulation */
        std::unique_ptr<ImpactXParticleContainer> mypc;
    };

} // namespace impactx

#endif // IMPACT_X_H
