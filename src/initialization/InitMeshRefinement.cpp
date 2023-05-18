/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell, Remi Lehe
 * License: BSD-3-Clause-LBNL
 */
#include "ImpactX.H"
#include "initialization/InitAmrCore.H"
#include "particles/ImpactXParticleContainer.H"
#include "particles/distribution/Waterbag.H"

#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>

#include <stdexcept>
#include <string>
#include <vector>


namespace impactx
{
    /** Tag cells for refinement.  TagBoxArray tags is built on level lev grids.
     *
     * @TODO handle ngrow argument?
     */
    void ImpactX::ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real /* time */, int /* ngrow */)
    {
        using namespace amrex;

        const auto problo = Geom(lev).ProbLoArray();
        const auto dx = Geom(lev).CellSizeArray();

        amrex::ParmParse pp_geometry("geometry");

        // TODO: this is hard-coded to one level right now
        Vector<Real> lo, hi;
        pp_geometry.getarr("fine_tag_lo", lo);
        pp_geometry.getarr("fine_tag_hi", hi);
        amrex::RealVect fine_tag_lo = RealVect{lo};
        amrex::RealVect fine_tag_hi = RealVect{hi};

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(tags); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.fabbox();
            const auto& fab = tags.array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                RealVect pos {AMREX_D_DECL((i+0.5_rt)*dx[0]+problo[0],
                                           (j+0.5_rt)*dx[1]+problo[1],
                                           (k+0.5_rt)*dx[2]+problo[2])};
                if (pos > fine_tag_lo && pos < fine_tag_hi) {
                    fab(i,j,k) = TagBox::SET;
                }
            });
        }
    }

    /** Make a new level from scratch using provided BoxArray and DistributionMapping.
     *
     * Only used during initialization.
     */
    void ImpactX::MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
                                          const amrex::DistributionMapping& dm)
    {
        amrex::ignore_unused(time);

        // set human-readable tag for each MultiFab
        auto const tag = [lev]( std::string tagname ) {
            tagname.append("[l=").append(std::to_string(lev)).append("]");
            return amrex::MFInfo().SetTag(std::move(tagname));
        };

        // charge (rho) mesh
        amrex::BoxArray const& cba = ba;
        // for MR levels (TODO):
        //cba.coarsen(refRatio(lev - 1));

        // staggering and number of charge components in the field
        auto const rho_nodal_flag = amrex::IntVect::TheNodeVector();
        int const num_components_rho = 1;

        // guard cells for charge deposition
        int const particle_shape = m_particle_container->GetParticleShape();
        int num_guards_rho = 0;
        if (particle_shape % 2 == 0)  // even shape orders
            num_guards_rho = particle_shape / 2 + 1;
        else  // odd shape orders
            num_guards_rho = (particle_shape + 1) / 2;

        m_rho.emplace(
            lev,
            amrex::MultiFab{amrex::convert(cba, rho_nodal_flag), dm, num_components_rho, num_guards_rho, tag("rho")});

        // scalar potential
        auto const phi_nodal_flag = rho_nodal_flag;
        int const num_components_phi = 1;
        int const num_guards_phi = num_guards_rho + 1; // todo: I think this just depends on max(MLMG, force calc)
        m_phi.emplace(
            lev,
            amrex::MultiFab{amrex::convert(cba, phi_nodal_flag), dm, num_components_phi, num_guards_phi, tag("phi")});

        // space charge force
        std::unordered_map<std::string, amrex::MultiFab> f_comp;
        for (std::string const comp : {"x", "y", "z"})
        {
            std::string const str_tag = "space_charge_field_" + comp;
            f_comp.emplace(
                comp,
                amrex::MultiFab{
                    amrex::convert(cba, rho_nodal_flag),
                    dm,
                    num_components_rho,
                    num_guards_rho,
                    tag(str_tag)
                }
            );
        }
        m_space_charge_field.emplace(lev, std::move(f_comp));
    }

    /** Make a new level using provided BoxArray and DistributionMapping and fill
     *  with interpolated coarse level data.
     *
     * @todo this function is not (yet) implemented.
     */
    void ImpactX::MakeNewLevelFromCoarse (int lev, amrex::Real time, const amrex::BoxArray& ba,
                                         const amrex::DistributionMapping& dm)
    {
        // todo
        amrex::ignore_unused(lev, time, ba, dm);
    }

    /** Remake an existing level using provided BoxArray and DistributionMapping
     *  and fill with existing fine and coarse data.
     *
     * @todo this function is not (yet) implemented.
     */
    void ImpactX::RemakeLevel (int lev, amrex::Real time, const amrex::BoxArray& ba,
                              const amrex::DistributionMapping& dm)
    {
        // todo
        amrex::ignore_unused(lev, time, ba, dm);
    }

    /** Delete level data
     */
    void ImpactX::ClearLevel (int lev)
    {
        m_rho.erase(lev);
        m_phi.erase(lev);
        m_space_charge_field.erase(lev);
    }

    void ImpactX::ResizeMesh ()
    {
        BL_PROFILE("ImpactX::ResizeMesh");

        // Extract the min and max of the particle positions
        auto const [x_min, y_min, z_min, x_max, y_max, z_max] = m_particle_container->MinAndMaxPositions();

        // guard for flat beams:
        //   https://github.com/ECP-WarpX/impactx/issues/44
        if (x_min == x_max || y_min == y_max || z_min == z_max)
            throw std::runtime_error("Flat beam detected. This is not yet supported: https://github.com/ECP-WarpX/impactx/issues/44");

        amrex::ParmParse pp_geometry("geometry");
        bool dynamic_size = true;
        pp_geometry.query("dynamic_size", dynamic_size);

        amrex::RealBox rb;
        if (dynamic_size)
        {
            // The box is expanded beyond the min and max of particles.
            // This controlled by the variable `frac` below.
            amrex::Real frac = 3.0;
            pp_geometry.query("prob_relative", frac);

            if (frac < 3.0)
                ablastr::warn_manager::WMRecordWarning(
                    "ImpactX::ResizeMesh",
                    "Dynamic resizing of the mesh uses a geometry.prob_relative "
                    "with less than 3x the beam size. This might result in boundary "
                    "artifacts for space charge calculation. "
                    "There is no minimum good value for this parameter, consider "
                    "doing a convergence test.",
                    ablastr::warn_manager::WarnPriority::high
                );

            if (frac < 1.0)
                throw std::runtime_error("geometry.prob_relative must be >= 1.0 (the beam size) on the coarsest level");

            amrex::RealVect const beam_min(x_min, y_min, z_min);
            amrex::RealVect const beam_max(x_max, y_max, z_max);
            amrex::RealVect const beam_width(beam_max - beam_min);

            amrex::RealVect const beam_padding = beam_width * (frac - 1.0) / 2.0;
            //                           added to the beam extent --^         ^-- box half above/below the beam
            rb.setLo(beam_min - beam_padding);
            rb.setHi(beam_max + beam_padding);
        }
        else
        {
            // note: we read and set the size again because an interactive /
            //       Python user might have changed it between steps
            amrex::Vector<amrex::Real> prob_lo;
            amrex::Vector<amrex::Real> prob_hi;
            pp_geometry.getarr("prob_lo", prob_lo);
            pp_geometry.getarr("prob_hi", prob_hi);

            rb = {prob_lo.data(), prob_hi.data()};
        }

        // updating geometry.prob_lo/hi for consistency
        amrex::Vector<amrex::Real> const prob_lo = {rb.lo()[0], rb.lo()[1], rb.lo()[2]};
        amrex::Vector<amrex::Real> const prob_hi = {rb.hi()[0], rb.hi()[1], rb.hi()[2]};
        pp_geometry.addarr("prob_lo", prob_lo);
        pp_geometry.addarr("prob_hi", prob_hi);

        /* lev==0 */
        {
            // Resize the domain size
            amrex::Geometry::ResetDefaultProbDomain(rb);
            int const lev = 0;

            amrex::Geometry g = Geom(lev);
            g.ProbDomain(rb);
            amrex::AmrMesh::SetGeometry(lev, g);
        }

        for (int lev = 1; lev <= this->max_level; ++lev) {
            // TODO: this is hard-coded to one level right now
            amrex::Vector<amrex::Real> lo, hi;
            pp_geometry.getarr("fine_tag_lo", lo);
            pp_geometry.getarr("fine_tag_hi", hi);
            amrex::RealVect fine_tag_lo = amrex::RealVect{lo};
            amrex::RealVect fine_tag_hi = amrex::RealVect{hi};
            amrex::RealBox rb_lvl1(fine_tag_lo.dataPtr(), fine_tag_hi.dataPtr());

            amrex::Geometry g = Geom(lev);
            g.ProbDomain(rb_lvl1);
            amrex::AmrMesh::SetGeometry(lev, g);
        }

        if (this->max_level > 1)
        {
            amrex::Abort("Did not implement resize for >1 level yet");
        }
    }
} // namespace impactx
