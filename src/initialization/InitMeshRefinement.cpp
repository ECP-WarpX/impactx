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

#include <limits>
#include <stdexcept>
#include <string>
#include <vector>


namespace impactx
{
namespace detail
{
    amrex::Vector<amrex::Real>
    read_mr_prob_relative ()
    {
        amrex::ParmParse pp_amr("amr");
        amrex::ParmParse pp_geometry("geometry");

        int max_level = 0;
        pp_amr.query("max_level", max_level);

        // The box is expanded beyond the min and max of the particle beam.
        amrex::Vector<amrex::Real> prob_relative(max_level + 1, 1.0);
        prob_relative[0] = 3.0;  // top/bottom pad the beam on the lowest level by default by its width
        pp_geometry.queryarr("prob_relative", prob_relative);

        if (prob_relative[0] < 3.0)
            ablastr::warn_manager::WMRecordWarning(
                    "ImpactX::read_mr_prob_relative",
                    "Dynamic resizing of the mesh uses a geometry.prob_relative "
                    "padding of less than 3 for level 0. This might result in boundary "
                    "artifacts for space charge calculation. "
                    "There is no minimum good value for this parameter, consider "
                    "doing a convergence test.",
                    ablastr::warn_manager::WarnPriority::high
            );

        // check that prob_relative[0] > prob_relative[1] > prob_relative[2] ...
        amrex::Real last_lev_rel = std::numeric_limits<amrex::Real>::max();
        for (int lev = 0; lev <= max_level; ++lev) {
            amrex::Real const prob_relative_lvl = prob_relative[lev];
            if (prob_relative_lvl <= 0.0)
                throw std::runtime_error("geometry.prob_relative must be strictly positive for all levels");
            if (prob_relative_lvl > last_lev_rel)
                throw std::runtime_error("geometry.prob_relative must be descending over refinement levels");

            last_lev_rel = prob_relative_lvl;
        }

        return prob_relative;
    }
}
    /** Tag cells on level for refinement.
     *
     * @param lev the current level to refine
     * @param tags the TagBoxArray tags are built on level lev grids
     * @param time current simulation time, unused
     * @param ngrow unused legacy argument that is always zero
     */
    void ImpactX::ErrorEst (
        int lev,
        amrex::TagBoxArray& tags,
        [[maybe_unused]] amrex::Real time,
        [[maybe_unused]] int ngrow
    )
    {
        // level zero is of size in meters (per dimension):
        //   rb_0 = beam_width * prob_relative[lvl=0]
        // level zero is of size in cells:
        //   nx_0 = amr.n_cell
        // level one is of size in meters:
        //   rb_1 = beam_width * prob_relative[lvl=1]
        // level one is of size in cells:
        //   nx_1 = amr.n_cell * rb_1                 / rb_0                 * amr.ref_ratio[lvl=0]
        //        = amr.n_cell * prob_relative[lvl=1] / prob_relative[lvl=0] * amr.ref_ratio[lvl=0]
        // level one lowest/highest cell index to refine:
        //   ref_1_lo = nx_1 / nx_0 / 2
        //   ref_1_hi = ref_1_lo + nx_1
        // ...
        // level n is of size in cells:
        //   nx_n = amr.n_cell   * rb_n                         / rb_0                         * amr.ref_ratio[lvl=n-1] * amr.ref_ratio[lvl=n-2] * ... * amr.ref_ratio[lvl=0]
        //        = amr.n_cell   * prob_relative[lvl=n] / prob_relative[lvl=0] * amr.ref_ratio[lvl=n-1] * amr.ref_ratio[lvl=n-2] * ... * amr.ref_ratio[lvl=0]
        //        = n_cell_{n-1} * prob_relative[lvl=n] / prob_relative[lvl=0] * amr.ref_ratio[lvl=n-1]
        // level n+1 is of size in cells:
        //   nx_nup = n_cell_n * prob_relative[lvl=n+1] / prob_relative[lvl=0] * amr.ref_ratio[lvl=n]
        // level n and n+1 have the following difference in cells, if coarsened to level n
        //   nx_n_diff = (nx_n - nx_nup / amr.ref_ratio[lvl=n])
        // level n lowest/highest cell index to refine:
        //   ref_n_lo = small_end_n + nx_n_diff / 2
        //   ref_n_hi = big_end_n   - nx_n_diff / 2

        const auto dom = Geom(lev).Domain();  // index space of the current level
        auto const r_cell_n = amrex::RealVect(dom.size());  // on the current level
        amrex::RealVect const r_ref_ratio = ref_ratio[lev];

        amrex::Vector<amrex::Real> const prob_relative = detail::read_mr_prob_relative();  // relative padding around beam width

        amrex::RealVect const n_cell_nup = prob_relative[lev+1] / prob_relative[0]
                                           * r_ref_ratio * amrex::RealVect(r_cell_n);
        amrex::RealVect const n_cell_diff = (r_cell_n - n_cell_nup / r_ref_ratio);

        amrex::RealVect const r_fine_tag_lo = amrex::RealVect(dom.smallEnd()) + n_cell_diff / 2.0;
        amrex::RealVect const r_fine_tag_hi = amrex::RealVect(dom.bigEnd())   - n_cell_diff / 2.0;

        amrex::IntVect const fine_tag_lo = {
            AMREX_D_DECL((int)std::ceil(r_fine_tag_lo[0]), (int)std::ceil(r_fine_tag_lo[1]), (int)std::ceil(r_fine_tag_lo[2]))
        };
        amrex::IntVect const fine_tag_hi = {
            AMREX_D_DECL((int)std::ceil(r_fine_tag_hi[0]), (int)std::ceil(r_fine_tag_hi[1]), (int)std::ceil(r_fine_tag_hi[2]))
        };

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(tags); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.fabbox();
            const auto& fab = tags.array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                amrex::IntVect const idx {AMREX_D_DECL(i, j, k)};
                if (idx >= fine_tag_lo && idx <= fine_tag_hi) {
                    fab(i,j,k) = amrex::TagBox::SET;
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

        amrex::Vector<amrex::RealBox> rb(this->finestLevel() + 1);  // extent per level
        if (dynamic_size)
        {
            // The box is expanded (or reduced) relative the min and max of particles.
            auto const prob_relative = detail::read_mr_prob_relative();

            if (prob_relative[0] < 3.0)
                ablastr::warn_manager::WMRecordWarning(
                    "ImpactX::ResizeMesh",
                    "Dynamic resizing of the mesh uses a geometry.prob_relative "
                    "with less than 3x the beam size for level=0. This might result "
                    "in boundary artifacts for space charge calculation. "
                    "There is no minimum good value for this parameter, consider "
                    "doing a convergence test.",
                    ablastr::warn_manager::WarnPriority::high
                );

            if (prob_relative[0] < 1.0)
                throw std::runtime_error("geometry.prob_relative must be >= 1.0 (the beam size) on the coarsest level");

            for (int lev = 0; lev <= this->finestLevel(); ++lev)
            {
                amrex::Real const frac = prob_relative[lev];
                amrex::RealVect const beam_min(x_min, y_min, z_min);
                amrex::RealVect const beam_max(x_max, y_max, z_max);
                amrex::RealVect const beam_width(beam_max - beam_min);

                amrex::RealVect const beam_padding = beam_width * (frac - 1.0) / 2.0;
                //                           added to the beam extent --^         ^-- box half above/below the beam
                rb[lev].setLo(beam_min - beam_padding);
                rb[lev].setHi(beam_max + beam_padding);
            }
        }
        else
        {
            // note: we read and set the size again because an interactive /
            //       Python user might have changed it between steps
            amrex::Vector<amrex::Real> prob_lo;
            amrex::Vector<amrex::Real> prob_hi;
            pp_geometry.getarr("prob_lo", prob_lo);
            pp_geometry.getarr("prob_hi", prob_hi);

            rb[0] = {prob_lo.data(), prob_hi.data()};

            if (this->max_level > 1)
                amrex::Abort("Did not implement ResizeMesh for static domains and >1 MR levels.");
        }

        // updating geometry.prob_lo/hi for consistency
        amrex::Vector<amrex::Real> const prob_lo = {rb[0].lo()[0], rb[0].lo()[1], rb[0].lo()[2]};
        amrex::Vector<amrex::Real> const prob_hi = {rb[0].hi()[0], rb[0].hi()[1], rb[0].hi()[2]};
        pp_geometry.addarr("prob_lo", prob_lo);
        pp_geometry.addarr("prob_hi", prob_hi);

        // Resize the domain size
        amrex::Geometry::ResetDefaultProbDomain(rb[0]);

        for (int lev = 0; lev <= this->finestLevel(); ++lev)
        {
            amrex::Geometry g = Geom(lev);
            g.ProbDomain(rb[lev]);
            SetGeometry(lev, g);

            m_particle_container->SetParticleGeometry(lev, g);
        }
    }
} // namespace impactx
