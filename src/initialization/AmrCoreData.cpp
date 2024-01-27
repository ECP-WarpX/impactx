/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell, Ji Qiang
 * License: BSD-3-Clause-LBNL
 */
#include "AmrCoreData.H"

#include "initialization/InitMeshRefinement.H"

#include <AMReX.H>


namespace impactx::initialization
{
    AmrCoreData::AmrCoreData (
        amrex::Geometry const& level_0_geom,
        amrex::AmrInfo const& amr_info
    )
        : amrex::AmrCore(level_0_geom, amr_info)
    {
    }

    AmrCoreData::AmrCoreData (
        amrex::RealBox const & rb,
        int max_level_in,
        amrex::Vector<int> const & n_cell_in,
        int coord,
        amrex::Vector<amrex::IntVect> const & ref_ratios,
        amrex::Array<int,AMREX_SPACEDIM> const & is_per
    )
        : amrex::AmrCore(rb, max_level_in, n_cell_in, coord, ref_ratios, is_per)
    {
    }

    void
    AmrCoreData::ErrorEst (
        int lev,
        amrex::TagBoxArray& tags,
        [[maybe_unused]] amrex::Real time,
        [[maybe_unused]] int ngrow)
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
        //   nx_n = amr.n_cell   * rb_n                 / rb_0                         * amr.ref_ratio[lvl=n-1] * amr.ref_ratio[lvl=n-2] * ... * amr.ref_ratio[lvl=0]
        //        = amr.n_cell   * prob_relative[lvl=n] / prob_relative[lvl=0] * amr.ref_ratio[lvl=n-1] * amr.ref_ratio[lvl=n-2] * ... * amr.ref_ratio[lvl=0]
        //        = n_cell_{n-1} * prob_relative[lvl=n] / prob_relative[lvl=0] * amr.ref_ratio[lvl=n-1]
        // level n+1 is of size in cells:
        //   nx_nup = n_cell_n * prob_relative[lvl=n+1] / prob_relative[lvl=0] * amr.ref_ratio[lvl=n]
        // level n and n+1 have the following difference in cells, if coarsened to level n
        //   nx_n_diff = (nx_n - nx_nup / amr.ref_ratio[lvl=n])
        // level n lowest/highest cell index to refine:
        //   ref_n_lo = small_end_n + nx_n_diff / 2
        //   ref_n_hi = big_end_n   - nx_n_diff / 2

        const auto dom = boxArray(lev).minimalBox();  // covered index space of the current level
        auto const r_cell_n = amrex::RealVect(dom.size());  // on the current level
        amrex::RealVect const r_ref_ratio = ref_ratio[lev];

        // level width relative to beam width
        amrex::Vector<amrex::Real> const prob_relative = read_mr_prob_relative();

        amrex::RealVect const n_cell_nup = amrex::RealVect(r_cell_n) * r_ref_ratio
                                           * prob_relative[lev+1] / prob_relative[0];
        amrex::RealVect const n_cell_diff = (r_cell_n - n_cell_nup / r_ref_ratio);

        amrex::RealVect const r_fine_tag_lo = amrex::RealVect(dom.smallEnd()) + n_cell_diff / 2.0;
        amrex::RealVect const r_fine_tag_hi = amrex::RealVect(dom.bigEnd())   - n_cell_diff / 2.0;

        amrex::IntVect const fine_tag_lo = {
                AMREX_D_DECL((int)std::ceil(r_fine_tag_lo[0]), (int)std::ceil(r_fine_tag_lo[1]), (int)std::ceil(r_fine_tag_lo[2]))
        };
        amrex::IntVect const fine_tag_hi = {
                AMREX_D_DECL((int)std::ceil(r_fine_tag_hi[0]), (int)std::ceil(r_fine_tag_hi[1]), (int)std::ceil(r_fine_tag_hi[2]))
        };

        amrex::Box const fine_tag = {fine_tag_lo, fine_tag_hi};

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
                if (fine_tag.contains(idx)) {
                    fab(i,j,k) = amrex::TagBox::SET;
                }
            });
        }
    }

    void
    AmrCoreData::MakeNewLevelFromScratch (
        [[maybe_unused]] int lev,
        [[maybe_unused]] amrex::Real time,
        [[maybe_unused]] const amrex::BoxArray& ba,
        [[maybe_unused]] const amrex::DistributionMapping& dm)
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

    void
    AmrCoreData::MakeNewLevelFromCoarse (
        [[maybe_unused]] int lev,
        [[maybe_unused]] amrex::Real time,
        [[maybe_unused]] const amrex::BoxArray& ba,
        [[maybe_unused]] const amrex::DistributionMapping& dm)
    {
        amrex::Abort("MakeNewLevelFromCoarse: Not Implemented Yet");
    }

    void
    AmrCoreData::RemakeLevel (
        [[maybe_unused]] int lev,
        [[maybe_unused]] amrex::Real time,
        [[maybe_unused]] const amrex::BoxArray& ba,
        [[maybe_unused]] const amrex::DistributionMapping& dm)
    {
        amrex::Abort("RemakeLevel: Not Implemented Yet");
    }

    void
    AmrCoreData::ClearLevel ([[maybe_unused]] int lev)
    {
        m_rho.erase(lev);
        m_phi.erase(lev);
        m_space_charge_field.erase(lev);
    }
} // namespace impactx::initialization
