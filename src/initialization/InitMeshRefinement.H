/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#pragma once

#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>

#include <limits>
#include <stdexcept>
#include <string>
#include <vector>


namespace impactx::initialization
{
    amrex::Vector<amrex::Real>
    read_mr_prob_relative ()
    {
        amrex::ParmParse pp_algo("algo");
        amrex::ParmParse pp_amr("amr");
        amrex::ParmParse pp_geometry("geometry");

        bool space_charge = false;
        pp_algo.queryAdd("space_charge", space_charge);

        std::string poisson_solver = "multigrid";
        pp_algo.queryAdd("poisson_solver", poisson_solver);

        int max_level = 0;
        pp_amr.query("max_level", max_level);

        if (max_level > 1 && !space_charge)
            throw std::runtime_error(
                "Mesh-refinement (amr.max_level>=0) is only supported with "
                "space charge modeling (algo.space_charge=1).");

        // The box is expanded beyond the min and max of the particle beam.
        amrex::Vector<amrex::Real> prob_relative(max_level + 1, 1.0);
        prob_relative[0] = 3.0;  // top/bottom pad the beam on the lowest level by default by its width
        pp_geometry.queryarr("prob_relative", prob_relative);

        if (prob_relative[0] < 3.0 && space_charge && poisson_solver == "multigrid")
            ablastr::warn_manager::WMRecordWarning(
                    "ImpactX::read_mr_prob_relative",
                    "Dynamic resizing of the mesh uses a geometry.prob_relative "
                    "padding of less than 3 for level 0. This might result in boundary "
                    "artifacts for space charge calculation. "
                    "There is no minimum good value for this parameter, consider "
                    "doing a convergence test.",
                    ablastr::warn_manager::WarnPriority::high
            );

        if (prob_relative[0] < 1.0)
            throw std::runtime_error("geometry.prob_relative must be >= 1.0 (the beam size) on the coarsest level");

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
} // namespace impactx::initialization
