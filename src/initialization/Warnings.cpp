/* Copyright 2022 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Marco Garten
 * License: BSD-3-Clause-LBNL
 */
#include "ImpactX.H"

#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX_BLProfiler.H>
#include <AMReX_ParmParse.H>

#include <optional>
#include <stdexcept>
#include <string>


namespace impactx
{
void ImpactX::init_warning_logger ()
{
    amrex::ParmParse pp_impactx("impactx");

    // Set the flag to control if ImpactX has to emit a warning message
    // as soon as a warning is recorded
    bool always_warn_immediately = false;
    pp_impactx.query("always_warn_immediately", always_warn_immediately);
    ablastr::warn_manager::GetWMInstance().
        SetAlwaysWarnImmediately(always_warn_immediately);

    // Set the WarnPriority threshold to decide if ImpactX has to abort
    // when a warning is recorded
    if(std::string str_abort_on_warning_threshold = "";
            pp_impactx.query("abort_on_warning_threshold", str_abort_on_warning_threshold)){
        std::optional<ablastr::warn_manager::WarnPriority> abort_on_warning_threshold = std::nullopt;
        if (str_abort_on_warning_threshold == "high")
            abort_on_warning_threshold = ablastr::warn_manager::WarnPriority::high;
        else if (str_abort_on_warning_threshold == "medium" )
            abort_on_warning_threshold = ablastr::warn_manager::WarnPriority::medium;
        else if (str_abort_on_warning_threshold == "low")
            abort_on_warning_threshold = ablastr::warn_manager::WarnPriority::low;
        else {
            throw std::runtime_error(str_abort_on_warning_threshold
                +"is not a valid option for impactx.abort_on_warning_threshold (use: low, medium or high)");
        }
        ablastr::warn_manager::GetWMInstance().
            SetAbortThreshold(abort_on_warning_threshold);
    }
}

bool ImpactX::early_param_check ()
{
    BL_PROFILE("ImpactX::early_param_check");

    amrex::Print() << "\n";
    amrex::ParmParse().QueryUnusedInputs();

    // Print the warning list right after the first step.
    amrex::Print() << ablastr::warn_manager::GetWMInstance()
            .PrintGlobalWarnings("FIRST STEP");

    return true;
}
} // namespace impactx
