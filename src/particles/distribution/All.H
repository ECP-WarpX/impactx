/* Copyright 2022 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_DISTRIBUTION_ALL_H
#define IMPACTX_DISTRIBUTION_ALL_H

#include "Gaussian.H"
#include "Kurth4D.H"
#include "Kurth6D.H"
#include "KVdist.H"
#include "None.H"
#include "Semigaussian.H"
#include "Waterbag.H"

#include <variant>


namespace impactx
{
namespace distribution
{
    using KnownDistributions = std::variant<
            None, /* must be first, so KnownDistributions creates a default constructor */
            Gaussian, Kurth4D, Kurth6D, KVdist, Semigaussian, Waterbag>;

} // namespace distribution
} // namespace impactx

#endif // IMPACTX_DISTRIBUTION_ALL_H