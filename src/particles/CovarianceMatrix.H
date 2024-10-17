/* Copyright 2022-2024 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_DISTRIBUTION_COVARIANCE_MATRIX_H
#define IMPACTX_DISTRIBUTION_COVARIANCE_MATRIX_H

#include <AMReX_REAL.H>
#include <AMReX_SmallMatrix.H>


namespace impactx
{
    /** this is a 6x6 matrix */
    using CovarianceMatrix = amrex::SmallMatrix<amrex::ParticleReal, 6, 6, amrex::Order::F, 1>;

} // namespace impactx::distribution

#endif // IMPACTX_DISTRIBUTION_COVARIANCE_MATRIX_H
