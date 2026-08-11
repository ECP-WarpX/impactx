/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_TRANSFORM_LATTICE_INSERT_H
#define IMPACTX_TRANSFORM_LATTICE_INSERT_H

#include "elements/All.H"
#include "Lattice.H"

#include <AMReX_REAL.H>

#include <list>


namespace impactx::elements::transformation
{
    /** Insert an element every s into a lattice
     *
     * Splits up every element that is on s = N * ds for N>0.
     *
     * @param[in] lattice the lattice to insert into
     * @param[in] ds spacing in meters along s to add an element
     * @param[in] element the extra element to add every s
     * @return a new lattice with the inserted elements
     */
    Lattice
    insert_element_every_ds (
        Lattice const & lattice,
        amrex::ParticleReal ds,
        elements::KnownElements element
    );

} // namespace impactx::elements::transformation

#endif // IMPACTX_TRANSFORM_LATTICE_INSERT_H
