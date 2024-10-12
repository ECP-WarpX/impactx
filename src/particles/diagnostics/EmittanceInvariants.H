/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Marco Garten, Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_EMITTANCE_INVARIANTS
#define IMPACTX_EMITTANCE_INVARIANTS

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_REAL.H>
#include <AMReX_SmallMatrix.H>

#include <string>
#include <tuple>


namespace impactx::diagnostics
{

    /** Returns the three independent kinematic moment invariants
     *  denoted I2, I4, and I6 as constructed from the 6x6
     *  beam covariance matrix.  These three quantities are invariant
     *  under any linear symplectic transport map, and are used in the
     *  calculation of the three eigenemittances, as described in:
     *
     *  G. Rangarajan, F. Neri, and A. Dragt, "Generalized Emittance
     *  Invariants," in Proc. 1989 IEEE Part. Accel. Conf., Chicago, IL,
     *  1989, doi:10.1109/PAC.1989.73422.
     *  A. Dragt, F. Neri, and G. Rangarajan, "General Moment Invariants
     *  for Linear Hamiltonian Systems," Phys. Rev. A 45, 2572-2585 (1992),
     *  doi:10.1103/PhysRevA.45.2572.
     *
     * @param[in] Sigma symmetric 6x6 covariance matrix
     * @returns tuple containing invariants I2, I4, and I6
     */
    std::tuple<
        amrex::ParticleReal,
        amrex::ParticleReal,
        amrex::ParticleReal
    >
    KineticInvariants (
        amrex::SmallMatrix<amrex::ParticleReal, 6, 6, amrex::Order::F, 1> const & Sigma
    );

    /** Returns the three eigenemittances
     *  denoted e1, e2, and e3 as constructed from the 6x6
     *  beam covariance matrix.  These three quantities are invariant
     *  under any linear symplectic transport map, and reduce to
     *  the projected normalized rms emittances in the limit of
     *  uncoupled transport.  These quantities are described in:
     *
     *  G. Rangarajan, F. Neri, and A. Dragt, "Generalized Emittance
     *  Invariants," in Proc. 1989 IEEE Part. Accel. Conf., Chicago, IL,
     *  1989, doi:10.1109/PAC.1989.73422.
     *  A. Dragt, F. Neri, and G. Rangarajan, "General Moment Invariants
     *  for Linear Hamiltonian Systems," Phys. Rev. A 45, 2572-2585 (1992),
     *  doi:10.1103/PhysRevA.45.2572.
     *
     * @param[in] Sigma symmetric 6x6 covariance matrix
     * @returns tuple containing eigenemittances e1, e2, and e3
     */
    std::tuple<
        amrex::ParticleReal,
        amrex::ParticleReal,
        amrex::ParticleReal
    >
    Eigenemittances (
        amrex::SmallMatrix<amrex::ParticleReal, 6, 6, amrex::Order::F, 1> const & Sigma
    );

} // namespace impactx::diagnostics

#endif // IMPACTX_EMITTANCE_INVARIANTS
