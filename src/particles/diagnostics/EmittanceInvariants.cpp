/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "EmittanceInvariants.H"
#include "CovarianceMatrixMath.H"

#include <AMReX_BLProfiler.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuComplex.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <stdexcept>
#include <tuple>
#include <vector>


namespace impactx::diagnostics
{

    /** This function returns the three independent kinetic invariants
     *  denoted I2, I4, and I6 as constructed from the 6x6
     *  beam covariance matrix.  These three quantities are invariant
     *  under any linear symplectic transport map, and are used in the
     *  calculation of the three eigenemittances.
     *
     * input - Sigma symmetric 6x6 covariance matrix
     * returns - tuple containing invariants I2, I4, and I6
     */
    std::tuple<
        amrex::ParticleReal,
        amrex::ParticleReal,
        amrex::ParticleReal
    >
    KineticInvariants (
        amrex::SmallMatrix<amrex::ParticleReal, 6, 6, amrex::Order::F, 1> const & Sigma
    )
    {
        using namespace amrex::literals;

        // Intermediate matrices used for storage.
        amrex::SmallMatrix<amrex::ParticleReal, 6, 6, amrex::Order::F, 1> S1{};

        // Construct the matrix S1 = Sigma*J.  This is a
        // permutation of the columns of Sigma with
        // a change of sign.
        for (int i = 1; i < 7; i++) {
            for (int j = 1; j < 7; j++) {
                if (j % 2 != 0) {
                   S1(i,j) = -Sigma(i,j+1); // if j is odd
                }
                else {
                   S1(i,j) = +Sigma(i,j-1); // if j is even
                }
            }
        }

        // Carry out necessary matrix multiplications (3 are needed).
        auto const S2 = S1 * S1;
        auto const S4 = S2 * S2;
        auto const S6 = S2 * S4;

        // Define the three kinematic invariants (should be nonnegative).
        amrex::ParticleReal const I2 = -S2.trace() / 2.0_prt;
        amrex::ParticleReal const I4 = +S4.trace() / 2.0_prt;
        amrex::ParticleReal const I6 = -S6.trace() / 2.0_prt;


        std::tuple<amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal> invariants = std::make_tuple(I2, I4, I6);
        return invariants;
    }


    /** This function returns the three eigenemittances
     *  denoted e1, e2, and e3 as constructed from the 6x6
     *  beam covariance matrix.  These three quantities are invariant
     *  under any linear symplectic transport map, and reduce to
     *  the projected normalized rms emittances in the limit of
     *  uncoupled transport.
     *
     * input - Sigma symmetric 6x6 covariance matrix
     * returns - tuple containing eigenemittances e1, e2, and e3
     */
    std::tuple<
            amrex::ParticleReal,
            amrex::ParticleReal,
            amrex::ParticleReal>
    Eigenemittances (
        amrex::SmallMatrix<amrex::ParticleReal, 6, 6, amrex::Order::F, 1> const & Sigma
    )
    {
        BL_PROFILE("impactx::diagnostics::Eigenemittances");

        using namespace amrex::literals;

        std::tuple<amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal> invariants;
        std::tuple<amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal> roots;
        std::tuple<amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal> emittances;

        // Get the invariants I2, I4, and I6 from the covariance matrix.
        invariants = KineticInvariants(Sigma);
        amrex::ParticleReal I2 = std::get<0>(invariants);
        amrex::ParticleReal I4 = std::get<1>(invariants);
        amrex::ParticleReal I6 = std::get<2>(invariants);

        // Construct the coefficients of the cubic polynomial.
        // This expression for the characteristic polynomial can be found in:
        // V. Balandin, W. Decking, and N. Golubeva, "Relations Between Projected
        // Emittances and Eigenemittances," in IPAC2013, Shanghai, China, 2013,
        // doi:10.48550/arXiv.1305.1532.
        amrex::ParticleReal a = 1.0_prt;
        amrex::ParticleReal b = -I2;
        amrex::ParticleReal c = (std::pow(I2, 2) - I4) / 2.0_prt;
        amrex::ParticleReal d = -std::pow(I2, 3) / 6.0_prt + I2 * I4 / 2.0_prt - I6 / 3.0_prt;

        // Return the cubic coefficients
        //std::cout << "Return a,b,c,d " << a << " " << b << " " << c << " " << d << "\n";

        // Solve for the roots to obtain the eigenemittances.
        // Caution: The order of e1,e2,e3 should be consistent with the
        // order ex,ey,et in the limit of uncoupled transport.
        // The order below is important and has been checked.
        roots = CubicRootsTrig(a, b, c, d);
        amrex::ParticleReal e1 = std::sqrt(std::abs(std::get<1>(roots)));
        amrex::ParticleReal e2 = std::sqrt(std::abs(std::get<2>(roots)));
        amrex::ParticleReal e3 = std::sqrt(std::abs(std::get<0>(roots)));

        emittances = std::make_tuple(e1,e2,e3);
        return emittances;
    }


} // namespace impactx::diagnostics
