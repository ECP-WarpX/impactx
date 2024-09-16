/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "CovarianceMatrixMath.H"

#include <stdexcept>
#include <AMReX_Extension.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuComplex.H>
#include <cmath>
#include <vector>
#include <tuple>


namespace impactx::diagnostics
{

    /** This function returns the three independent kinetic invariants
     *  denoted I2, I4, and I6 as constructed from the 6x6
     *  beam covariance matrix.  These three quantities are invariant
     *  under any linear symplectic transport map, and are used in the
     *  calculation of the three eigenemittances.
     *
     * input - Sigma symmetric 6x6 covariance matrix
     * returns - tuple containing invarants I2, I4, and I6
     */
    std::tuple<
            amrex::ParticleReal,
            amrex::ParticleReal,
            amrex::ParticleReal>
    KineticInvariants (
        amrex::Array2D<amrex::ParticleReal, 1, 6, 1, 6> Sigma
    )
    {
        using namespace amrex::literals;

        std::tuple <amrex::ParticleReal,amrex::ParticleReal,amrex::ParticleReal> invariants;
        amrex::ParticleReal I2 = 0.0_prt;
        amrex::ParticleReal I4 = 0.0_prt;
        amrex::ParticleReal I6 = 0.0_prt;

        // Intermediate matrices used for storage.
        amrex::Array2D<amrex::ParticleReal, 1, 6, 1, 6> S1;
        amrex::Array2D<amrex::ParticleReal, 1, 6, 1, 6> S2;
        amrex::Array2D<amrex::ParticleReal, 1, 6, 1, 6> S4;
        amrex::Array2D<amrex::ParticleReal, 1, 6, 1, 6> S6;

        // Construct the matrix S1 = Sigma*J.  This is a
        // permutation of the columns of Sigma with
        // a change of sign.
        for (int i = 1; i < 7; i++) {
            for (int j = 1; j < 7; j++) {
                if (j % 2) {
                   S1(i,j) = -Sigma(i,j+1); // if j is odd
                }
                else {
                   S1(i,j) = +Sigma(i,j-1); // if j is even
                }
            }
        }

        // Carry out necessary matrix multiplications (3 are needed).
        S2 = impactx::diagnostics::MultiplyMat(S1,S1);
        S4 = impactx::diagnostics::MultiplyMat(S2,S2);
        S6 = impactx::diagnostics::MultiplyMat(S2,S4);

        // Define the three kinematic invariants (should be nonnegative).
        I2 = -impactx::diagnostics::TraceMat(S2)/2.0_prt;
        I4 = +impactx::diagnostics::TraceMat(S4)/2.0_prt;
        I6 = -impactx::diagnostics::TraceMat(S6)/2.0_prt;


        invariants = std::make_tuple(I2,I4,I6);
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
        amrex::Array2D<amrex::ParticleReal, 1, 6, 1, 6> Sigma
    )
    {
        using namespace amrex::literals;

        std::tuple <amrex::ParticleReal,amrex::ParticleReal,amrex::ParticleReal> invariants;
        std::tuple <amrex::ParticleReal,amrex::ParticleReal,amrex::ParticleReal> roots;
        std::tuple <amrex::ParticleReal,amrex::ParticleReal,amrex::ParticleReal> emittances;

        // Get the invariants I2, I4, and I6 from the covariance matrix.
        invariants = KineticInvariants(Sigma);
        amrex::ParticleReal I2 = std::get<0>(invariants);
        amrex::ParticleReal I4 = std::get<1>(invariants);
        amrex::ParticleReal I6 = std::get<2>(invariants);

        // Construct the coefficients of the cubic polynomial
        amrex::ParticleReal a = 1.0_prt;
        amrex::ParticleReal b = -I2;
        amrex::ParticleReal c = (pow(I2,2)-I4)/2.0_prt;
        amrex::ParticleReal d = -pow(I2,3)/6.0_prt + I2*I4/2.0_prt - I6/3.0_prt;

        // Return the cubic coefficients
        //std::cout << "Return a,b,c,d " << a << " " << b << " " << c << " " << d << "\n";

        // Solve for the roots to obtain the eigenemittances.
        // Caution: The order of e1,e2,e3 should be consistent with the
        // order ex,ey,et in the limit of uncoupled transport.
        // The order below is important and has been checked.
        roots = CubicRootsTrig(a,b,c,d);
        amrex::ParticleReal e1 = sqrt(std::get<1>(roots));
        amrex::ParticleReal e2 = sqrt(std::get<2>(roots));
        amrex::ParticleReal e3 = sqrt(std::get<0>(roots));

        emittances = std::make_tuple(e1,e2,e3);
        return emittances;
    }


} // namespace impactx::diagnostics
