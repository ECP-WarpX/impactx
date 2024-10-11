/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef COVARIANCE_MATRIX_MATH_H
#define COVARIANCE_MATRIX_MATH_H

#include <ablastr/constant.H>
#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX_REAL.H>
#include <AMReX_GpuComplex.H>

#include <cmath>
#include <tuple>


namespace impactx::diagnostics
{

    /** Function to return the roots of a cubic polynomial ax^3 + bx^2 + cx + d.
     *  The trigonometric form of Cardano's formula is used.
     *  This implementation expects three real roots, which is verified
     *  by checking the sign of the discriminant.
     *
     * @param[in] a coefficient of cubic term
     * @param[in] b coefficient of quadratic term
     * @param[in] c coefficient of linear term
     * @param[in] d coefficient of constant term
     * @returns tuple of three real roots
     */
    std::tuple<
        amrex::ParticleReal,
        amrex::ParticleReal,
        amrex::ParticleReal
    >
    CubicRootsTrig (
        amrex::ParticleReal a,
        amrex::ParticleReal b,
        amrex::ParticleReal c,
        amrex::ParticleReal d
    )
    {
        using namespace amrex::literals;
        using ablastr::constant::math::pi;

        std::tuple<amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal> roots;
        amrex::ParticleReal x1 = 0.0_prt;
        amrex::ParticleReal x2 = 0.0_prt;
        amrex::ParticleReal x3 = 0.0_prt;

        amrex::ParticleReal Q = (3.0_prt*a*c - std::pow(b,2))/(9.0_prt * std::pow(a,2));
        amrex::ParticleReal R = (9.0_prt*a*b*c - 27_prt*std::pow(a,2)*d - 2.0_prt*std::pow(b,3))/(54.0_prt*std::pow(a,3));
        amrex::ParticleReal discriminant = std::pow(Q,3) + std::pow(R,2);

        // Discriminant should be < 0.  Otherwise, keep theta at default and throw an error.
        amrex::ParticleReal tol = 1.0e-12;  //allow for roundoff error
        if (discriminant > tol) {

            ablastr::warn_manager::WMRecordWarning(
                    "Impactx::diagnostics::CubicRootsTrig",
                    "Polynomial appearing in CubicRootsTrig has one or more complex "
                    "(non-real) roots.  Only the real part is returned.  This "
                    "suggests a loss of numerical precision in computation of the "
                    "eigenemittances.  Treat eigenemittance values with caution.",
                    ablastr::warn_manager::WarnPriority::medium
            );

           std::cout << "Polynomial in CubicRoots has one or more complex roots." << "\n";

        } else if (Q == 0.0_prt) {  // Special case of a triple root

           x1 = - b/(3.0_prt*a);
           x2 = - b/(3.0_prt*a);
           x3 = - b/(3.0_prt*a);

        } else {

           //Three real roots in trigonometric form.
           amrex::ParticleReal theta = std::acos(R/std::sqrt(-std::pow(Q,3)));
           x1 = 2.0_prt*std::sqrt(-Q)*std::cos(theta/3.0_prt) - b/(3.0_prt*a);
           x2 = 2.0_prt*std::sqrt(-Q)*std::cos(theta/3.0_prt + 2.0_prt*pi/3.0_prt) - b/(3.0_prt*a);
           x3 = 2.0_prt*std::sqrt(-Q)*std::cos(theta/3.0_prt + 4.0_prt*pi/3.0_prt) - b/(3.0_prt*a);

        }

        roots = std::make_tuple(x1,x2,x3);
        return roots;
    }


    /** Function to return the roots of a cubic polynomial ax^3 + bx^2 + cx + d.
     *  The algebraic form of Cardano's formula is used.
     *  This implementation expects three real roots, which is verified
     *  by checking the sign of the discriminant.
     *
     * @param[in] a coefficient of cubic term
     * @param[in] b coefficient of quadratic term
     * @param[in] c coefficient of linear term
     * @param[in] d coefficient of constant term
     * @returns tuple of three real roots
     */
    std::tuple<
        amrex::ParticleReal,
        amrex::ParticleReal,
        amrex::ParticleReal
    >
    CubicRootsAlg (
        amrex::ParticleReal a,
        amrex::ParticleReal b,
        amrex::ParticleReal c,
        amrex::ParticleReal d
    )
    {
        using namespace amrex::literals;
        using Complex = amrex::GpuComplex<amrex::ParticleReal>;

        std::tuple<amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal> roots;
        amrex::ParticleReal x1 = 0.0_prt;
        amrex::ParticleReal x2 = 0.0_prt;
        amrex::ParticleReal x3 = 0.0_prt;

        amrex::ParticleReal Q = (3.0_prt*a*c - std::pow(b,2))/(9.0_prt * std::pow(a,2));
        amrex::ParticleReal R = (9.0_prt*a*b*c - 27_prt*std::pow(a,2)*d - 2.0_prt*std::pow(b,3))/(54.0_prt*std::pow(a,3));
        amrex::ParticleReal discriminant = std::pow(Q,3) + std::pow(R,2);

        // Define complex variable C
        Complex Qc(Q,0.0_prt);
        Complex Rc(R,0.0_prt);
        Complex Dc(discriminant,0.0_prt);
        Complex C = amrex::pow(-Rc + amrex::sqrt(Dc),1.0/3.0);

        // Define a cubic root of unity xi
        amrex::ParticleReal xire = -1.0/2.0;
        amrex::ParticleReal xiim = std::sqrt(3.0)/2.0;
        Complex xi(xire,xiim);

        //Three roots in algebraic form.

        if (C.m_real == 0.0_prt && C.m_imag == 0.0_prt) {  // Special case of a triple root

           x1 = - b/(3.0_prt*a);
           x2 = - b/(3.0_prt*a);
           x3 = - b/(3.0_prt*a);

        } else {

           Complex z1 = Qc/C - C;
           Complex z2 = Qc/(xi*C) - xi*C;
           Complex z3 = Qc/(amrex::pow(xi,2)*C) - amrex::pow(xi,2)*C;
           x1 = z2.m_real - b/(3.0*a);
           x2 = z1.m_real - b/(3.0*a);
           x3 = z3.m_real - b/(3.0*a);
        }

        roots = std::make_tuple(x1,x2,x3);
        return roots;
    }

} // namespace impactx::diagnostics

#endif // COVARIANCE_MATRIX_MATH_H
