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
#include <AMReX_Array.H>
#include <AMReX_REAL.H>
#include <tuple>

namespace impactx::diagnostics
{

    /** Function to return the roots of a cubic polynomial ax^3 + bx^2 + cx + d.
     *  The trigonometric form of Cardano's formula is used.
     *  This implementation expects three real roots, which is verified
     *  by checking the sign of the discriminant.
     *
     * @param[in] coefficient a
     * @param[in] coefficient b
     * @param[in] coefficient c
     * @param[in] coefficient d
     * @returns tuple of three real roots
     */
    std::tuple<
            amrex::ParticleReal,
            amrex::ParticleReal,
            amrex::ParticleReal>
    CubicRoots (
        amrex::ParticleReal a,
        amrex::ParticleReal b,
        amrex::ParticleReal c,
        amrex::ParticleReal d
    )
    {
        using namespace amrex::literals;
        using ablastr::constant::math::pi;

        std::tuple <amrex::ParticleReal,amrex::ParticleReal,amrex::ParticleReal> roots;

        amrex::ParticleReal Q = (3.0_prt*a*c - pow(b,2))/(9.0_prt * pow(a,2));
        amrex::ParticleReal R = (9.0_prt*a*b*c - 27_prt*pow(a,2)*d - 2.0_prt*pow(b,3))/(54.0_prt*pow(a,3));
        amrex::ParticleReal discriminant = pow(Q,3) + pow(R,2);

        amrex::ParticleReal theta = 0.0_prt;

        // Discriminant should be < 0.  Otherwise, keep theta at default and throw an error.
        if(discriminant > 0.0_prt){ 
           std::cout << "Polynomial in CubicRoots has one or more complex roots." << "\n";
        } else {
           theta = acos(R/sqrt(-pow(Q,3)));
        }

        //Three real roots in trigonometric form.
        amrex::ParticleReal x1 = 2.0_prt*sqrt(-Q)*cos(theta/3.0_prt) - b/(3.0_prt*a);
        amrex::ParticleReal x2 = 2.0_prt*sqrt(-Q)*cos(theta/3.0_prt + 2.0_prt*pi/3.0_prt) - b/(3.0_prt*a);
        amrex::ParticleReal x3 = 2.0_prt*sqrt(-Q)*cos(theta/3.0_prt + 4.0_prt*pi/3.0_prt) - b/(3.0_prt*a);

        std::cout << "Discriminant, Q, R " << discriminant << " " << Q << " " << R << "\n";
        std::cout << "Return x1, x2, x3 " << x1 << " " << x2 << " " << x3 << "\n";

        roots = std::make_tuple(x1,x2,x3);        
        return roots;
    }


    /** Function to take the trace of a square 6x6 matrix.
     *  
     * @param[in] matrix A
     * @returns the trace of A
     */
    amrex::ParticleReal
    TraceMat (
        amrex::Array2D<amrex::ParticleReal, 1, 6, 1, 6> A
    )
    {
        int const dim = 6;              
        amrex::ParticleReal trA = 0.0;

        for (int i = 1; i < dim+1; i++) {
            trA += A(i,i);         
        }
        return trA;
    }


    /** Function to multiply two square matrices of dimension 6.
     *
     * @param[in] matrix A
     * @param[in] matrix B
     * @returns matrix C = AB
     */
    amrex::Array2D<amrex::ParticleReal, 1, 6, 1, 6>
    MultiplyMat (
        amrex::Array2D<amrex::ParticleReal, 1, 6, 1, 6> A,
        amrex::Array2D<amrex::ParticleReal, 1, 6, 1, 6> B
    )
    {
        amrex::Array2D<amrex::ParticleReal, 1, 6, 1, 6> C;
        int const dim = 6;

        for (int i = 1; i < dim+1; i++) {
            for (int j = 1; j < dim+1; j++) {
                C(i,j) = 0;

                for (int k = 1; k < dim+1; k++) {
                    C(i,j) += A(i,k) * B(k,j);
                }

            }
             
        }
        return C;   
    }


    void TestMult ()
    {
        amrex::Array2D<amrex::ParticleReal, 1, 6, 1, 6> A;
        amrex::Array2D<amrex::ParticleReal, 1, 6, 1, 6> B;
        amrex::Array2D<amrex::ParticleReal, 1, 6, 1, 6> C;

        for (int i = 1; i < 7; i++) {
            for (int j = 1; j < 7; j++) {
                A(i,j) = 0;
                B(i,j) = 0;
                C(i,j) = 0;
            }
        }
        for (int i = 1; i < 7; i++) {
            A(i,i) = 1.0;
            B(i,i) = 1.0;
        }
        A(2,3) = -0.3;
        B(4,5) = 0.5;
        C = MultiplyMat(A,B);
        std::cout << "Matrix multiplication gives " << "\n";
        for (int i = 1; i < 7; i++) {
            for (int j = 1; j < 7; j++) {
                std::cout << i << " " << j << " " << C(i,j) << "\n";
            }
        }


    }



} // namespace impactx::diagnostics

#endif // COVARIANCE_MATRIX_MATH_H
