/*****************************************************************************
 *   Copyright (C) 2004-2018 The pykep development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://gitter.im/esa/pykep                                             *
 *   https://github.com/esa/pykep                                            *
 *                                                                           *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.               *
 *****************************************************************************/

#include "lambert.hpp"

#include <cmath>
#include <iostream>
#include "../math/const.hpp"

namespace astro_cpp
{
    // Computes the derivatives of T (non-dimensional time of flight) with respect to x
    void LambertProblem::dT_over_dx(double x, double T,
                                    double& DT, double& DDT, double& DDDT)
    {
        // Semi-major axis normalized with respect to one corresponding to the min energy ellipse
        double a_m_over_a = 1.0 - x * x,
               a_over_a_m = 1.0 / a_m_over_a,
               inv_y = 1.0 / std::sqrt(1.0 - lambda2 * a_m_over_a),
               inv_y2 = inv_y * inv_y,
               inv_y3 = inv_y2 * inv_y;
        DT = a_over_a_m * (3.0 * T * x - 2.0 + 2.0 * lambda3 * x * inv_y);
        DDT = a_over_a_m * (3.0 * T + 5.0 * x * DT + 2.0 * (1.0 - lambda2) * lambda3 * inv_y3);
        DDDT = a_over_a_m * (7.0 * x * DDT + 8.0 * DT - 6.0 * (1.0 - lambda2) * lambda5 * x * inv_y3 * inv_y2);
    }

    // When computing T(x) in the single revolution case, a loss of precision is encountered when x -> 1.
    // So in this case we use Battin's expression which involves ordinary hypergeometric function
    double LambertProblem::hypergeometric_function(double z, double eps)
    {
        double Sj = 1.0, Cj = 1.0,
               Cj1 = 0.0, Sj1 = 0.0,
               err = 1.0;
        int j = 0;
        while (err > eps)
        {
            Cj1 = Cj * (3.0 + j) * (1.0 + j) / (2.5 + j) * z / (j + 1);
            Sj1 = Sj + Cj1;
            err = std::abs(Cj1);
            Sj = Sj1;
            Cj = Cj1;
            ++j;
        }
        return Sj;
    }

    // Computes T from x and number of revolutions using Lagrange's expression
    // (0.01 < |x - 1| < 0.2)
    double LambertProblem::T_from_x_Lagrange(double x, int revs)
    {
        // a_over_a_m - semi-major axis normalized with respect to one corresponding to the min energy ellipse
        double a_m_over_a = 1.0 - x * x,
               a_over_a_m = 1.0 / a_m_over_a;
        if (a_over_a_m > 0)
        {
            // ellipse
            double alpha = 2.0 * std::acos(x),
                   beta = 2.0 * std::asin(lambda * std::sqrt(a_m_over_a));
            if (lambda < 0.0)
            {
                beta = -beta;
            }
            return 0.5 * (a_over_a_m * std::sqrt(a_over_a_m) * ((alpha - std::sin(alpha)) - (beta - std::sin(beta)) + math::TAU * revs));
        }
        else
        {
            // hyperbola
            double alpha = 2.0 * std::acosh(x),
                   beta = 2.0 * std::asinh(lambda * std::sqrt(-a_m_over_a));
            if (lambda < 0.0)
            {
                beta = -beta;
            }
            return -0.5 * a_over_a_m * std::sqrt(-a_over_a_m) * ((beta - std::sinh(beta)) - (alpha - std::sinh(alpha)));
        }
    }

    double LambertProblem::T_from_x(double x, int revs)
    {
        double battin = 0.01, lagrange = 0.2,
               dist = std::abs(x - 1.0);


        // |x - 1| > 0.2
        if (dist > lagrange)
        {
            // We use Lancaster's expression for T
            double a_m_over_a = 1.0 - x * x,
                   abs_a_m_over_a = std::abs(a_m_over_a),
                   z = std::sqrt(1 - lambda2 * a_m_over_a);
            double y = std::sqrt(abs_a_m_over_a),
                   g = x * z + lambda * a_m_over_a,
                   d;
            if (a_m_over_a >= 0)
            {
                d = std::acos(g);
            }
            else
            {
                double f = y * (z - lambda * x);
                d = std::log(f + g);
            }
            return (x - lambda * z - d / y) / (-a_m_over_a);
        }

        // 0.01 < |x - 1| < 0.2
        if (dist > battin)
        {
            // We use Lagrange's expression for T
            return T_from_x_Lagrange(x, revs);
        }

        // |x - 1| < 0.01
        // We use Battin's series expression for T
        double a_m_over_a = 1.0 - x * x,
               abs_a_m_over_a = std::abs(a_m_over_a),
               z = std::sqrt(1 - lambda2 * a_m_over_a);
        double eta = z - lambda * x,
               S1 = (1.0 - lambda - x * eta) * 0.5,
               Q = hypergeometric_function(S1, 1e-11);
        Q *= (4.0 / 3.0);
        return (eta * eta * eta * Q + 4.0 * lambda * eta) * 0.5 + revs * math::PI / (abs_a_m_over_a * std::sqrt(abs_a_m_over_a));
    }

    double LambertProblem::x_from_T_Householder(double T, double x0, int revs, double eps, int iter_max)
    {
        int iter_count = 0;
        double err = 1.0,
               x = 0.0,
               T_of_x0 = 0.0, delta_T = 0.0,
               DT = 0.0, DDT = 0.0, DDDT = 0.0;
        while ((err > eps) && (iter_count < iter_max))
        {
            T_of_x0 = T_from_x(x0, revs);
            dT_over_dx(x0, T_of_x0, DT, DDT, DDDT);
            delta_T = T_of_x0 - T;
            double DT2 = DT * DT;
            x = x0 - delta_T * (DT2 - delta_T * DDT * 0.5) / (DT * (DT2 - delta_T * DDT) + DDDT * delta_T * delta_T / 6.0);
            err = std::abs(x - x0);
            x0 = x;
            ++iter_count;
        }
        return x;
    }

    LambertProblem::LambertProblem(const Vector3D& r1, const Vector3D& r2,
                                   double tof, double mu, bool is_retrograde, int max_revs)
    {
        // 1 - Getting lambda and T
        const double c = (r2 - r1).mag(), // chord
                     r1_mag = r1.mag(),
                     r2_mag = r2.mag(),
                     inv_r1_mag = 1.0 / r1_mag,
                     inv_r2_mag = 1.0 / r2_mag,
                     s = 0.5 * (c + r1_mag + r2_mag); // semi-perimeter
        const Vector3D r1_hat = inv_r1_mag * r1,
                       r2_hat = inv_r2_mag * r2,
                       h_hat = Vector3D::cross_product(r1, r2).unit_vector();

        const double inv_s = 1.0 / s;
        lambda2 = 1.0 - c * inv_s;
        lambda = std::sqrt(lambda2);

        // tangential unit vectors
        Vector3D t1_hat = Vector3D::cross_product(h_hat, r1_hat);
        Vector3D t2_hat = Vector3D::cross_product(h_hat, r2_hat);
        // if transfer angle is (180; 360) degrees
        // xor orbit is retrograde
        if ((h_hat.z < 0.0) ^ is_retrograde)
        {
            lambda = -lambda;
            t1_hat.invert();
            t2_hat.invert();
        }
        lambda3 = lambda2 * lambda;
        lambda5 = lambda3 * lambda2;
        const double T = tof * std::sqrt(2.0 * mu * inv_s) * inv_s;

        // 2 - Find and adjust max number of revolutions to compute
                     // maximum possible number of revolutions with this T
        int max_possible_revs = static_cast<int>(T * math::INV_PI);
                     // T for x = 0 (minimum energy ellipse) and zero revolutions
        const double T00 = std::acos(lambda) + lambda * std::sqrt(1.0 - lambda2),
                     // T for x = 0 and maximum revolutions
                     T0 = T00 + max_possible_revs * math::PI,
                     // T for x = 1 (parabolic orbit)
                     T1 = (2.0 / 3.0) * (1.0 - lambda3);

        // 2.1 - adjust max_possible_revs using T_min(x_min) found by Halley method
        if ((max_possible_revs <= max_revs) && (max_possible_revs > 0) && (T < T0))
        {
            int it = 0;
            double err = 1.0,
                   T_min = T0,
                   x_old = 0.0, x_new = 0.0,
                   DT, DDT, DDDT;
            while (1)
            {
                dT_over_dx(x_old, T_min, DT, DDT, DDDT);
                if (DT > 1e-14)
                {
                    x_new = x_old - DT * DDT / (DDT * DDT - 0.5 * DT * DDDT);
                }
                err = std::abs(x_old - x_new);
                if ((err < 1e-13) || (it > 12))
                {
                    break;
                }
                T_min = T_from_x(x_new, max_possible_revs);
                x_old = x_new;
                ++it;
            }
            if (T_min > T)
            {
                --max_possible_revs;
            }
        }
        max_revs = std::min(max_possible_revs, max_revs);

        // 2.1 We now reserve the memory for the output variables
        const int vector_size = 2*max_revs + 1;
        std::vector<double> x_vector;
        x_vector.reserve(vector_size);
        velocity.reserve(vector_size);

        // 3 - We may now find all solutions in x,y
        // 3.1 - 0 rev solution
        // 3.1.1 initial guess
        double x0;
        if (T > T00)
        {
            x0 = (T00 - T) / (T - T00 + 4.0); // paper uses another guess
        }
        else if (T < T1)
        {
            x0 = 5.0 * T1 * (T1 - T) / (2.0 * T * (1.0 - lambda5)) + 1.0;
        }
        else
        {
            x0 = std::pow((T / T00), math::LN_2 / std::log(T1 / T00)) - 1.0;
        }

        // 3.1.2 Householder iterations for 0 revolutions
        x_vector.emplace_back(x_from_T_Householder(T, x0, 0, 1e-5, 15));

        // 3.2 multi rev solutions
        for (int i = 1; i <= max_revs; ++i)
        {
            // 3.2.1 left Householder iterations
            x0 = (i * math::PI + math::PI) / (8.0 * T);
            x0 = std::cbrt(x0 * x0);
            x0 = (x0 - 1.0) / (x0 + 1.0);
            x_vector.emplace_back(x_from_T_Householder(T, x0, i, 1e-8, 15));
            // 3.2.2 right Householder iterations
            x0 = (8.0 * T) / (i * math::PI);
            x0 = std::cbrt(x0 * x0);
            x0 = (x0 - 1.0) / (x0 + 1.0);
            x_vector.emplace_back(x_from_T_Householder(T, x0, i, 1e-8, 15));
        }

        // 4 - For each found x value we reconstruct the terminal velocities
        double gamma = std::sqrt(0.5 * mu * s);
        double rho = (r1_mag - r2_mag) / c;
        double sigma = std::sqrt(1 - rho * rho);
        double v1_r, v1_t, v2_r, v2_t, y;
        for (auto x: x_vector)
        {
            y = std::sqrt(1.0 - lambda2 + lambda2 * x * x);
            v1_r =  gamma * ((lambda * y - x) - rho * (lambda * y + x)) * inv_r1_mag;
            v2_r = -gamma * ((lambda * y - x) + rho * (lambda * y + x)) * inv_r2_mag;
            const double v_t = gamma * sigma * (y + lambda * x);
            v1_t = v_t * inv_r1_mag;
            v2_t = v_t * inv_r2_mag;
            velocity.emplace_back(v1_r*r1_hat + v1_t*t1_hat, v2_r*r2_hat + v2_t*t2_hat);
        }
    }
} // namespace astro_cpp
