#include "state_vectors.hpp"

namespace astro_cpp
{
    void StateVectors::print() const
    {
        std::cout << "\nRadius vector:";
        r.print();
        std::cout << "Velocity vector:";
        v.print();
    }

    double StateVectors::delta_E_from_delta_M(double delta_M, double e_cos_E1, double e_sin_E1) const
    {
        double delta_E = delta_M,
               delta_E0;
        auto f = [&e_cos_E1, &e_sin_E1, &delta_M](double delta_E)
        {
            return delta_E - e_cos_E1 * std::sin(delta_E) + e_sin_E1 * (1.0 - std::cos(delta_E)) - delta_M;
        };
        auto df = [&e_cos_E1, &e_sin_E1, &delta_M](double delta_E)
        {
            return 1.0 - e_cos_E1 * std::cos(delta_E) + e_sin_E1 * std::sin(delta_E);
        };

        int iters = 0;
        do
        {
            delta_E0 = delta_E;
            delta_E = delta_E0 - f(delta_E0) / df(delta_E0);
            ++iters;
        } while (std::abs(delta_E0 - delta_E) > 1e-14);
        if (iters > 5)
        {
            std::cout << iters << '\n';
        }

        return delta_E;
    }

    StateVectors StateVectors::propagated(double std_grav_param, double elapsed_time) const
    {
        // Shortcuts for function arguments
        double& mu = std_grav_param;
        double& delta_t = elapsed_time;

        // |r|, |v|^2, semi-major axis
        const double r1_mag = r.mag(),
                     a = (mu * r1_mag) / (2.0 * mu - r1_mag * v.mag_sqr());

        double f, g, f_dot, g_dot; // Lagrange coefficients
        if (a > 0)
        { // Elliptic motion
            // Some intermediate quantities
            const double root = std::sqrt(mu * a),
                         inv_root = 1 / std::sqrt(mu * a),
                         e_cos_E1 = 1 - r1_mag / a,
                         e_sin_E1 = Vector3D::dot_product(r, v) * inv_root;

            // Difference in mean and eccentric anomalies
            const double delta_M = std::sqrt(mu / (a * a * a)) * delta_t,
                         delta_E = delta_E_from_delta_M(delta_M, e_cos_E1, e_sin_E1);

            // Compute sin(delta_E) and cos(delta_E) only one time
            const double sin_delta_E = std::sin(delta_E),
                         cos_delta_E = std::cos(delta_E);

            // |r2|
            const double r2_mag = a * (1.0 - e_cos_E1 * cos_delta_E + e_sin_E1 * sin_delta_E);

            // Finally, calculate Lagrange coefficients
            f = 1 - (a / r1_mag) * (1 - cos_delta_E);
            g = a * inv_root * (r1_mag * sin_delta_E + a * e_sin_E1 * (1 - cos_delta_E));
            f_dot = -root * sin_delta_E / (r1_mag * r2_mag);
            g_dot = 1 - (a / r2_mag) * (1 - cos_delta_E);
        }
        else
        { // Hyperbolic motion
            // W.I.P.
            f = 0;
            g = 0;
            f_dot = 0;
            g_dot = 0;
        }

        // [r2, v2] = linear combination of [r1, v1] with Lagrange coefficients
        return { {f * r + g * v}, {f_dot * r + g_dot * v} };
    }
} // namespace astro_cpp
