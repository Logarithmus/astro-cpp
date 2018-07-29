#include "astro/orbit.hpp"

#include <iostream>
#include <cmath>
#include "../math/const.hpp"
#include "../math/vector3d.hpp"

namespace astro_cpp
{
    void StateVectors::print()
    {
        std::cout << "\nRadius vector:";
        r.print();
        std::cout << "Velocity vector:";
        v.print();
    }

    EllipticOrbit::EllipticOrbit(double std_grav_param, double semi_major_axis, double eccentricity,
                                 double inclination,         //in radians
                                 double longitude_of_AN,     //in radians
                                 double arg_of_periapsis,    //in radians
                                 double mean_anomaly_at_t0,  //in radians
                                 double t0):
        mu(std_grav_param),
        a(semi_major_axis),
        e(eccentricity),
        i(inclination),
        W(longitude_of_AN),
        w(arg_of_periapsis),
        M0(mean_anomaly_at_t0),
        t0(t0)
    {
        const double inv_a = 1.0 / a;
        speed_root = std::sqrt(mu * inv_a);
        mean_motion = speed_root * inv_a;
        T = math::TAU / mean_motion;
        e_root = std::sqrt(1.0 - e*e);
        a_e_root = a * e_root;

        using std::sin;
        using std::cos;
        const double sin_W = sin(W),
                     cos_W = cos(W),
                     sin_w = sin(w),
                     cos_w = cos(w),
                     sin_i = sin(i),
                     cos_i = cos(i),
                     cos_W_cos_w = cos_W * cos_w,
                     sin_W_sin_w = sin_W * sin_w,
                     sin_W_cos_w = sin_W * cos_w,
                     cos_W_sin_w = cos_W * sin_w;

        orb_to_ecl[0] =  cos_W_cos_w - sin_W_sin_w * cos_i;
        orb_to_ecl[1] = -cos_W_sin_w - sin_W_cos_w * cos_i;
        orb_to_ecl[2] =  sin_W_cos_w + cos_W_sin_w * cos_i;
        orb_to_ecl[3] = -sin_W_sin_w + cos_W_cos_w * cos_i;
        orb_to_ecl[4] =  sin_w * sin_i;
        orb_to_ecl[5] =  cos_w * sin_i;
    }

    void EllipticOrbit::print() const
    {
        std::cout << "\nstd. grav. param. = " << mu
                  << "\nsemi-major axis = " << a
                  << "\neccentricity = " << e
                  << "\ninclination = " << i * math::RAD_TO_DEG << " deg"
                  << "\nlongitude of the AN = " << W * math::RAD_TO_DEG << " deg"
                  << "\narg. of periapsis = " << w * math::RAD_TO_DEG << " deg"
                  << "\nmean anomaly at t0 = " << M0 << " rad"
                  << "t0 = " << t0 << '\n';
    }

    // Kepler's equation solver (about 40% faster than naive implementation with Newton's method)
    // Based on the paper by F. Landis Markley "Kepler equation solver" (1995)
    double EllipticOrbit::E_from_M(double M) const
    {
        // Wrap M from [0; TAU] to [-PI; +PI]
        if (M > math::PI)
        {
            M -= math::TAU;
        }
        else if (M < -math::PI)
        {
            M += math::TAU;
        }

        // Solving cubic equation based on Pade approximant for sin(E) (more info in the paper)
        const double alpha = (3.0 * math::PI_SQR + 1.6 * math::PI * (math::PI - std::abs(M)) / (1.0 + e)) * alpha_mult,
                     d = 3.0 * (1.0 - e) + alpha * e,
                     M_sqr = M * M,
                     alpha_d = alpha * d,
                     q = 2.0 * alpha_d * (1.0 - e) - M_sqr,
                     q_sqr = q * q,
                     r = 3.0 * alpha_d * (d - 1.0 + e) * M + M_sqr * M,
                     pre_w = std::abs(r) + std::sqrt(q_sqr*q + r*r),
                     w = std::cbrt(pre_w * pre_w);

        // Finally, get our initial guess
        double E = ((2.0 * r * w) / (w*w + w*q + q_sqr) + M) / d;

        // 3rd order Householder's method
        double sin_E = std::sin(E),
               cos_E = std::cos(E),
               f = E - e * sin_E - M,
               Df = 1.0 - e * cos_E,
               DDf = e * sin_E,
               DDDf = e * cos_E,
               Df_Df = Df * Df,
               f_f = f * f;
        E = E - (6.0 * f * Df_Df - 3.0 * f_f * DDf) /
                (6.0 * Df_Df * Df - 6.0 * f * Df * DDf + f_f * DDDf);

        return E;
    }

    StateVectors EllipticOrbit::state_vectors_from_time(double t) const
    {
        double delta_t = std::fmod(t - t0, period());
        double M = mean_motion * delta_t + M0, // Mean anomaly
               E = E_from_M(M);

        // Recompute sin(E) and cos(E) for final value of E
        const double sin_E = std::sin(E),
                     cos_E = std::cos(E);

        // Radius vector in orbital basis
        const double r_orb_x = a * (cos_E - e),
                     r_orb_y = a_e_root * sin_E;

        // Velocity in orbital basis
        const double v_c = speed_root / (1 - e*cos_E),
                     v_orb_x = -v_c * sin_E,
                     v_orb_y =  v_c * e_root * cos_E;

        // Transform to the ecliptical basis
        return {{orb_to_ecl[0]*r_orb_x + orb_to_ecl[1]*r_orb_y,
                 orb_to_ecl[2]*r_orb_x + orb_to_ecl[3]*r_orb_y,
                 orb_to_ecl[4]*r_orb_x + orb_to_ecl[5]*r_orb_y},

                {orb_to_ecl[0]*v_orb_x + orb_to_ecl[1]*v_orb_y,
                 orb_to_ecl[2]*v_orb_x + orb_to_ecl[3]*v_orb_y,
                 orb_to_ecl[4]*v_orb_x + orb_to_ecl[5]*v_orb_y}};
    }
} // namespace astro_cpp

