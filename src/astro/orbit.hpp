#ifndef ASTROCPP_ASTRO_ORBIT_HPP
#define ASTROCPP_ASTRO_ORBIT_HPP

#include <cstdio>
#include <cmath>
#include "math/const.hpp"
#include "astro/state_vectors.hpp"

namespace Astro
{
    class EllipticOrbit
    {
        public:
            EllipticOrbit(const double& std_grav_param,
                          const double& semi_major_axis,
                          const double& eccentricity,
                          const double& inclination,         //in radians
                          const double& longitude_of_AN,     //in radians
                          const double& arg_of_periapsis,    //in radians
                          const double& mean_anomaly_at_t0,  //in radians
                          const double& t0);

            void print() const;
            StateVectors state_vectors_from_time(const double& t) const;


        private:
            double mu, a, e, i, W, w, M0, t0,

                   //some pre-computed consts
                   speed_root,      //sqrt(mu / a)
                   mean_motion,     //sqrt(mu / a^3)
                   e_root,          //sqrt(1 - e^2)
                   a_e_root,        //a * e_root
                   orb_to_ecl[6];   //transformation matrix
                                    //from orbital basis to
                                    //ecliptic one
    };

    //EllipticOrbit class methods' implementation
    EllipticOrbit::EllipticOrbit(const double& std_grav_param,
                                 const double& semi_major_axis,
                                 const double& eccentricity,
                                 const double& inclination,         //in radians
                                 const double& longitude_of_AN,     //in radians
                                 const double& arg_of_periapsis,    //in radians
                                 const double& mean_anomaly_at_t0,  //in radians
                                 const double& t0):
        mu(std_grav_param),
        a(semi_major_axis),
        e(eccentricity),
        i(inclination),
        W(longitude_of_AN),
        w(arg_of_periapsis),
        M0(mean_anomaly_at_t0),
        t0(t0)
    {
        const double inv_a = 1 / a;
        speed_root = std::sqrt(mu * inv_a);
        mean_motion = speed_root * inv_a;
        e_root = sqrt(1 - e*e);
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
        orb_to_ecl[1] =  cos_W_sin_w + sin_W_cos_w * cos_i;
        orb_to_ecl[2] = -sin_W_cos_w - cos_W_sin_w * cos_i;
        orb_to_ecl[3] = -sin_W_sin_w + cos_W_cos_w * cos_i;
        orb_to_ecl[4] = -sin_w * sin_i;
        orb_to_ecl[5] =  cos_w * sin_i;
    }

    inline void EllipticOrbit::print() const
    {
        std::printf("\nstd. grav. param. = %lf\n"
                    "semi-major axis = %lf\n"
                    "eccentricity = %lf\n"
                    "inclination = %lf deg\n"
                    "longitude of the AN = %lf deg\n"
                    "arg. of periapsis = %lf deg\n"
                    "mean anomaly at t0 = %lf rad\n"
                    "t0 = %lf\n",
                    mu, a, e, i*Math::RAD_TO_DEG,
                    W*Math::RAD_TO_DEG, w*Math::RAD_TO_DEG, M0, t0);
    }

    StateVectors EllipticOrbit::state_vectors_from_time(const double& t) const
    {
        const double M = mean_motion*(t - t0) + M0; //Mean anomaly

        //Solving Kepler's equation for eccentric anomaly using Newton's method
        const double EPSILON = 1e-14; //precision
        double E0,
               E = M + 0.5*e; //initial guess
        do
        {
            E0 = E;
            E = E0 - (E0 - e*std::sin(E0) - M) / (1 - e*std::cos(E0));
        } while (std::abs(E-E0) > EPSILON);

        //For the sake of performance
        const double sin_E = sin(E),
                     cos_E = cos(E);

        //Radius vector in orbital basis
        const double r_orb_x = a * (cos_E - e),
                     r_orb_y = a_e_root * sin_E;

        //Velocity in orbital basis
        const double v_c = speed_root / (1 - e*cos_E),
                     v_orb_x = -v_c * sin_E,
                     v_orb_y =  v_c * e_root * cos_E;

        //Transform to ecliptic basis
        StateVectors state = {{orb_to_ecl[0]*r_orb_x + orb_to_ecl[1]*r_orb_y,
                               orb_to_ecl[2]*r_orb_x + orb_to_ecl[3]*r_orb_y,
                               orb_to_ecl[4]*r_orb_x + orb_to_ecl[5]*r_orb_y},

                              {orb_to_ecl[0]*v_orb_x + orb_to_ecl[1]*v_orb_y,
                               orb_to_ecl[2]*v_orb_x + orb_to_ecl[3]*v_orb_y,
                               orb_to_ecl[4]*v_orb_x + orb_to_ecl[5]*v_orb_y}};
        return state;
    }
} //namespace Astro

#endif // ASTROCPP_ASTRO_ORBIT_HPP
