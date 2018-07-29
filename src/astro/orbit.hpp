#ifndef ASTROCPP_ORBIT_HPP
#define ASTROCPP_ORBIT_HPP

#include "math/vector3d.hpp"
#include "math/const.hpp"

namespace astro_cpp
{
    struct StateVectors
    {
        Vector3D r, v;
        void print();
    };

    class EllipticOrbit
    {
        public:
            EllipticOrbit(double std_grav_param, double semi_major_axis, double eccentricity,
                          double inclination,         //in radians
                          double longitude_of_AN,     //in radians
                          double arg_of_periapsis,    //in radians
                          double mean_anomaly_at_t0,  //in radians
                          double t0);
            void print() const;
            double period() const;
            StateVectors state_vectors_from_time(double t) const;

        private:
            static const constexpr double alpha_mult = 1.0 / (math::PI_SQR - 6);
            double mu, a, e, i, W, w, M0, t0,

                   // some pre-computed consts
                   speed_root,      // sqrt(mu / a)
                   mean_motion,     // sqrt(mu / a^3)
                   T,               // orbital period
                   e_root,          // sqrt(1 - e^2)
                   a_e_root,        // a * e_root
                   orb_to_ecl[6];   // transformation matrix
                                    // from the orbital basis to
                                    // the ecliptic one

            double E_from_M(double M) const;
    };

    inline double EllipticOrbit::period() const
    {
        return T;
    }
} // namespace astro_cpp

#endif // ASTROCPP_ORBIT_HPP
