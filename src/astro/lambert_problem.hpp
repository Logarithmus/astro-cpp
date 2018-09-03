#ifndef ASTROCPP_LAMBERT_HPP
#define ASTROCPP_LAMBERT_HPP

#include <vector>
#include <utility>
#include "../math/vector3d.hpp"

namespace astro_cpp
{
    class LambertProblem
    {
        public:
            LambertProblem(const Vector3D& r1, const Vector3D& r2,
                           double tof, double mu, bool is_retrograde, int max_revs);
            std::vector<std::pair<Vector3D, Vector3D>> velocity;

        private:
            double lambda, lambda2, lambda3, lambda5;

            void dT_over_dx(double x, double T,
                            double& DT, double& DDT, double& DDDT);
            double hypergeometric_function(double z, double eps);
            double T_from_x_Lagrange(double x, int revs);
            double T_from_x(double x, int revs);
            double x_from_T_Householder(double T, double x0, int revs, double eps, int iter_max);
    };
} // namespace astro_cpp

#endif // ASTROCPP_LAMBERT_HPP
