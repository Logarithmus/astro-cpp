#ifndef ASTROCPP_STATE_VECTORS_HPP
#define ASTROCPP_STATE_VECTORS_HPP

#include <iostream>
#include "math/vector3d.hpp"


namespace astro_cpp
{
    class StateVectors
    {
        public:
            Vector3D r, v;
            void print() const;
            StateVectors propagated(double std_grav_param, double elapsed_time) const;
        private:
            double delta_E_from_delta_M(double delta_M, double e_cos_E1, double e_sin_E1) const;
    };
}

#endif // ASTROCPP_STATE_VECTORS_HPP
