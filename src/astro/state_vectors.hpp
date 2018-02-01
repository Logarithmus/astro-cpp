#ifndef ASTROCPP_ASTRO_STATE_VECTORS_HPP
#define ASTROCPP_ASTRO_STATE_VECTORS_HPP

#include "math/vector3d.hpp"

namespace Astro
{
    struct StateVectors
    {
        Math::Vector3D r, v;
    };
}//namespace Astro

#endif // ASTROCPP_ASTRO_STATE_VECTORS_HPP
