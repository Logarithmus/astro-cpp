#ifndef ASTROCPP_PROPAGATE_STATE_VECTORS_HPP
#define ASTROCPP_PROPAGATE_STATE_VECTORS_HPP

#include "orbit/state_vectors.hpp"

namespace astro_cpp
{
    StateVectors propagate_state_vectors(const StateVectors& state1, double delta_t, double mu);
}

#endif // ASTROCPP_PROPAGATE_STATE_VECTORS_HPP
