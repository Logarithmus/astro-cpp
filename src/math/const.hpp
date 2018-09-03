#ifndef ASTROCPP_MATH_CONST_HPP
#define ASTROCPP_MATH_CONST_HPP

#include <cmath>

namespace astro_cpp::math
{
    constexpr double PI = 3.14159265358979,
                     INV_PI = 1 / PI,
                     TAU = 2 * PI,
                     PI_OVER_2 = PI / 2,
                     PI_SQR = PI * PI,
                     LN_2 = std::log(2),
                     RAD_TO_DEG = 180 / PI,
                     DEG_TO_RAD = PI / 180;

    // Stolen from StackOverflow
    template<typename T>
    inline int sgn(T val)
    {
        return (T(0) < val) - (val < T(0));
    }
} // namespace astro_cpp::math

#endif // ASTROCPP_MATH_CONST_HPP
