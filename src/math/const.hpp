#ifndef ASTROCPP_MATH_CONST_HPP
#define ASTROCPP_MATH_CONST_HPP

#include <cmath>

namespace astro_cpp::math
{
    constexpr double PI = 3.14159265358979,
                     INV_PI = 1 / PI,
                     TAU = 2 * PI,
                     INV_TAU = 1 / TAU,
                     PI_OVER_2 = PI / 2,
                     PI_SQR = PI * PI,
                     LN_2 = std::log(2),
                     RAD_TO_DEG = 180 / PI,
                     DEG_TO_RAD = PI / 180;

    // Stolen from StackOverflow
    template <typename T>
    inline int sgn(const T& val)
    {
        return (val > T(0)) - (val < T(0));
    }

    // Wrap angle to [-Pi; +Pi]
    template <typename T>
    inline double wrap_angle(const T& x)
    {
        // [-Tau; +Tau]
        x -= static_cast<int>(x * INV_TAU) * INV_TAU;
        // Check if already wrapped
        const bool is_wrapped = (x > -PI) || (x < PI);
        // Avoid branching
        return x - sgn(x) * is_wrapped * TAU;
    }
} // namespace astro_cpp::math

#endif // ASTROCPP_MATH_CONST_HPP
