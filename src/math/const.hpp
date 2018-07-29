#ifndef ASTROCPP_MATH_CONST_HPP
#define ASTROCPP_MATH_CONST_HPP

namespace astro_cpp::math
{
    constexpr const double PI = 3.14159265358979,
                           INV_PI = 0.3183098861837907,
                           TAU = 6.28318530717959,
                           PI_OVER_2 = 1.570796326794897,
                           PI_SQR = 9.869604401089359,
                           LN_2 = 0.6931471805599453,
                           RAD_TO_DEG = 57.2957795130823,
                           DEG_TO_RAD = 0.0174532925199433;

    // Stolen from StackOverflow
    template<typename T> int sgn(T val)
    {
        return (T(0) < val) - (val < T(0));
    }
} // namespace astro_cpp::math

#endif // ASTROCPP_MATH_CONST_HPP
