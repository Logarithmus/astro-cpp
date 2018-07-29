#ifndef ASTROCPP_ACOS_TEST_HPP
#define ASTROCPP_ACOS_TEST_HPP

#include <iostream>
#include <cmath>
#include "math/const.hpp"

namespace astro_cpp::tests
{
    template<typename T> int sgn(T val)
    {
        return (T(0) < val) - (val < T(0));
    }

    inline double fast_acos_nvidia(double x)
    {
          double negate = double(x < 0);
          x = std::abs(x);
          double ret = -0.0187293;
          ret *= x;
          ret += 0.0742610;
          ret *= x;
          ret -= 0.2121144;
          ret *= x;
          ret += 1.5707288;
          ret *= std::sqrt(1.0 - x);
          ret = ret - 2 * negate * ret;
          return negate * 3.14159265358979 + ret;
    }

    inline double fast_acos(double x)
    {
        static const int size = 10000;
        static const double step = 1.0 / size;
        static const struct LookupTable
        {
            double value[size + 1];
            LookupTable()
            {
                for (int i = 0; i <= size; ++i)
                {
                    value[i] = std::asin(step * i);
                }
            }
        } lut;

        double abs_x = std::abs(x);
        int i = static_cast<int>(abs_x * size);
        //std::cout << i << std::endl;
        double asin_x = sgn(x) * (lut.value[i + 1] - lut.value[i]) * size * (abs_x - i * step) + lut.value[i];
        return M_PI_2 - asin_x;
    }

    inline double fast_acos_arora_russel(double x)
    {
        double abs_x = std::abs(x);
        double tmp;
        if (abs_x < 0.6)
        {
            tmp = (0.000014773722 + (1.1782782 - 0.52020038 * abs_x) * abs_x) /
                  (1.1793469 + (-0.53277664 - 0.14454764 * abs_x) * abs_x);
        }
        else if (abs_x < 0.97)
        {
            tmp = (0.011101554f + (8.9810074 + (-14.816468 + 5.9249913 * abs_x) * abs_x) * abs_x) /
                  (9.2299851f + (-16.001036 + 6.8381053 * abs_x) * abs_x);
        }
        else if (abs_x < 0.99)
        {
            tmp = (-35.750586f + (107.24325 - 70.780244 * abs_x) * abs_x) /
                  (27.105764 - 26.638535 * abs_x);
        }
        else
        {
            tmp = std::asin(abs_x);
        }

        return math::PI_OVER_2 - math::sgn(x) * tmp;
    }

    void acos_test(const int N);
} // namespace astro_cpp::tests

#endif // ASTROCPP_ACOS_TEST_HPP
