#ifndef ASTROCPP_ASTRO_ORBIT_HPP
#define ASTROCPP_ASTRO_ORBIT_HPP

#include <cstdio>

namespace Astro
{
    class Orbit
    {
        public:
            Orbit(const double& mu, const double& a, const double& e,
                  const double& i, const double& W, const double& w,
                  const double& M0, const double& t0):
                mu(mu), a(a), e(e), i(i),
                long_of_an(W), arg_of_pe(w),
                M0(M0), t0(t0)
            {

            }

            void print()
            {

            }

        private:
            double mu, a, e, i,
                   long_of_an, arg_of_pe,
                   M0, t0;

    }
} //namespace Astro

#endif // ASTROCPP_ASTRO_ORBIT_HPP
