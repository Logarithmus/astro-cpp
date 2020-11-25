#include "acos_test.hpp"

#include <chrono>
#include <cmath>
#include <iostream>
#include <random>

namespace astro_cpp::tests {
void acos_test(const int N) {
  // std::random_device rand_device;
  // std::mt19937 rand_gen(rand_device());
  // std::uniform_real_distribution<double> cosine(0.99999999, 1.0);

  std::cout.precision(15);
  const double cos_value = 0.999999999262492;
  const double std_acos_value = std::acos(cos_value);
  const double my_acos_value = std::sqrt(2.0 * (1.0 - cos_value));
  std::cout << "cos       = " << cos_value << "\nstd::acos = " << std_acos_value
	    << "\nmy::acos  = " << my_acos_value << '\n';
}
} // namespace astro_cpp::tests
