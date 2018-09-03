#ifndef ASTROCPP_ANOMALIES_HPP
#define ASTROCPP_ANOMALIES_HPP

namespace astro_cpp::anomalies
{
    // Eccentric from mean
    double E_from_M(double M, double e);
    double H_from_M(double M, double e);
    double delta_E_from_delta_M(double delta_E, double e);
    double delta_H_from_delta_M(double delta_E, double e);
}

#endif // ASTROCPP_ANOMALIES_HPP
