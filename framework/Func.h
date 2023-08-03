#pragma once

#include <cmath>

/**
   Polynomial
   Coefficients should be listed from lowest to highest order
*/
template <typename COEFF> double poly(const double &s, const COEFF &c) {
  return c;
  (void)s; // Supresses compiler warnings
}

template <typename COEFF_FIRST, typename... COEFF_REST>
double poly(const double &s, const COEFF_FIRST &c_first,
            const COEFF_REST &...c_rest) {
  return c_first + s * poly(s, c_rest...);
}

// Limited integral
inline double norm_expcosh(const double &tmin, const double &tmax,
                           const double &tau, const double &DG) {
  const double term1_num =
      exp(-tmax / tau) *
      (2.0 * cosh(tmax * DG / 2.0) + DG * tau * sinh(tmax * DG / 2.0));

  const double term2_num =
      exp(-tmin / tau) *
      (2.0 * cosh(tmin * DG / 2.0) + DG * tau * sinh(tmin * DG / 2.0));

  return (2.0 * tau * term1_num - term2_num) / (-4.0 + std::pow(DG * tau, 2));
}

// Limited integral
inline double norm_expsinh(const double &tmin, const double &tmax,
                           const double &tau, const double &DG) {
  const double term1_num =
      exp(-tmax / tau) *
      (DG * tau * cosh(tmax * DG / 2.0) + 2.0 * sinh(tmax * DG / 2.0));

  const double term2_num =
      exp(-tmin / tau) *
      (DG * tau * cosh(tmin * DG / 2.0) + 2.0 * sinh(tmin * DG / 2.0));

  return (2.0 * tau * term1_num - term2_num) / (-4.0 + std::pow(DG * tau, 2));
}
