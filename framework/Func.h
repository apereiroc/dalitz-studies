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

// Time integrals
double norm_expcosh(const double &t, const double &tau, const double &DGon2) {
  const double tau2 = std::pow(tau, 2.0);

  return (((((tau2 * DGon2) + tau) * exp(2.0 * DGon2 * t)) - (tau2 * DGon2) +
           tau) *
          exp((-DGon2 * t) - (t / tau))) /
         ((2.0 * tau2 * std::pow(DGon2, 2.0)) - 2.0);
}

double norm_expcosh(const double &t_min, const double &t_max, const double &tau,
                    const double &DeltaGamma) {
  return norm_expcosh(t_max, tau, DeltaGamma / 2.0) -
         norm_expcosh(t_min, tau, DeltaGamma / 2.0);
}

double norm_expsinh(const double &t, const double &tau, const double &DGon2) {
  return (tau *
          ((((tau * DGon2) + 1.0) * exp(2.0 * DGon2 * t)) + (tau * DGon2) -
           1.0) *
          exp((-DGon2 * t) - (t / tau))) /
         (2.0 * (std::pow(tau * DGon2, 2.0) - 1.0));
}

double norm_expsinh(const double &t_min, const double &t_max, const double &tau,
                    const double &DeltaGamma) {
  return norm_expsinh(t_max, tau, DeltaGamma / 2.0) -
         norm_expsinh(t_min, tau, DeltaGamma / 2.0);
}
