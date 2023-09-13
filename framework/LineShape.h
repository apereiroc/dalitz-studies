#pragma once

#include "BarrierFactors.h"
#include <complex>

/**
   Breit-Wigner
 */
inline std::complex<double> bw(const double &m, const double &m0,
                               const double &Gamma0, const double &q,
                               const double &q0, const double &Rrad,
                               const unsigned int &L) {
  const double Gamma = Gamma0 * m0 / m * std::pow(q / q0, (2 * L) + 1) *
                       std::pow(bwbf(q, Rrad, L), 2.0) /
                       std::pow(bwbf(q0, Rrad, L), 2.0);

  return m0 * Gamma0 / (std::complex<double>(m0 * m0 - m * m, -m0 * Gamma));
}

/**
  * LASS
*/
inline std::complex<double>
lass(const double &m, const double &m0, const double &Gamma0, const double &q,
     const double q0, const double &a, const double &r) {
  const double Gamma = Gamma0 * m0 / m * q / q0;

  const double deltaR = std::atan(m0 * Gamma / (m0 * m0 - m * m));

  const double inv_arg_atan = 1.0 / a / q + 0.5 * r * q;

  const double deltaF = std::atan(1.0 / inv_arg_atan);

  return std::polar(sin(deltaF), deltaF) +
         std::polar(sin(deltaR), deltaR + 2.0 * deltaF);
}