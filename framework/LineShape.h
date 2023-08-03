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
