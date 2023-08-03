#pragma once

#include <cmath>

/**
   Calculate daughter momentum in the 12 reference frame
 */
inline double get_q(const double &m12, const double &m1, const double &m2) {
  const double m12_2 = m12 * m12;
  const double m1_p2_m2 = (m1 + m2) * (m1 + m2);
  const double m1_m2_m2 = (m1 - m2) * (m1 - m2);

  return std::sqrt((m12_2 - m1_p2_m2) * (m12_2 - m1_m2_m2)) / (2.0 * m12);
}
