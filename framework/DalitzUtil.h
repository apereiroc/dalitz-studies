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

/**
 * Calculate jacobian transformation from phase space
 */

inline double get_detJ(const double &m12, const double &m34,
                       const double &m1234 = mass_Bs,
                       const double &m1 = mass_Kp, const double &m2 = mass_pip,
                       const double &m3 = mass_Kp,
                       const double &m4 = mass_pip) {
  return
          get_q(m1234, m12, m34) *
          get_q(m12, m1, m2) *
          get_q(m34, m3, m4);
}

