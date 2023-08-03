#pragma once

#include "Func.h"
#include <cmath>
#include <iostream>

/**
   Blatt-Weisskopf barrier factors
   p:  Momentum
   r:  Meson "radius"
   L:  Relative angular momentum between resonance daughters
 */
inline double bwbf(const double &p, const double &r, const unsigned int &L) {
  const double z = p * p * r * r;

  switch (L) {
  case 0:
    return 1.0;
  case 1:
    return std::sqrt(1.0 / poly(z, 1.0, 1.0));
  case 2:
    return std::sqrt(9.0 / poly(z, 9.0, 3.0, 1.0));
  case 3:
    return std::sqrt(225.0 / poly(z, 225.0, 45.0, 6.0, 1.0));
  case 4:
    return std::sqrt(11025.0 / poly(z, 11025.0, 1575.0, 135.0, 10.0, 1.0));
  case 5:
    return std::sqrt(893025.0 /
                     poly(z, 893025.0, 99225.0, 6300.0, 315.0, 15.0, 1.0));
  default:
    std::cout << "ERROR: L>5" << std::endl;
    exit(1);
  }
}
