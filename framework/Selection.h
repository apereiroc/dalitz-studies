#pragma once

#include "Constant.h"

inline bool Selection_masses(const double &mass1, const double &mass2) {
  return mass1 < kpi_max and mass2 < kpi_max;
}

// vim: fdm=marker ts=2 sw=2 sts=2 sr noet
