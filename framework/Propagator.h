#pragma once

#include "Par.h"
#include <complex>
#include <vector>

/*
 * Base class to implement 2-body mass propagator
 */

class Propagator {
private:
public:
  Propagator() = default;

  virtual std::complex<double>
  evaluate(const double &x, const std::vector<double> &par) const = 0;
  virtual void add_parameters(std::vector<Par> &par_idx,
                              std::vector<double> &par) const = 0;
};
