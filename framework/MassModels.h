#pragma once

#include <complex>
#include <memory>
#include <vector>

#include "Constant.h"
#include "DalitzUtil.h"
#include "LineShape.h"
#include "Par.h"
#include "Propagator.h"

class BreitWigner : public Propagator {
public:
  BreitWigner() = default;

  using Propagator::evaluate;

  [[nodiscard]] inline std::complex<double>
  evaluate(const double &x, const std::vector<double> &par) const override {
    const double mass0 = par[Par::mass_kst892];
    const double width0 = par[Par::width_kst892];
    const double q = get_q(x, mass_Kp, mass_pip);
    const double q0 = get_q(mass0, mass_Kp, mass_pip);

    const std::complex<double> amp = bw(x, mass0, width0, q, q0, Rrad, 1);

    return amp;
  }

  inline void add_parameters(std::vector<Par> &par_idx) const override {
    par_idx.push_back(Par::mass_kst892);
    par_idx.push_back(Par::width_kst892);
  }
};

class LASS : public Propagator {
public:
  LASS() = default;

  [[nodiscard]] inline std::complex<double>
  evaluate(const double &x, const std::vector<double> &par) const override {
    const double mass0 = par[Par::mass_kst1430];
    const double width0 = par[Par::width_kst1430];
    const double a = par[Par::a_lass];
    const double r = par[Par::r_lass];

    const double q = get_q(x, mass_Kp, mass_pip);
    const double q0 = get_q(mass0, mass_Kp, mass_pip);

    const std::complex<double> amp =
            lass(x, mass0, width0, q, q0, a, r);

    return amp;
  }

  inline void add_parameters(std::vector<Par> &par_idx) const override {
    par_idx.push_back(Par::mass_kst1430);
    par_idx.push_back(Par::width_kst1430);
    par_idx.push_back(Par::a_lass);
    par_idx.push_back(Par::r_lass);
  }
};
