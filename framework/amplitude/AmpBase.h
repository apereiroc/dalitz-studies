#pragma once

#include <complex>
#include <memory>
#include <string>
#include <vector>

#include "CPUtil.h"
#include "Event.h"
#include "Par.h"
#include "Propagator.h"

class AmpBase {

private:
  const std::string name;
  const Par idx_abs, idx_arg, idx_absLambda, idx_argLambda;

  mutable bool par_change = true;
  std::vector<Par> par_idx;
  mutable std::vector<double> par;

public:
  AmpBase(const std::string &name, const Par &idx_abs, const Par &idx_arg,
          const Par &idx_absLambda, const Par &idx_argLambda)
      : name(name), idx_abs(idx_abs), idx_arg(idx_arg),
        idx_absLambda(idx_absLambda), idx_argLambda(idx_argLambda) {
    // Store parameter indices
    par_idx.push_back(idx_abs);
    par_idx.push_back(idx_arg);
    par_idx.push_back(idx_absLambda);
    par_idx.push_back(idx_argLambda);

    // Used to find out if the parameters have changed
    par =
        std::vector<double>(par_idx.size(), std::numeric_limits<double>::max());
  }

  // Pure virtual method to obtain the dynamical component
  inline std::complex<double> GetAmp(const Event &event,
                                     const std::vector<double> &par,
                                     const CPConf &CP_conf) = 0;

  // Method to obtain the strong+weak coupling
  // Common to all derived classes
  inline std::complex<double> GetCoupling(const std::vector<double> &par,
                                          const CPConf &CP_conf) {
    std::complex<double> weak_coupling = {1.0, 0.0};
    if (CP_conf == CPConf::Abar)
      weak_coupling = std::polar(par[idx_absLambda], par[idx_argLambda]);

    return std::polar(par[idx_abs], par[idx_arg]) * weak_coupling;
  }

  inline void UpdatePar(const std::vector<double> &par) const {
    this->par_change = false;
    for (unsigned int i = 0; i < par_idx.size(); i++) {
      if (par[this->par_idx[i]] != this->par[i]) {
        this->par[i] = par[this->par_idx[i]];
        this->par_change = true;
      }
    }
  }

  inline bool ChangePar() const { return this->par_change; }
};
