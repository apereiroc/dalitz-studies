//
// Created by Asier Pereiro Castro on 11/9/23.
//

#pragma once

#include <complex>
#include <iostream>
#include <memory>

#include "AmpBase.h"
#include "Event.h"
#include "Minimal4Vector.h"
#include "PropUtil.h"
#include "Propagator.h"

class AmpSS : public AmpBase {
private:
    std::shared_ptr<Propagator> propagator = nullptr;

public:
    AmpSS(std::string name, const Par &idx_abs, const Par &idx_arg,
          const Par &idx_absLambda, const Par &idx_argLambda)
            : AmpBase(std::move(name), idx_abs, idx_arg,
                      idx_absLambda, idx_argLambda) {}

    inline void SetPropagator(const PropConf &conf) {
        propagator = NewPropagator(conf);

        propagator->add_parameters(this->par_idx);

        // Used to find out if the parameters have changed
        par = std::vector<double>(par_idx.size(),
                                  std::numeric_limits<double>::max());
    }

    inline std::complex<double> GetAmp(const Event &event,
                                       const std::vector<double> &par,
                                       const CPConf &CP_conf) const override {

        const std::array<Minimal4Vector, nBody> &p4s = event.Getp4(CP_conf);
        const double mass1 = (p4s[0] + p4s[1]).M();
        const double mass2 = (p4s[2] + p4s[3]).M();

        const std::complex<double> propS1{propagator->evaluate(mass1, par)};
        const std::complex<double> propS2{propagator->evaluate(mass2, par)};

        return propS1 * propS2;
    }
};
