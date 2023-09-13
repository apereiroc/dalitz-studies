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

enum VSConf {
    Plus, Minus
};

class AmpVS : public AmpBase {
private:
    const VSConf VS_conf;

    std::shared_ptr<Propagator> propagatorV = nullptr;
    std::shared_ptr<Propagator> propagatorS = nullptr;

public:
    AmpVS(std::string name, const VSConf &VS_conf, const Par &idx_abs,
          const Par &idx_arg, const Par &idx_absLambda,
          const Par &idx_argLambda)
            : AmpBase(std::move(name), idx_abs, idx_arg, idx_absLambda,
                      idx_argLambda),
              VS_conf(VS_conf) {}

    inline void SetPropagatorV(const PropConf &conf) {
        propagatorV = NewPropagator(conf);

        propagatorV->add_parameters(this->par_idx);

        // Used to find out if the parameters have changed
        par = std::vector<double>(par_idx.size(),
                                  std::numeric_limits<double>::max());
    }

    inline void SetPropagatorS(const PropConf &conf) {
        propagatorS = NewPropagator(conf);

        propagatorS->add_parameters(this->par_idx);

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

        const std::complex<double> propV1{propagatorV->evaluate(mass1, par)};
        const std::complex<double> propS1{propagatorS->evaluate(mass1, par)};
        const std::complex<double> propV2{propagatorV->evaluate(mass2, par)};
        const std::complex<double> propS2{propagatorS->evaluate(mass2, par)};

        const std::complex<double> amp_VS =
                event.GetVS().GetAmpBarrier(TopoVSConf::VS, CP_conf) *
                event.GetVS().GetAmpSpin(TopoVSConf::VS, CP_conf) *
                propV1 * propS2;

        const std::complex<double> amp_SV =
                event.GetVS().GetAmpBarrier(TopoVSConf::SV, CP_conf) *
                event.GetVS().GetAmpSpin(TopoVSConf::SV, CP_conf) *
                propS1 * propV2;

        switch (VS_conf) {
            case VSConf::Plus:
                return (amp_VS + amp_SV) / M_SQRT2;
            case VSConf::Minus:
                return (amp_VS - amp_SV) / M_SQRT2;
            default:
                std::cout << "ERROR: Invalid _VS_conf = " << VS_conf << " set"
                          << std::endl;
                std::exit(1);
        }
    }
};
