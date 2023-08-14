#include <complex>
#include <iostream>
#include <memory>

#include "AmpBase.h"
#include "Event.h"
#include "Minimal4Vector.h"
#include "PropUtil.h"
#include "Propagator.h"

#include "FitParameters.h"

class AmpVV : public AmpBase {
private:
  const unsigned int L_conf;
  std::shared_ptr<Propagator> propagator = nullptr;

public:
  AmpVV(const std::string &name, const unsigned int &L_conf, const Par &idx_abs,
        const Par &idx_arg, const Par &idx_absLambda, const Par &idx_argLambda)
      : AmpBase(name, idx_abs, idx_arg, idx_absLambda, idx_argLambda),
        L_conf(L_conf) {}

  inline void SetPropagator(const PropConf &conf) {
    propagator = NewPropagator(conf);

    propagator->add_parameters(this->par_idx);

    // Used to find out if the parameters have changed
    par =
        std::vector<double>(par_idx.size(), std::numeric_limits<double>::max());
  }

  using AmpBase::GetAmp;
  virtual inline std::complex<double> GetAmp(const Event &event,
                                             const std::vector<double> &par,
                                             const CPConf &CP_conf) const {

    const std::array<Minimal4Vector, nBody> &p4s = event.Getp4(CP_conf);
    const double mass1 = (p4s[0] + p4s[1]).M();
    const double mass2 = (p4s[2] + p4s[3]).M();

    const std::complex<double> propV1{propagator->evaluate(mass1, par)};
    const std::complex<double> propV2{propagator->evaluate(mass2, par)};

    return event.GetVV().GetAmpBarrier(L_conf, CP_conf) *
           event.GetVV().GetAmpSpin(L_conf, CP_conf) * propV1 * propV2;
  }
};
