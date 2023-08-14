#pragma once

#include "FitParameters.h"
#include "Par.h"
#include "TRandom3.h"

inline void RandomiseCovariantCoeffs(FitParameters &mn_param,
                                     const unsigned int &seed) {
  const auto rand = std::make_unique<TRandom3>(seed);

  for (unsigned int idx = Par::abs_VV_S; idx <= Par::arg_VV_D; idx++) {
    if (mn_param.Parameter(idx).IsFixed() == false) {
      const double ll = mn_param.Parameter(idx).LowerLimit();
      const double ul = mn_param.Parameter(idx).UpperLimit();

      const double val = rand->Uniform(ll, ul);

      mn_param.SetValue(idx, val);
    }
  }
}
