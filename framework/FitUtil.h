#pragma once

#include "FitParameters.h"
#include "Par.h"
#include "TRandom3.h"

inline void RandomiseCovariantCoeffs(FitParameters &mn_param,
                                     const unsigned int &seed) {
  const auto rand = std::make_unique<TRandom3>(seed);

  for (unsigned int idx = Par::x_VV_S; idx <= Par::y_SS; idx++) {
    if (!mn_param.Parameter(idx).IsFixed()) {
      const double ll = mn_param.Parameter(idx).LowerLimit();
      const double ul = mn_param.Parameter(idx).UpperLimit();

      const double val = rand->Uniform(ll, ul);

      mn_param.SetValue(idx, val);
    }
  }
}

inline void print_minimum_info(const ROOT::Minuit2::FunctionMinimum &min) {

  enum MinType { YesIsGood, YesIsBad };

  const auto print_info = [](const std::string &message, const bool &response,
                             const MinType &type) {
    switch (type) {
    case YesIsGood:
      if (response == true)
        std::cout << message << "yes ✅\n";
      else
        std::cout << message << "no  ❗️\n";
      break;
    case YesIsBad:
      if (response == true)
        std::cout << message << "yes ❗️\n";
      else
        std::cout << message << "no  ✅\n";
    }
  };

  print_info("    Is valid:            ", min.IsValid(), MinType::YesIsGood);
  print_info("    HasValidParameters:  ", min.HasValidParameters(),
             MinType::YesIsGood);
  print_info("    HasValidCovariance:  ", min.HasValidCovariance(),
             MinType::YesIsGood);
  print_info("    HasAccurateCovar:    ", min.HasAccurateCovar(),
             MinType::YesIsGood);
  print_info("    HasPosDefCovar:      ", min.HasPosDefCovar(),
             MinType::YesIsGood);
  print_info("    HasMadePosDefCovar:  ", min.HasMadePosDefCovar(),
             MinType::YesIsBad);
  print_info("    HesseFailed:         ", min.HesseFailed(), MinType::YesIsBad);
  print_info("    HasCovariance:       ", min.HasCovariance(),
             MinType::YesIsGood);
  print_info("    IsAboveMaxEdm:       ", min.IsAboveMaxEdm(),
             MinType::YesIsBad);
  print_info("    HasReachedCallLimit: ", min.HasReachedCallLimit(),
             MinType::YesIsBad);
}
