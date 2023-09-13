#pragma once

#include "MassModels.h"
#include "Propagator.h"
#include <iostream>
#include <memory>

enum PropConf { BW, LASS };

inline std::shared_ptr<Propagator> NewPropagator(const PropConf &conf) {
  switch (conf) {
  case PropConf::BW:
    return std::make_shared<class BreitWigner>();
  case PropConf::LASS:
      return std::make_shared<class LASS>();
  default: {
    std::cout << "ERROR: propagator not implemented" << std::endl;
    std::exit(1);
  }
  }
}
