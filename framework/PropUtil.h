#pragma once

#include "MassModels.h"
#include "Propagator.h"
#include <iostream>
#include <memory>

enum PropConf { BW };

inline std::shared_ptr<Propagator> NewPropagator(const PropConf &conf) {
  switch (conf) {
  case PropConf::BW:
    return std::make_shared<class BreitWigner>();
  default: {
    std::cout << "ERROR: propagator not implemented" << std::endl;
    std::exit(1);
  }
  }
}
