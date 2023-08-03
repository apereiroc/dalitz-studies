#pragma once

#include "Propagator.h"
#include <memory>

enum PropConf { BW };

inline std::shared_ptr<Propagator> NewPropagator(const PropConf &conf);
