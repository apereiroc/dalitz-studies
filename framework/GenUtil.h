#pragma once

#include "TGenPhaseSpace.h"
#include "TRandom.h"

/**
   Generate flat event in phase space
 */
inline void generate_flat_event(TGenPhaseSpace &generator) {
  while (true) {
    const double weight = generator.Generate();
    const double height = gRandom->Uniform(1.0 / weight);
    if (height < 1.0)
      break;
  }
}
