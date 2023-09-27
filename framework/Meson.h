//
// Created by Asier Pereiro Castro on 27/9/23.
//

#pragma once

#include <iostream>
#include "Par.h"
#include "Constant.h"

enum class MesonConf {
  Bs,
  Du,
};

inline MesonConf parseMeson(const std::string &meson) {
  if (meson == "Bs")
    return MesonConf::Bs;
  else if (meson == "Du")
    return MesonConf::Du;
  else {
    std::cout << "ERROR: could not parse meson\n";
    std::exit(1);
  }
}

class Meson {
private:
  Par idx_tau, idx_DG, idx_Dm;
  double mass;

public:
  Meson(const MesonConf &conf) {
    switch (conf) {
      case MesonConf::Bs:
        idx_tau = Par::tau_Bs;
        idx_DG = Par::DG_Bs;
        idx_Dm = Par::Dm_Bs;
        mass = mass_Bs;
        break;
      case MesonConf::Du:
        idx_tau = Par::tau_Du;
        idx_DG = Par::DG_Du;
        idx_Dm = Par::Dm_Du;
        mass = mass_Du;
        break;
      default:
        std::cout << "ERROR: wrong meson\n";
        std::exit(1);
    }
  }

  Meson(const std::string &meson) :
          Meson(parseMeson(meson)) {}

  inline Par getTauIdx() const { return idx_tau; }

  inline Par getDGIdx() const { return idx_DG; }

  inline Par getDmIdx() const { return idx_Dm; }

  inline double getMass() const { return mass; }
};
