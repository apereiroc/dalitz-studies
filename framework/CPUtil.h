#pragma once

#include "NbodyUtil.h"
#include <array>
#include <iostream>

constexpr unsigned int nCPConf = 2;

/**
   Amplitude ordering
   0: B/D
   1: B-bar/D-bar
 */
enum CPConf { A, Abar };

/**
   Particle ordering for charge conjugate combinations
   A: K+ pi- K- pi+
   Abar: K- pi+ K+ pi-
 */
constexpr std::array<unsigned int, nBody>
CPMap(const unsigned int &CP_conf, std::array<unsigned int, nBody> sequence) {
  switch (CP_conf) {
  case CPConf::A:
    return sequence;
  case CPConf::Abar:
    std::rotate(sequence.begin(), sequence.begin() + 2, sequence.end());
    return sequence;
  default:
    std::cout << "ERROR: Invalid CP_conf = " << CP_conf << " set" << std::endl;
    exit(1);
  }
}
