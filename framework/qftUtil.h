#pragma once

#include "qft++.h"

#include "Minimal4Vector.h"

/**
   Convert 4-vector type between qft++ and ROOT
 */
TLorentzVector qfttoROOT(const Vector4<double> &p4qft) {
  const TLorentzVector p4ROOT(p4qft.Px(), p4qft.Py(), p4qft.Pz(), p4qft.E());

  return p4ROOT;
}

Vector4<double> ROOTtoqft(const TLorentzVector &pROOT) {
  const Vector4<double> pqft(pROOT.E(), pROOT.Px(), pROOT.Py(), pROOT.Pz());

  return pqft;
}

/**
   Convert 4-vector type between Minimal4Vector and qft++
 */
Vector4<double> Minimaltoqft(const Minimal4Vector &p4Min) {
  const Vector4<double> p4qft(p4Min.E(), p4Min.Px(), p4Min.Py(), p4Min.Pz());

  return p4qft;
}

template <std::size_t SIZE>
std::array<Vector4<double>, SIZE>
Minimaltoqft(const std::array<Minimal4Vector, SIZE> &p4Mins) {
  std::array<Vector4<double>, SIZE> p4qfts;
  for (unsigned int i = 0; i < p4Mins.size(); ++i)
    p4qfts[i] = Minimaltoqft(p4Mins[i]);

  return p4qfts;
}
