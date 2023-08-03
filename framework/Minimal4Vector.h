#pragma once

#include "TLorentzVector.h"
#include <cmath>

/**
   Minimal 4-vector class designed to save memory and time
 */
class Minimal4Vector {
public:
  Minimal4Vector() = default;

  Minimal4Vector(const TLorentzVector &p4)
      : Minimal4Vector(p4.E(), p4.Px(), p4.Py(), p4.Pz()) {}

  Minimal4Vector(const double &E, const double &px, const double &py,
                 const double &pz)
      : _E(E), _px(px), _py(py), _pz(pz), _m(SetMass(E, px, py, pz)) {}

  Minimal4Vector(const double &px, const double &py, const double &pz)
      : _px(px), _py(py), _pz(pz) {}

  double E() const { return _E; }

  double Px() const { return _px; }

  double Py() const { return _py; }

  double Pz() const { return _pz; }

  TLorentzVector Getp4ROOT() const { return TLorentzVector(_px, _py, _pz, _E); }

  double M() const { return _m; }

  double s() const { return _m * _m; }

  void SetMassAndEnergy(const double &mass) {
    _m = mass;
    _E = std::sqrt(mass * mass + _px * _px + _py * _py + _pz * _pz);
  }

  Minimal4Vector FlipParity() const {
    return Minimal4Vector(_E, -_px, -_py, -_pz);
  }

  Minimal4Vector operator+(const Minimal4Vector &other) const {
    const Minimal4Vector p4(this->_E + other._E, this->_px + other._px,
                            this->_py + other._py, this->_pz + other._pz);

    return p4;
  }

private:
  double SetMass(const double &E, const double &px, const double &py,
                 const double &pz) {
    return sqrt(E * E - px * px - py * py - pz * pz);
  }

  double _E, _px, _py, _pz, _m;
};
