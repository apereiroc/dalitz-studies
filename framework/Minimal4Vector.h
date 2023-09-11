#pragma once

#include "TLorentzVector.h"
#include <cmath>

/**
   Minimal 4-vector class designed to save memory and time
 */
class Minimal4Vector {
private:
    static double SetMass(const double &E, const double &px, const double &py,
                          const double &pz) {
        return sqrt(E * E - px * px - py * py - pz * pz);
    }

    double _px, _py, _pz, _m, _E;
public:
  Minimal4Vector() = default;

  explicit Minimal4Vector(const TLorentzVector &p4)
      : Minimal4Vector(p4.E(), p4.Px(), p4.Py(), p4.Pz()) {}

  Minimal4Vector(const double &E, const double &px, const double &py,
                 const double &pz)
      : _px(px), _py(py), _pz(pz), _m(SetMass(E, px, py, pz)), _E(E) {}

  Minimal4Vector(const double &px, const double &py, const double &pz)
      : _px(px), _py(py), _pz(pz), _m(0), _E(P()) {}

  [[nodiscard]] double E() const { return _E; }

  [[nodiscard]] double P() const {
      return std::sqrt(_px * _px + _py * _py + _pz * _pz);
  }

  [[nodiscard]] double Px() const { return _px; }

  [[nodiscard]] double Py() const { return _py; }

  [[nodiscard]] double Pz() const { return _pz; }

  [[nodiscard]] TLorentzVector Getp4ROOT() const { return {_px, _py, _pz, _E}; }

  [[nodiscard]] double M() const { return _m; }

  [[nodiscard]] double s() const { return _m * _m; }

  void SetMassAndEnergy(const double &mass) {
    _m = mass;
      _E = std::sqrt(mass * mass + _px * _px + _py * _py + _pz * _pz);
  }

  [[nodiscard]] Minimal4Vector FlipParity() const {
    return {_E, -_px, -_py, -_pz};
  }

  Minimal4Vector operator+(const Minimal4Vector &other) const {
    const Minimal4Vector p4(this->_E + other._E, this->_px + other._px,
                            this->_py + other._py, this->_pz + other._pz);

    return p4;
  }


};
