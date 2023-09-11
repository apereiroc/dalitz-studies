#pragma once

#include <array>

#include "CPUtil.h"
#include "Minimal4Vector.h"
#include "NbodyUtil.h"
#include "TLorentzVector.h"
#include "TopoVV.h"
#include <Eigen/Core>

class Event {
private:
  double _gen_pdf = 1.0;

  double _time;
  int _qtag;

  std::array<Minimal4Vector, nBody> _p4{}, _p4_CP{};
  TopoVV VV;
  Eigen::VectorXcd _a, _abar;
  double _theta1{}, _theta2{}, _chi{};
  double _mKpPim, _mKmPip;

  double _eff = 1.0;

public:
  /**
     Set has_cache to false if barrier factors and spin amplitudes not already
     calculated
  */
  Event(const TLorentzVector &p4Kp, const TLorentzVector &p4pim,
        const TLorentzVector &p4Km, const TLorentzVector &p4pip,
        const double &time, const int &qtag, const bool &has_cache = false)
      : _time(time), _qtag(qtag) {
    _p4[0] = Minimal4Vector(p4Kp);
    _p4[1] = Minimal4Vector(p4pim);
    _p4[2] = Minimal4Vector(p4Km);
    _p4[3] = Minimal4Vector(p4pip);

    for (unsigned int i = 0; i < nBody; ++i)
      _p4_CP[CPMap(CPConf::Abar, ordering)[i]] = _p4[i].FlipParity();

    if (has_cache == false) {
      VV.AmpsBarrier(_p4, _p4_CP);
      VV.AmpsSpin(_p4, _p4_CP);
    }

    Helicities();

    _mKpPim = (p4Kp + p4pim).M();
    _mKmPip = (p4Km + p4pip).M();
  }

  Event(const double &Kp_PX, const double &Kp_PY, const double &Kp_PZ,
        const double &Kp_E, const double &pim_PX, const double &pim_PY,
        const double &pim_PZ, const double &pim_E, const double &Km_PX,
        const double &Km_PY, const double &Km_PZ, const double &Km_E,
        const double &pip_PX, const double &pip_PY, const double &pip_PZ,
        const double &pip_E, const double &time, const int &qtag,
        const bool &has_cache = false)
      : _time(time), _qtag(qtag) {
    _p4[0] = Minimal4Vector(Kp_E, Kp_PX, Kp_PY, Kp_PZ);
    _p4[1] = Minimal4Vector(pim_E, pim_PX, pim_PY, pim_PZ);
    _p4[2] = Minimal4Vector(Km_E, Km_PX, Km_PY, Km_PZ);
    _p4[3] = Minimal4Vector(pip_E, pip_PX, pip_PY, pip_PZ);

    for (unsigned int i = 0; i < nBody; ++i)
      _p4_CP[CPMap(CPConf::Abar, ordering)[i]] = _p4[i].FlipParity();

    if (has_cache == false) {
      VV.AmpsBarrier(_p4, _p4_CP);
      VV.AmpsSpin(_p4, _p4_CP);
    }

    Helicities();
    _mKpPim = (_p4[0] + _p4[1]).M();
    _mKmPip = (_p4[2] + _p4[3]).M();
  }

  Event(const double &Kp_PX, const double &Kp_PY, const double &Kp_PZ,
        const double &pim_PX, const double &pim_PY, const double &pim_PZ,
        const double &Km_PX, const double &Km_PY, const double &Km_PZ,
        const double &pip_PX, const double &pip_PY, const double &pip_PZ,
        const double &time, const int &qtag, const bool &has_cache = false)
      : _time(time), _qtag(qtag) {
    _p4[0] = Minimal4Vector(Kp_PX, Kp_PY, Kp_PZ);
    _p4[1] = Minimal4Vector(pim_PX, pim_PY, pim_PZ);
    _p4[2] = Minimal4Vector(Km_PX, Km_PY, Km_PZ);
    _p4[3] = Minimal4Vector(pip_PX, pip_PY, pip_PZ);

    _p4[0].SetMassAndEnergy(mass_Kp);
    _p4[1].SetMassAndEnergy(mass_pip);
    _p4[2].SetMassAndEnergy(mass_Kp);
    _p4[3].SetMassAndEnergy(mass_pip);

    for (unsigned int i = 0; i < nBody; ++i)
      _p4_CP[CPMap(CPConf::Abar, ordering)[i]] = _p4[i].FlipParity();

    if (has_cache == false) {
      VV.AmpsBarrier(_p4, _p4_CP);
      VV.AmpsSpin(_p4, _p4_CP);
    }

    Helicities();
    _mKpPim = (_p4[0] + _p4[1]).M();
    _mKmPip = (_p4[2] + _p4[3]).M();
  }

  [[nodiscard]] const double &GetGenPDF() const { return _gen_pdf; }

  void SetGenPDF(const double &gen_pdf) { _gen_pdf = gen_pdf; }

  [[nodiscard]] const double &GetTime() const { return _time; }

  void SetTime(const double &time) { _time = time; }

  [[nodiscard]] const int &GetQtag() const { return _qtag; }

  void SetQtag(const int &qtag) { _qtag = qtag; }

  [[nodiscard]] const std::array<Minimal4Vector, nBody> &
  Getp4(const unsigned int &CP_conf) const {
    switch (CP_conf) {
    case CPConf::A:
      return _p4;
    case CPConf::Abar:
      return _p4_CP;
    default:
      std::cout << "ERROR: Invalid CP_conf = " << CP_conf << " set"
                << std::endl;
      exit(1);
    }
  }

  [[nodiscard]] const Minimal4Vector &Getp4Kp(const unsigned int &CP_conf) const {
    switch (CP_conf) {
    case CPConf::A:
      return _p4[0];
    case CPConf::Abar:
      return _p4_CP[0];
    default:
      std::cout << "ERROR: Invalid CP_conf = " << CP_conf << " set"
                << std::endl;
      exit(1);
    }
  }

  [[nodiscard]] const Minimal4Vector &Getp4pim(const unsigned int &CP_conf) const {
    switch (CP_conf) {
    case CPConf::A:
      return _p4[1];
    case CPConf::Abar:
      return _p4_CP[1];
    default:
      std::cout << "ERROR: Invalid CP_conf = " << CP_conf << " set"
                << std::endl;
      exit(1);
    }
  }

  [[nodiscard]] const Minimal4Vector &Getp4Km(const unsigned int &CP_conf) const {
    switch (CP_conf) {
    case CPConf::A:
      return _p4[2];
    case CPConf::Abar:
      return _p4_CP[2];
    default:
      std::cout << "ERROR: Invalid CP_conf = " << CP_conf << " set"
                << std::endl;
      exit(1);
    }
  }

  [[nodiscard]] const Minimal4Vector &Getp4pip(const unsigned int &CP_conf) const {
    switch (CP_conf) {
    case CPConf::A:
      return _p4[3];
    case CPConf::Abar:
      return _p4_CP[3];
    default:
      std::cout << "ERROR: Invalid CP_conf = " << CP_conf << " set"
                << std::endl;
      exit(1);
    }
  }

  TopoVV &GetToSetVV() { return VV; }

  [[nodiscard]] const TopoVV &GetVV() const { return VV; }

  [[nodiscard]] const Eigen::VectorXcd &GetAmps(const unsigned int &CP_conf) const {
    switch (CP_conf) {
    case CPConf::A:
      return _a;
    case CPConf::Abar:
      return _abar;
    default:
      exit(1);
    }
  }

  Eigen::VectorXcd &GetToSetAmps(const unsigned int &CP_conf) {
    switch (CP_conf) {
    case CPConf::A:
      return _a;
    case CPConf::Abar:
      return _abar;
    default:
      exit(1);
    }
  }

  [[nodiscard]] double GetTheta1() const { return _theta1; }

  [[nodiscard]] double GetTheta2() const { return _theta2; }

  [[nodiscard]] double GetChi() const { return _chi; }

  [[nodiscard]] double GetCosTheta1() const { return cos(_theta1); }

  [[nodiscard]] double GetCosTheta2() const { return cos(_theta2); }

  [[nodiscard]] double GetSinTheta1() const { return sin(_theta1); }

  [[nodiscard]] double GetSinTheta2() const { return sin(_theta2); }

  [[nodiscard]] double GetSin2Theta1() const { return sin(2.0 * _theta1); }

  [[nodiscard]] double GetSin2Theta2() const { return sin(2.0 * _theta2); }

  [[nodiscard]] double GetSinChi() const { return sin(_chi); }

  [[nodiscard]] double GetCosChi() const { return cos(_chi); }

  [[nodiscard]] double GetSin2Chi() const { return sin(2.0 * _chi); }

  [[nodiscard]] double GetMassKpPim() const { return (_p4[0] + _p4[1]).M(); }

  [[nodiscard]] double GetMassKmPip() const { return (_p4[2] + _p4[3]).M(); }

  void Helicities() {
    const TLorentzVector p4Kp = _p4[0].Getp4ROOT();
    const TLorentzVector p4pim = _p4[1].Getp4ROOT();
    const TLorentzVector p4Km = _p4[2].Getp4ROOT();
    const TLorentzVector p4pip = _p4[3].Getp4ROOT();

    const TVector3 V1_boost = (p4Kp + p4pim).BoostVector();
    const TVector3 V2_boost = (p4Km + p4pip).BoostVector();
    const TVector3 B_boost = (p4Kp + p4Km + p4pip + p4pim).BoostVector();

    TLorentzVector pKp_V1(p4Kp);
    TLorentzVector pV2_V1(p4Km + p4pip);
    TLorentzVector pKm_V2(p4Km);
    TLorentzVector pV1_V2(p4Kp + p4pim);

    pKp_V1.Boost(-V1_boost);
    pV2_V1.Boost(-V1_boost);
    pKm_V2.Boost(-V2_boost);
    pV1_V2.Boost(-V2_boost);

    _theta1 = atan2((pKp_V1.Vect().Cross(-pV2_V1.Vect())).Mag(),
                    pKp_V1.Vect().Dot(-pV2_V1.Vect()));
    _theta2 = atan2((pKm_V2.Vect().Cross(-pV1_V2.Vect())).Mag(),
                    pKm_V2.Vect().Dot(-pV1_V2.Vect()));

    TLorentzVector p4Kp_B(p4Kp);
    TLorentzVector p4Km_B(p4Km);
    TLorentzVector p4pip_B(p4pip);
    TLorentzVector p4pim_B(p4pim);

    p4Kp_B.Boost(-B_boost);
    p4Km_B.Boost(-B_boost);
    p4pim_B.Boost(-B_boost);
    p4pip_B.Boost(-B_boost);

    const TVector3 V1_norm = p4Kp_B.Vect().Cross(p4pim_B.Vect());
    const TVector3 V2_norm = p4Km_B.Vect().Cross(p4pip_B.Vect());
    const TVector3 chi_norm =
        (p4Km_B + p4pip_B).Vect() * (1.0 / (p4Km_B + p4pip_B).Vect().Mag());

    _chi = atan2((V1_norm.Cross(V2_norm)).Dot(chi_norm), V1_norm.Dot(V2_norm));
  }

  [[nodiscard]] const double &GetEff() const { return _eff; }

  void SetEff(const double &eff) { _eff = eff; }

  void WeightEff(const double &x) { _eff *= x; }
};
