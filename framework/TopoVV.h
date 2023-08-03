#pragma once

#include "BarrierFactors.h"
#include "CPUtil.h"
#include "Constant.h"
#include "DalitzUtil.h"
#include "NbodyUtil.h"
#include "qftUtil.h"

/**
   Orbital angular momentum configurations
 */
constexpr unsigned int nVVConf = 3;

enum SpinVV { S, P, D };

/**
   P0 -> V1 -> [P1 P2] V2 -> [P3 P4]
 */
class TopoVV {
public:
  typedef std::array<double, nVVConf> AmpBarrierMap;
  typedef std::array<std::complex<double>, nVVConf> AmpSpinMap;

  double GetAmpBarrier(const unsigned int &L_conf,
                       const unsigned int &CP_conf) const {
    switch (CP_conf) {
    case CPConf::A:
      return _amp_barrier[L_conf];
    case CPConf::Abar:
      return _amp_barrier_CP[L_conf];
    default:
      std::cout << "ERROR: Invalid CP_conf = " << CP_conf << " set"
                << std::endl;
      exit(1);
    }
  }

  void SetAmpsBarrier(const AmpBarrierMap &amp_barrier,
                      const AmpBarrierMap &amp_barrier_CP) {
    _amp_barrier = amp_barrier;
    _amp_barrier_CP = amp_barrier_CP;
  }

  void AmpsBarrier(const std::array<Minimal4Vector, nBody> &p4s,
                   const std::array<Minimal4Vector, nBody> &p4s_CP) {
    _amp_barrier = AmpsBarrier(p4s);
    _amp_barrier_CP = AmpsBarrier(p4s_CP);
  }

  AmpBarrierMap AmpsBarrier(const std::array<Minimal4Vector, nBody> &p4s) {
    return AmpsBarrier(p4s[0], p4s[1], p4s[2], p4s[3]);
  }

  AmpBarrierMap AmpsBarrier(const Minimal4Vector &p1, const Minimal4Vector &p2,
                            const Minimal4Vector &p3,
                            const Minimal4Vector &p4) const {
    const Minimal4Vector pV1 = p1 + p2;
    const Minimal4Vector pV2 = p3 + p4;
    const Minimal4Vector p0 = pV1 + pV2;

    const double qV1 = get_q(pV1.M(), p1.M(), p2.M());
    const double qV2 = get_q(pV2.M(), p3.M(), p4.M());
    const double qP0 = get_q(p0.M(), pV1.M(), pV2.M());

    return {bwbf(qV1, Rrad, 1) * bwbf(qV2, Rrad, 1),
            bwbf(qP0, Srad, 1) * bwbf(qV1, Rrad, 1) * bwbf(qV2, Rrad, 1),
            bwbf(qP0, Srad, 2) * bwbf(qV1, Rrad, 1) * bwbf(qV2, Rrad, 1)};
  }

  std::complex<double> GetAmpSpin(const unsigned int &L_conf,
                                  const unsigned int &CP_conf) const {
    switch (CP_conf) {
    case CPConf::A:
      return _amp_spin[L_conf];
    case CPConf::Abar:
      return _amp_spin_CP[L_conf];
    default:
      std::cout << "ERROR: Invalid CP_conf = " << CP_conf << " set"
                << std::endl;
      exit(1);
    }
  }

  void SetAmpsSpin(const AmpSpinMap &amp_spin, const AmpSpinMap &amp_spin_CP) {
    _amp_spin = amp_spin;
    _amp_spin_CP = amp_spin_CP;
  }

  void AmpsSpin(const std::array<Minimal4Vector, nBody> &p4s,
                const std::array<Minimal4Vector, nBody> &p4s_CP) {
    _amp_spin = AmpsSpin(Minimaltoqft(p4s));
    _amp_spin_CP = AmpsSpin(Minimaltoqft(p4s_CP));
  }

  AmpSpinMap AmpsSpin(const std::array<Vector4<double>, nBody> &p4s) {
    return AmpsSpin(p4s[0], p4s[1], p4s[2], p4s[3]);
  }

  AmpSpinMap AmpsSpin(const Vector4<double> &p1, const Vector4<double> &p2,
                      const Vector4<double> &p3,
                      const Vector4<double> &p4) const {
    const Vector4<double> pV1 = p1 + p2;
    const Vector4<double> pV2 = p3 + p4;

    const Vector4<double> p0 = pV1 + pV2;

    const LeviCivitaTensor lct;

    OrbitalTensor L1_V1(1);
    OrbitalTensor L1_V2(1);

    OrbitalTensor L1_P0(1);
    OrbitalTensor L2_P0(2);

    L1_V1.SetP4(p1, p2);
    L1_V2.SetP4(p3, p4);

    L1_P0.SetP4(pV1, pV2);
    L2_P0.SetP4(pV1, pV2);

    return {(L1_V1 * L1_V2)(0, 0),
            (lct * L1_P0 * L1_V1 * L1_V2 * p0)(0, 0),
            ((L2_P0 * L1_V1) * L1_V2)(0, 0)};
  }

private:
  AmpBarrierMap _amp_barrier, _amp_barrier_CP;
  AmpSpinMap _amp_spin, _amp_spin_CP;
};
