//
// Created by Asier Pereiro Castro on 11/9/23.
//

#pragma once

#include <complex>

#include "qftUtil.h"
#include "CPUtil.h"
#include "Constant.h"
#include "DalitzUtil.h"
#include "BarrierFactors.h"

/**
   Orbital angular momentum configurations
 */
constexpr unsigned int nTopoVSConf = 2;

enum TopoVSConf {
    VS, SV
};

/**
   P0 -> V -> [P1 P2] S -> [P3 P4]
 */
class TopoVS {
public:
    typedef std::array<double, nTopoVSConf> AmpBarrierMap;
    typedef std::array<std::complex<double>, nTopoVSConf> AmpSpinMap;
private:
    AmpBarrierMap _amp_barrier{}, _amp_barrier_CP{};
    AmpSpinMap _amp_spin, _amp_spin_CP;
public:
    [[nodiscard]] double GetAmpBarrier(const TopoVSConf &topo_conf,
                                       const CPConf &CP_conf) const {
        switch (CP_conf) {
            case CPConf::A:
                return _amp_barrier[topo_conf];
            case CPConf::Abar:
                return _amp_barrier_CP[topo_conf];
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

    static AmpBarrierMap AmpsBarrier(
            const std::array<Minimal4Vector, nBody> &p4s) {
        return AmpsBarrier(p4s[0], p4s[1], p4s[2], p4s[3]);
    }

    static AmpBarrierMap
    AmpsBarrier(const Minimal4Vector &p1, const Minimal4Vector &p2,
                const Minimal4Vector &p3, const Minimal4Vector &p4) {
        const Minimal4Vector pV1 = p1 + p2;
        const Minimal4Vector pV2 = p3 + p4;
        const Minimal4Vector p0 = pV1 + pV2;

        const double qV1 = get_q(pV1.M(), p1.M(), p2.M());
        const double qV2 = get_q(pV2.M(), p3.M(), p4.M());
        const double qP0 = get_q(p0.M(), pV1.M(), pV2.M());

        return
                {bwbf(qP0, Srad, 1) * bwbf(qV1, Rrad, 1),
                 bwbf(qP0, Srad, 1) * bwbf(qV2, Rrad, 1)};
    }

    [[nodiscard]] std::complex<double> GetAmpSpin(const TopoVSConf &topo_conf,
                                                  const CPConf &CP_conf) const {
        switch (CP_conf) {
            case CPConf::A:
                return _amp_spin[topo_conf];
            case CPConf::Abar:
                return _amp_spin_CP[topo_conf];
            default:
                std::cout << "ERROR: Invalid CP_conf = " << CP_conf << " set"
                          << std::endl;
                exit(1);
        }
    }

    void SetAmpsSpin(const AmpSpinMap &amp_spin,
                     const AmpSpinMap &amp_spin_CP) {
        _amp_spin = amp_spin;
        _amp_spin_CP = amp_spin_CP;
    }

    void AmpsSpin(const std::array<Minimal4Vector, nBody> &p4s,
                  const std::array<Minimal4Vector, nBody> &p4s_CP) {
        _amp_spin = AmpsSpin(Minimaltoqft(p4s));
        _amp_spin_CP = AmpsSpin(Minimaltoqft(p4s_CP));
    }

    static AmpSpinMap AmpsSpin(const std::array<Vector4<double>, nBody> &p4s) {
        return AmpsSpin(p4s[0], p4s[1], p4s[2], p4s[3]);
    }

    static AmpSpinMap
    AmpsSpin(const Vector4<double> &p1, const Vector4<double> &p2,
             const Vector4<double> &p3, const Vector4<double> &p4) {
        const Vector4<double> pV1 = p1 + p2;
        const Vector4<double> pV2 = p3 + p4;

        const Vector4<double> p0 = pV1 + pV2;

        OrbitalTensor L1_V1(1);
        OrbitalTensor L1_V2(1);

        OrbitalTensor L1_P0(1);

        L1_V1.SetP4(p1, p2);
        L1_V2.SetP4(p3, p4);

        L1_P0.SetP4(pV1, pV2);

        return
                {(L1_P0 * L1_V1)(0, 0),
                 (L1_P0 * L1_V2)(0, 0)};
    }
};
