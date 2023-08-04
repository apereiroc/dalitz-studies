#pragma once

#include "Event.h"
#include "Func.h"
#include "Par.h"

template <typename... AMPS> class SigPDF {
private:
  const std::tuple<AMPS...> _Amps;
  static constexpr unsigned int _nAmps =
      std::tuple_size<std::tuple<AMPS...>>::value;

  const Par idx_tau, idx_DG, idx_Dm;

  double normExpcosh, normExpsinh;

  mutable double _norm_amp2, _norm_astarabarRe;
  mutable double _fit_frac[_nAmps * _nAmps] = {};

  mutable std::complex<double> _c[_nAmps] = {};
  mutable std::complex<double> _cbar[_nAmps] = {};

  mutable double _a2Re[_nAmps * _nAmps] = {};
  mutable double _a2Im[_nAmps * _nAmps] = {};
  mutable double _abar2Re[_nAmps * _nAmps] = {};
  mutable double _abar2Im[_nAmps * _nAmps] = {};
  mutable double _astarabarRe[_nAmps * _nAmps] = {};
  mutable double _astarabarIm[_nAmps * _nAmps] = {};

public:
  SigPDF(const std::tuple<AMPS...> &Amps, const Par &idx_tau, const Par &idx_DG,
         const Par &idx_Dm)
      : _Amps(Amps), idx_tau(idx_tau), idx_DG(idx_DG), idx_Dm(idx_Dm),
        normExpcosh(0.0), normExpsinh(0.0) {}

  inline void NormTime(const std::vector<double> &par) {
    this->normExpcosh =
        norm_expcosh(time_min, time_max, par[idx_tau], par[idx_DG]);
    this->normExpsinh =
        norm_expsinh(time_min, time_max, par[idx_tau], par[idx_DG]);
  }

  double GetPDF(const Event &event, const double &time, const int &qtag,
                const std::vector<double> &par) const {
    const double q = static_cast<double>(qtag);

    const std::complex<double> A = GetAmp(event, par, CPConf::A);
    const std::complex<double> Abar = GetAmp(event, par, CPConf::Abar);

    return std::exp(-time / par[idx_tau]) *
           (((std::norm(A) + std::norm(Abar)) *
             std::cosh(par[idx_DG] / 2.0 * time)) -
            (2.0 * std::real(std::conj(A) * Abar) *
             std::sinh(par[idx_DG] / 2.0 * time)) +
            (q * (std::norm(A) - std::norm(Abar)) *
             std::cos(par[idx_Dm] * time)) -
            (2.0 * q * std::imag(std::conj(A) * Abar) *
             std::sin(par[idx_Dm] * time)));
  }

  double GetTimeIntegratedPDF(const Event &event,
                              const std::vector<double> &par) const {
    const std::complex<double> A = GetAmp(event, par, CPConf::A);
    const std::complex<double> Abar = GetAmp(event, par, CPConf::Abar);

    const double pdf_amp2 = std::norm(A) + std::norm(Abar);
    const double pdf_astarabar = std::real(std::conj(A) * Abar);

    return pdf_amp2 * normExpcosh - 2.0 * pdf_astarabar * normExpsinh;
  }

  double GetTimeIntegratedNorm(const std::vector<double> &par) const {
    return this->GetNormAmp2() * normExpcosh -
           2.0 * this->GetNormReAstarAbar() * normExpsinh;
    (void)par;
  }

  /**
   * Returns the total amplitude A = a_S * S + a_P * P + ...
   */
  std::complex<double> GetAmp(const Event &event,
                              const std::vector<double> &par,
                              const CPConf &CP_conf) const {
    std::complex<double> amp = {0.0, 0.0};
    const auto sum_amps =
        [&amp, &_event = std::as_const(event), &_par = std::as_const(par),
         &_CP_conf = std::as_const(CP_conf)](const auto &...Amps) {
          std::size_t i = 0;
          ((amp +=
            Amps.GetCoupling(_par, _CP_conf) * _event.GetAmps(_CP_conf)(i++)),
           ...);
        };

    std::apply(sum_amps, _Amps);

    return amp;
  }

  void Amps(Event &event, const std::vector<double> &par) const {
    const auto lambda_amps = [&_event = event,
                              &_par = std::as_const(par)](const auto &...Amps) {
      std::size_t i = 0;
      const auto lambda_amp = [&_event = _event, &_par = std::as_const(_par),
                               &_amp_idx = i](const auto &Amp) {
        if (Amp.ChangePar() == false) {
          _amp_idx++;
          return;
        }

        const std::complex<double> a = Amp.GetAmp(_event, _par, CPConf::A);
        const std::complex<double> abar =
            Amp.GetAmp(_event, _par, CPConf::Abar);

        _event.GetToSetAmps(CPConf::A)(_amp_idx) = a;
        _event.GetToSetAmps(CPConf::Abar)(_amp_idx) = abar;

        _amp_idx++;
      };

      (lambda_amp(Amps), ...);
    };
    std::apply(lambda_amps, _Amps);
  }

  /**
     Normalise PDF
   */
  void UpdateCouplings(const std::vector<double> &par) const {
    std::apply(
        [&_c = _c, &_cbar = _cbar,
         &_par = std::as_const(par)](const auto &...Amps) {
          std::size_t i = 0;
          const auto lambda_amp = [&_c, &_cbar, &_par = std::as_const(_par),
                                   &_amp_idx = i](const auto &Amp) {
            _c[_amp_idx] = Amp.GetCoupling(_par, CPConf::A);
            _cbar[_amp_idx] = Amp.GetCoupling(_par, CPConf::Abar);

            _amp_idx++;
          };

          (lambda_amp(Amps), ...);
        },
        _Amps);
  }

  double GetNormAmp2() const { return _norm_amp2; }

  double GetNormReAstarAbar() const { return _norm_astarabarRe; }

  void Norm(std::vector<Event> &norm, const std::vector<double> &par) const {
    // Determine which amplitudes should be updated
    bool par_change = false;
    std::apply(
        [&_par = std::as_const(par),
         &_par_change = par_change](const auto &...Amps) {
          const auto lambda_amp = [&_par = std::as_const(_par),
                                   &_par_change =
                                       _par_change](const auto &Amp) {
            Amp.UpdatePar(_par);
            if (Amp.ChangePar() == true)
              _par_change = true;
          };

          (lambda_amp(Amps), ...);
        },
        _Amps);

    if (par_change == true)
      NormAmp2(norm, par);

    UpdateCouplings(par);

    double norm_amp2 = 0.0, norm_astarabarRe = 0.0;
#pragma omp parallel for reduction(+ : norm_amp2, norm_astarabarRe)
    for (unsigned int i = 0; i < _nAmps; ++i)
      for (unsigned int j = 0; j < _nAmps; ++j) {
        _fit_frac[i * _nAmps + j] =
            std::real(_c[i] * std::conj(_c[j])) * _a2Re[i * _nAmps + j] -
            std::imag(_c[i] * std::conj(_c[j])) * _a2Im[i * _nAmps + j] +
            std::real(_cbar[i] * std::conj(_cbar[j])) *
                _abar2Re[i * _nAmps + j] -
            std::imag(_cbar[i] * std::conj(_cbar[j])) *
                _abar2Im[i * _nAmps + j];
        norm_amp2 += _fit_frac[i * _nAmps + j];

        norm_astarabarRe += std::real(std::conj(_c[i]) * _cbar[j]) *
                                _astarabarRe[i * _nAmps + j] -
                            std::imag(std::conj(_c[i]) * _cbar[j]) *
                                _astarabarIm[i * _nAmps + j];
      }

    _norm_amp2 = norm_amp2;
    _norm_astarabarRe = norm_astarabarRe;
  }

  void NormAmp2(std::vector<Event> &norm,
                const std::vector<double> &par) const {
    std::fill(_a2Re, _a2Re + _nAmps * _nAmps, 0.0);
    std::fill(_a2Im, _a2Im + _nAmps * _nAmps, 0.0);
    std::fill(_abar2Re, _abar2Re + _nAmps * _nAmps, 0.0);
    std::fill(_abar2Im, _abar2Im + _nAmps * _nAmps, 0.0);
    std::fill(_astarabarRe, _astarabarRe + _nAmps * _nAmps, 0.0);
    std::fill(_astarabarIm, _astarabarIm + _nAmps * _nAmps, 0.0);

    const unsigned int nEvents = norm.size();

#pragma omp parallel for reduction(                                            \
        + : _a2Re[ : _nAmps * _nAmps], _a2Im[ : _nAmps * _nAmps],              \
            _abar2Re[ : _nAmps * _nAmps], _abar2Im[ : _nAmps * _nAmps],        \
            _astarabarRe[ : _nAmps * _nAmps],                                  \
            _astarabarIm[ : _nAmps * _nAmps])

    for (unsigned int i = 0; i < nEvents; ++i) {
      Event &event = norm[i];

      std::apply(
          [&_event = event, &_par = std::as_const(par)](const auto &...Amps) {
            std::size_t j = 0;
            const auto lambda_amp = [&_event = _event,
                                     &_par = std::as_const(_par),
                                     &_amp_idx = j](const auto &Amp) {
              if (Amp.ChangePar() == false) {
                _amp_idx++;
                return;
              }

              const std::complex<double> a =
                  sqrt(_event.GetEff()) * Amp.GetAmp(_event, _par, CPConf::A);
              const std::complex<double> abar =
                  sqrt(_event.GetEff()) *
                  Amp.GetAmp(_event, _par, CPConf::Abar);

              _event.GetToSetAmps(CPConf::A)(_amp_idx) = a;
              _event.GetToSetAmps(CPConf::Abar)(_amp_idx) = abar;
              _amp_idx++;
            };

            (lambda_amp(Amps), ...);
          },
          _Amps);

      for (unsigned int j = 0; j < _nAmps; ++j)
        for (unsigned int k = 0; k < _nAmps; ++k) {
          const std::complex<double> a2 =
              event.GetAmps(CPConf::A)(j) *
              std::conj(event.GetAmps(CPConf::A)(k));
          const std::complex<double> abar2 =
              event.GetAmps(CPConf::Abar)(j) *
              std::conj(event.GetAmps(CPConf::Abar)(k));

          _a2Re[j * _nAmps + k] += std::real(a2) / event.GetGenPDF();
          _a2Im[j * _nAmps + k] += std::imag(a2) / event.GetGenPDF();
          _abar2Re[j * _nAmps + k] += std::real(abar2) / event.GetGenPDF();
          _abar2Im[j * _nAmps + k] += std::imag(abar2) / event.GetGenPDF();

          const std::complex<double> astarabar =
              std::conj(event.GetAmps(CPConf::A)(j)) *
              event.GetAmps(CPConf::Abar)(k);

          _astarabarRe[j * _nAmps + k] +=
              std::real(astarabar) / event.GetGenPDF();
          _astarabarIm[j * _nAmps + k] +=
              std::imag(astarabar) / event.GetGenPDF();
        }
    }
  }

  unsigned int GetNamps() const { return _nAmps; }

  /**
     Return fit fractions calculated from normalisation
     for a given amplitude index
   */
  const double &GetFitFrac(const unsigned int &idx_amp) const {
    return _fit_frac[idx_amp];
  }

  /**
     Resize amplitude container in Event class to match
     number of amplitudes
   */
  void ResizeEvents(Event &event) {
    event.GetToSetAmps(CPConf::A).resize(_nAmps);
    event.GetToSetAmps(CPConf::Abar).resize(_nAmps);
  }

  void ResizeEvents(std::vector<Event> &data) {
    const unsigned int nEvents = data.size();

#pragma omp parallel for
    for (unsigned int i = 0; i < nEvents; ++i)
      ResizeEvents(data[i]);
  }

  void ResizeEvents(std::vector<Event> &data, std::vector<Event> &norm) {
    ResizeEvents(data);
    ResizeEvents(norm);
  }
};
