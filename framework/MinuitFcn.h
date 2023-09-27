#pragma once

#include "Event.h"
#include "Minuit2/FCNBase.h"
#include "Par.h"
#include "SigPDF.h"
#include <iomanip>
#include <iostream>

template<typename SIGPDF>
class MinuitFcn : public ROOT::Minuit2::FCNBase {
private:
  std::vector<Event> &data, &norm;

  const SIGPDF &pdf_sig;

  const double _error_def = 1.0;
  double _best_min = 0.0;

public:
  MinuitFcn(std::vector<Event> &data, std::vector<Event> &norm,
            const SIGPDF &pdf_sig)
          : data(data), norm(norm), pdf_sig(pdf_sig) {
    std::cout << std::setprecision(std::numeric_limits<double>::digits10);
  }

  virtual double operator()(const std::vector<double> &par) const {

    pdf_sig.Norm(norm, par);

    const std::size_t nevents = data.size();

    const double norm_sig = pdf_sig.GetTimeIntegratedNorm(par);

    if (fabs(norm_sig) < 1.0e-15)
      throw std::runtime_error("norm_sig = 0");

    double logL = 0.0;

#pragma omp parallel for reduction(+ : logL)
    for (std::size_t i = 0; i < nevents; ++i) {
      Event &event = data[i];

      // Update amplitudes if necessary
      pdf_sig.updateAmpsInEvent(event, par);

      const double pdf = pdf_sig.GetTimeIntegratedPDF(event, par) / norm_sig;

      if (pdf > 0.0)
        logL += std::log(pdf); //* event.GetSWeight();
    }

    // std::cout << "-2logL = " << (-2.0 * logL) - _best_min << std::endl;

    return (-2.0 * logL) - _best_min;
  }

  virtual double Up() const { return _error_def; }

  [[nodiscard]] double GetBestMin() const { return _best_min; }

  void SetBestMin(const double &best_min) { _best_min += best_min; }
};
