#pragma once

#include <cmath>
#include <iostream>
#include <vector>

#include <Eigen/Dense>

#include "EigenUtil.h"
#include "FitCovariance.h"
#include "FitParameters.h"

template <typename PDF, typename EVENT> class Results {
protected:
  const PDF &_pdf;
  std::vector<EVENT> &_norm;

private:
  const FitParameters &_mn_param;

  bool _cov_loaded = false;
  FitCovariance _cov;
  std::vector<FitParameters> _mn_param_gen;

  Eigen::MatrixXd _ff, _ff_err;
  double _ffsum, _ffsum_err;

  double _ff_vv_SD = 0.0;
  double _ff_vv_P = 0.0;
  double _ff_vs_plus = 0.0;
  double _ff_vs_minus = 0.0;
  double _ff_ss = 0.0;
  double _ff_scalar = 0.0;

  double _ff_vv_SD_err = 0.0;
  double _ff_vv_P_err = 0.0;
  double _ff_vs_plus_err = 0.0;
  double _ff_vs_minus_err = 0.0;
  double _ff_ss_err = 0.0;
  double _ff_scalar_err = 0.0;

  Eigen::VectorXd _phases, _phases_err;

  double _arg_vv_D = 0.0;
  double _arg_vs_plus = 0.0;
  double _arg_vs_minus = 0.0;
  double _arg_ss = 0.0;

  double _arg_vv_D_err = 0.0;
  double _arg_vs_plus_err = 0.0;
  double _arg_vs_minus_err = 0.0;
  double _arg_ss_err = 0.0;

public:
  Results(const PDF &pdf, std::vector<EVENT> &norm,
          const FitParameters &mn_param)
      : _pdf(pdf), _norm(norm), _mn_param(mn_param) {}

  /**
     Load covariance matrix and generate pseudoexperiments of the free
     parameters
   */
  void LoadCov(const std::string &filename, const unsigned int nexp = 10000) {
    if (_cov_loaded == true) {
      std::cout << "ERROR: Covariance matrix can only be loaded once"
                << std::endl;
      exit(1);
    }

    FitCovariance cov(filename);
    _cov = cov;

    _mn_param_gen = GenParam(_mn_param, _cov, nexp);
    std::cout << nexp << " pseudoexperiments generated" << std::endl;

    _cov_loaded = true;
  }

  void LoadCov(const FitCovariance &cov, const unsigned int nexp = 10000) {
    if (_cov_loaded == true) {
      std::cout << "ERROR: Covariance matrix can only be loaded once"
                << std::endl;
      exit(1);
    }

    _cov = cov;

    _mn_param_gen = GenParam(_mn_param, _cov, nexp);
    std::cout << nexp << " pseudoexperiments generated" << std::endl;

    _cov_loaded = true;
  }

  template <typename FUNCTION> Eigen::VectorXd GetObsErr(const FUNCTION &func) {
    const Eigen::VectorXd obs = func(_mn_param);

    Eigen::VectorXd obs_err(obs.size());
    obs_err.setZero();
    for (unsigned int i = 0; i < _mn_param_gen.size(); ++i) {
      const Eigen::VectorXd obs_gen = func(_mn_param_gen[i]);
      const Eigen::VectorXd obs_diff2 = (obs_gen - obs).array().pow(2.0);
      obs_err += obs_diff2;
    }

    return obs_err.array().sqrt() /
           sqrt(static_cast<double>(_mn_param_gen.size()));
  }

  Eigen::MatrixXd GetFitFrac() const { return _ff; }

  double GetFitFracSum() const { return _ffsum; }

  void FitFrac() {
    _ff = FitFrac(_mn_param);
    _ffsum = _ff.trace();

    const auto _ff_sumcols = _ff.colwise().sum().array();
    _ff_vv_SD = _ff_sumcols(0) + _ff_sumcols(2);
    _ff_vv_P = _ff_sumcols(1);
    _ff_vs_plus = _ff_sumcols(3);
    _ff_vs_minus = _ff_sumcols(4);
    _ff_ss = _ff_sumcols(5);
    _ff_scalar = _ff_vs_plus + _ff_vs_minus + _ff_ss;
  }

  void Phases() {
    _phases = Phases(_mn_param);

    // _arg_vv_S = _phases(0);
    // _arg_vv_P = _phases(1);
    _arg_vv_D = _phases(2);
    _arg_vs_plus = _phases(3);
    _arg_vs_minus = _phases(4);
    _arg_ss = _phases(5);
  }

  Eigen::MatrixXd FitFrac(const FitParameters &mn_param) {
    const unsigned int nAmps = _pdf.GetNamps();

    const std::vector<double> &par = mn_param.Params();

    _pdf.Norm(_norm, par);

    // Exclude B-Bbar interference. Don't need it for fit fractions
    const double norm_sig = _pdf.GetNormAmp2();

    Eigen::MatrixXd fit_frac(nAmps, nAmps);
    fit_frac.setZero();

    // Double count the number of off-diagonal elements
    for (unsigned int i = 0; i < nAmps; ++i)
      for (unsigned int j = i; j < nAmps; ++j) {
        fit_frac(i, j) = 100.0 * _pdf.GetFitFrac(i * nAmps + j) / norm_sig;

        if (i != j)
          fit_frac(i, j) *= 2.0;
      }

    return fit_frac;
  }

  Eigen::VectorXd Phases(const FitParameters &mn_param) {
    const std::vector<double> &par = mn_param.Params();
    _pdf.ComputePhases(par);

    const unsigned int nAmps = _pdf.GetNamps();
    Eigen::VectorXd phases(nAmps);
    phases.setZero();

    for (unsigned int i = 0; i < nAmps; ++i) {
      phases(i) = _pdf.GetPhase(i);
    }

    return phases;
  }

  Eigen::MatrixXd GetFitFracErr() const { return _ff_err; }

  double GetFitFracSumErr() const { return _ffsum_err; }

  /* double GetFFVector() const { return _ff_vector; } */
  /**/
  /* double GetFFVectorErr() const { return _ff_vector_err; } */
  /**/
  double GetFF_Scalar() const { return _ff_scalar; }

  double GetFF_Scalar_Err() const { return _ff_scalar_err; }

  double GetFF_VV_SD() const { return _ff_vv_SD; }
  double GetFF_VV_SD_Err() const { return _ff_vv_SD_err; }

  double GetFF_VV_P() const { return _ff_vv_P; }
  double GetFF_VV_P_Err() const { return _ff_vv_P_err; }

  double GetFF_VS_Plus() const { return _ff_vs_plus; }
  double GetFF_VS_Plus_Err() const { return _ff_vs_plus_err; }

  double GetFF_VS_Minus() const { return _ff_vs_minus; }
  double GetFF_VS_Minus_Err() const { return _ff_vs_minus_err; }

  double GetFF_SS() const { return _ff_ss; }
  double GetFF_SS_Err() const { return _ff_ss_err; }

  double GetArg_VV_D() const { return _arg_vv_D; }
  double GetArg_VV_D_Err() const { return _arg_vv_D_err; }

  double GetArg_VS_Plus() const { return _arg_vs_plus; }
  double GetArg_VS_Plus_Err() const { return _arg_vs_plus_err; }

  double GetArg_VS_Minus() const { return _arg_vs_minus; }
  double GetArg_VS_Minus_Err() const { return _arg_vs_minus_err; }

  double GetArg_SS() const { return _arg_ss; }
  double GetArg_SS_Err() const { return _arg_ss_err; }

  void FitFracErr(const unsigned int &nexp = 1000) {
    if (nexp > _mn_param_gen.size()) {
      std::cout << "ERROR: Number of pseudoexperiments exceeds that generated "
                << "from covariance matrix" << std::endl;
      exit(1);
    }

    // Fit fractions
    _ff_err.resize(_ff.rows(), _ff.cols());
    _ff_err.setZero();
    _ffsum_err = 0.0;

    _ff_vv_SD_err = 0.0;
    _ff_vv_P_err = 0.0;
    _ff_vs_plus_err = 0.0;
    _ff_vs_minus_err = 0.0;
    _ff_ss_err = 0.0;
    _ff_scalar_err = 0.0;

    // Phases
    // _phases_err.resize(_phases.size());
    // _phases_err.setZero();

    _arg_vv_D_err = 0.0;
    _arg_vs_plus_err = 0.0;
    _arg_vs_minus_err = 0.0;
    _arg_ss_err = 0.0;

    // Parallel reductions have been already explored but
    // no substantial gains are observed, since the PDF
    // is being normalised in parallel as well
    for (unsigned int i = 0; i < nexp; ++i) {
      std::cout << "Fit Fraction Exp: " << i << '\n';

      const FitParameters &mn_param_gen = _mn_param_gen[i];

      // Fit fractions
      const Eigen::MatrixXd ffgen = FitFrac(mn_param_gen);
      const double ffsum = ffgen.trace();

      const auto ff_sumcols_gen = ffgen.colwise().sum().array();
      const double ff_vv_SD_gen = ff_sumcols_gen(0) + ff_sumcols_gen(2);
      const double ff_vv_P_gen = ff_sumcols_gen(1);
      const double ff_vs_plus_gen = ff_sumcols_gen(3);
      const double ff_vs_minus_gen = ff_sumcols_gen(4);
      const double ff_ss_gen = ff_sumcols_gen(5);
      const double ff_scalar_gen = ff_vs_plus_gen + ff_vs_minus_gen + ff_ss_gen;

      const Eigen::MatrixXd ff_diff2 = (ffgen - _ff).array().pow(2.0);
      _ff_err += ff_diff2;
      _ffsum_err += std::pow(ffsum - _ffsum, 2.0);

      /* _ff_vector_err += std::pow(ff_vector_gen - _ff_vector, 2); */
      /* _ff_scalar_err += std::pow(ff_scalar_gen - _ff_scalar, 2); */

      _ff_vv_SD_err += std::pow(ff_vv_SD_gen - _ff_vv_SD, 2);
      _ff_vv_P_err += std::pow(ff_vv_P_gen - _ff_vv_P, 2);
      _ff_vs_plus_err += std::pow(ff_vs_plus_gen - _ff_vs_plus, 2);
      _ff_vs_minus_err += std::pow(ff_vs_minus_gen - _ff_vs_minus, 2);
      _ff_ss_err += std::pow(ff_ss_gen - _ff_ss, 2);
      _ff_scalar_err += std::pow(ff_scalar_gen - _ff_scalar, 2);

      // Phases
      const Eigen::VectorXd phasesgen = Phases(mn_param_gen);

      // const double arg_vv_S_gen = phasesgen(0);
      // const double arg_vv_P_gen = phasesgen(1);
      const double arg_vv_D_gen = phasesgen(2);
      const double arg_vs_plus_gen = phasesgen(3);
      const double arg_vs_minus_gen = phasesgen(4);
      const double arg_ss_gen = phasesgen(5);

      // _arg_vv_S_err += std::pow(arg_vv_S_gen - _arg_vv_S, 2);
      // _arg_vv_P_err += std::pow(arg_vv_P_gen - _arg_vv_P, 2);
      _arg_vv_D_err += std::pow(arg_vv_D_gen - _arg_vv_D, 2);
      _arg_vs_plus_err += std::pow(arg_vs_plus_gen - _arg_vs_plus, 2);
      _arg_vs_minus_err += std::pow(arg_vs_minus_gen - _arg_vs_minus, 2);
      _arg_ss_err += std::pow(arg_ss_gen - _arg_ss, 2);
    }

    // Fit fractions
    _ff_err = _ff_err.array().sqrt();
    _ff_err /= sqrt(static_cast<double>(nexp));
    _ffsum_err = sqrt(_ffsum_err / static_cast<double>(nexp));

    /* _ff_vector_err = sqrt(_ff_vector_err / static_cast<double>(nexp)); */
    /* _ff_scalar_err = sqrt(_ff_scalar_err / static_cast<double>(nexp)); */

    _ff_vv_SD_err = sqrt(_ff_vv_SD_err / static_cast<double>(nexp));
    _ff_vv_P_err = sqrt(_ff_vv_P_err / static_cast<double>(nexp));
    _ff_vs_plus_err = sqrt(_ff_vs_plus_err / static_cast<double>(nexp));
    _ff_vs_minus_err = sqrt(_ff_vs_minus_err / static_cast<double>(nexp));
    _ff_ss_err = sqrt(_ff_ss_err / static_cast<double>(nexp));
    _ff_scalar_err = sqrt(_ff_scalar_err / static_cast<double>(nexp));

    // Phases
    // _arg_vv_S_err = sqrt(_arg_vv_S_err / static_cast<double>(nexp));
    // _arg_vv_P_err = sqrt(_arg_vv_P_err / static_cast<double>(nexp));
    _arg_vv_D_err = sqrt(_arg_vv_D_err / static_cast<double>(nexp));
    _arg_vs_plus_err = sqrt(_arg_vs_plus_err / static_cast<double>(nexp));
    _arg_vs_minus_err = sqrt(_arg_vs_minus_err / static_cast<double>(nexp));
    _arg_ss_err = sqrt(_arg_ss_err / static_cast<double>(nexp));
  }
};
