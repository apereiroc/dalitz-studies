#pragma once

#include <Eigen/Dense>

#include "TRandom3.h"

#include "FitCovariance.h"
#include "FitParameters.h"
#include "Minuit2Util.h"

/**
   Collection of useful utilities
 */

/**
   Convert Eigen parameter vector to std vector
 */
inline std::vector<double> Eigen_to_std(const Eigen::VectorXd &par_eigen) {
  const unsigned int nrow = par_eigen.size();

  std::vector<double> par_std(nrow);
  for (unsigned int i = 0; i < nrow; i++)
    par_std.at(i) = par_eigen(i);

  return par_std;
}

/**
   Convert Minuit2 parameter vector to Eigen vector
 */
inline Eigen::VectorXd
Minuit2_to_Eigen(const std::vector<double> &par_minuit2) {
  const unsigned int nrow = par_minuit2.size();

  Eigen::VectorXd par_eigen(nrow);
  for (unsigned int i = 0; i < nrow; i++)
    par_eigen(i) = par_minuit2.at(i);

  return par_eigen;
}

/**
   Convert Minuit2 covariance matrix to Eigen matrix
 */
inline Eigen::MatrixXd
Minuit2_to_Eigen(const ROOT::Minuit2::MnUserCovariance &cov_minuit2) {
  const unsigned int nrow = cov_minuit2.Nrow();

  Eigen::MatrixXd cov_eigen(nrow, nrow);
  for (unsigned int i = 0; i < nrow; i++)
    for (unsigned int j = 0; j < nrow; j++) {
      cov_eigen(i, j) = cov_minuit2(i, j);
    }

  return cov_eigen;
}

/**
   Reduce parameter vector with initialiser list of external indices
 */
inline Eigen::VectorXd
reduce_par(const ROOT::Minuit2::MnUserParameters &par_minuit2,
           const std::initializer_list<unsigned int> &indices) {
  const std::vector<unsigned int> &v_indices(indices);

  Eigen::VectorXd par_eigen(indices.size());
  for (unsigned int i = 0; i < indices.size(); ++i)
    par_eigen(i) = par_minuit2.Value(v_indices.at(i));

  return par_eigen;
}

/**
   Reduce covariance matrix with initialiser list of internal indices
 */
inline Eigen::MatrixXd
reduce_cov(const ROOT::Minuit2::MnUserCovariance &cov_minuit2,
           const std::initializer_list<unsigned int> &indices) {
  const std::vector<unsigned int> &v_indices(indices);

  Eigen::MatrixXd cov_eigen(indices.size(), indices.size());
  for (unsigned int i = 0; i < indices.size(); ++i)
    for (unsigned int j = 0; j < indices.size(); ++j)
      cov_eigen(i, j) = cov_minuit2(v_indices.at(i), v_indices.at(j));

  return cov_eigen;
}

/**
   Generate random set of parameters from a given set of parameters and its
   associated covariance matrix
   D.R. Barr and N.L. Slezak, Commun. ACM 15(12), p1048 (1972)
 */
inline Eigen::VectorXd GenGaussianXd(const Eigen::VectorXd &par,
                                     const Eigen::MatrixXd &cov) {
  const unsigned int ndim = par.size();

  // Cholesky decomposition
  const Eigen::MatrixXd L = cov.llt().matrixL();

  // Independent standard normal variates
  Eigen::VectorXd z(ndim);
  for (unsigned int i = 0; i < ndim; ++i)
    z(i) = gRandom->Gaus();

  const Eigen::VectorXd par_gen = par + (L * z);

  return par_gen;
}

/**
   Generate set of parameters from fit result
 */

inline FitParameters GenParam(const FitParameters &mn_param,
                              const FitCovariance &cov) {
  // Map internal parameter indices onto external parameter indices
  const std::map<unsigned int, unsigned int> int_to_ext_map =
      int_to_ext(mn_param);

  const Eigen::VectorXd par_gen_int = GenGaussianXd(
      Minuit2_to_Eigen(ext_to_int(mn_param)), Minuit2_to_Eigen(cov));

  FitParameters mn_param_gen = mn_param;
  for (unsigned int j = 0; j < par_gen_int.size(); ++j)
    mn_param_gen.SetValue(int_to_ext_map.at(j), par_gen_int(j));

  return mn_param_gen;
}

inline std::vector<FitParameters> GenParam(const FitParameters &mn_param,
                                           const FitCovariance &cov,
                                           const unsigned int &nexp) {
  std::vector<FitParameters> mn_params;
  mn_params.reserve(nexp);

  for (unsigned int i = 0; i < nexp; ++i) {
    FitParameters mn_param_gen = GenParam(mn_param, cov);

    mn_params.push_back(mn_param_gen);
  }

  return mn_params;
}
