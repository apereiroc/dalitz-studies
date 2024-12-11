#pragma once

#include <fstream>
#include <iostream>

#include "Minuit2/MnUserCovariance.h"

/**
   Enhance MnUserCovariance functionality
 */
class FitCovariance : public ROOT::Minuit2::MnUserCovariance {
public:
  FitCovariance() = default;

  FitCovariance(const ROOT::Minuit2::MnUserCovariance &cov)
      : ROOT::Minuit2::MnUserCovariance(cov) {}

  FitCovariance(const std::string &filename)
      : ROOT::Minuit2::MnUserCovariance(LoadCovData(filename),
                                        LoadCovNrow(filename)) {}

  /**
     Save covariance matrix to text file
   */
  void SaveCov(const std::string &filename) const {
    std::ofstream ofs(filename.c_str(), std::ofstream::out);

    ofs << this->Nrow() << std::endl;
    for (unsigned int i = 0; i < this->Data().size(); ++i)
      ofs << this->Data()[i] << std::endl;

    ofs.close();
  }

private:
  /**
     Load covariance matrix from text file
   */
  static unsigned int LoadCovNrow(const std::string &filename) {
    std::ifstream ifs(filename.c_str());

    if (ifs.is_open() == 0) {
      std::cout << "ERROR: " << filename << " could not be found" << std::endl;
      exit(1);
    }

    unsigned int nrow;
    ifs >> nrow;
    ifs.close();

    return nrow;
  }

  static std::vector<double> LoadCovData(const std::string &filename) {
    std::ifstream ifs(filename.c_str());

    if (ifs.is_open() == 0) {
      std::cout << "ERROR: " << filename << " could not be found" << std::endl;
      exit(1);
    }

    std::vector<double> data;

    unsigned int dummy;
    double cov;
    ifs >> dummy;
    ifs >> cov;
    while (ifs.eof() == 0) {
      data.push_back(cov);
      ifs >> cov;
    }
    ifs.close();

    return data;
  }
};
