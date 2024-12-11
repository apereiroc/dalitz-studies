#pragma once

#include "Minuit2/MnUserParameters.h"
#include <map>
#include <vector>

/**
   Get vector of variable parameters from full set
 */
inline std::vector<double>
ext_to_int(const ROOT::Minuit2::MnUserParameters &mn_param) {
  std::vector<double> par_int;

  const std::vector<double> &par = mn_param.Params();
  for (unsigned int i = 0; i < mn_param.Trafo().VariableParameters(); ++i) {
    par_int.push_back(par.at(mn_param.Trafo().ExtOfInt(i)));
  }

  return par_int;
}

/**
   Map internal parameter indices onto external parameter indices
*/
inline std::map<unsigned int, unsigned int>
int_to_ext(const ROOT::Minuit2::MnUserParameters &mn_param) {
  std::map<unsigned int, unsigned int> int_to_ext;
  for (unsigned int i = 0; i < mn_param.Trafo().VariableParameters(); ++i)
    int_to_ext.insert(
        std::pair<unsigned int, unsigned int>(i, mn_param.Trafo().ExtOfInt(i)));

  return int_to_ext;
}
