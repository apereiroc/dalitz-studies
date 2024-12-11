#include "AmpSS.h"
#include "AmpVS.h"
#include "AmpVV.h"
#include "Data.h"
#include "FitCovariance.h"
#include "FitParameters.h"
#include "Meson.h"
#include "Par.h"
#include "Results.h"
#include "SigPDF.h"
#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

namespace po = boost::program_options;

// inline void AddMassParametersToJson(const FitParameters &mn_param, json &out,
//                                     bool fit_failed = false) {
//   for (unsigned int i = Par::mass_kst0892; i <= Par::c4_arg_pelaez; i++) {
//     if (mn_param.Parameter(i).IsFixed())
//       continue;
//
//     // Create new parameter
//     json this_param;
//
//     // Fill info
//     this_param["val"] = fit_failed ? 999999.0 : mn_param.Value(i);
//     this_param["err"] = fit_failed ? 999999.0 : mn_param.Error(i);
//     this_param["fixed"] = false;
//     if (mn_param.Parameter(i).HasLowerLimit())
//       this_param["ll"] = mn_param.Parameter(i).LowerLimit();
//     else
//       this_param["ll"] = nullptr;
//
//     if (mn_param.Parameter(i).HasUpperLimit())
//       this_param["ul"] = mn_param.Parameter(i).UpperLimit();
//     else
//       this_param["ul"] = nullptr;
//
//     this_param["latex"] = mn_param.GetLatex(i);
//
//     // Add to output json file
//     out[mn_param.GetName(i)] = this_param;
//   }
// }

int main(const int argc, const char *argv[]) {

  std::string inMeson, normFile;
  std::string inPars, inCov;
  std::string outResults;
  unsigned int nexperiments;

  po::options_description desc{"Options"};

  auto op = desc.add_options();
  op("help,h", "Display usage");
  op("meson,m", po::value<std::string>(&inMeson), "Input meson (Bs/Bd)");
  op("normfile,F", po::value<std::string>(&normFile),
     "Input normalisation ROOT file");
  op("input-pars,i", po::value<std::string>(&inPars), "Input parameter file");
  op("input-cov,c", po::value<std::string>(&inCov),
     "Input covariance matrix file");
  op("nexperiments,n",
     po::value<unsigned int>(&nexperiments)->default_value(1000),
     "Number of experiments");
  op("output-results,o", po::value<std::string>(&outResults),
     "Output results file");

  po::variables_map args;
  po::store(po::parse_command_line(argc, argv, desc), args);
  po::notify(args);
  if (args.count("help")) {
    std::cout << desc << std::endl;
    std::exit(0);
  }

  FitParameters mn_param;
  std::cout << "Loading parameters..." << std::endl;
  mn_param.LoadParFromJSON(inPars);
  const std::vector<double> &par = mn_param.Params();

  Meson mother(inMeson);

  std::vector<Event> norm;

  std::cout << "Loading integration events..." << std::endl;
  LoadToy(normFile, "fitTree", norm, true);

  std::cout << "Parameter state: " << mn_param << std::endl;

  // Load amps

  // AmpVVSum ampVV_SD("K*(892)0 K*(892)0b [S+D]", SpinVV::S, SpinVV::D,
  //                   x_kst0892kst0892bar_S, y_kst0892kst0892bar_S,
  //                   x_kst0892kst0892bar_D, y_kst0892kst0892bar_D,
  //                   absL_kst0892kst0892bar_S, argL_kst0892kst0892bar_S,
  //                   absL_kst0892kst0892bar_D, argL_kst0892kst0892bar_D);

  AmpVV ampVV_S("K*(892)0 K*(892)0b [S]", SpinVV::S, Par::x_VV_S, Par::y_VV_S,
                Par::absLambda_VV_S, Par::argLambda_VV_S);

  AmpVV ampVV_P("K*(892)0 K*(892)0b [P]", SpinVV::P, Par::x_VV_P, Par::y_VV_P,
                Par::absLambda_VV_S, Par::argLambda_VV_S);

  AmpVV ampVV_D("K*(892)0 K*(892)0b [D]", SpinVV::D, Par::x_VV_D, Par::y_VV_D,
                Par::absLambda_VV_S, Par::argLambda_VV_S);

  AmpVS ampVS_p("K*(892) (Kpi)0 CP-odd", VSConf::Plus, Par::x_VS_plus,
                Par::y_VS_plus, Par::absLambda_VS_plus, Par::argLambda_VS_plus);

  AmpVS ampVS_m("K*(892) (Kpi)0 CP-even", VSConf::Minus, Par::x_VS_minus,
                Par::y_VS_minus, Par::absLambda_VS_minus,
                Par::argLambda_VS_minus);

  AmpSS ampSS("(Kpi)0 (Kpi)0bar", Par::x_SS, Par::y_SS, Par::absLambda_SS,
              Par::argLambda_SS);

  ampVV_S.SetPropagator(PropConf::BW);
  ampVV_P.SetPropagator(PropConf::BW);
  ampVV_D.SetPropagator(PropConf::BW);
  ampVS_p.SetPropagatorV(PropConf::BW);
  ampVS_p.SetPropagatorS(PropConf::LASS);
  ampVS_m.SetPropagatorV(PropConf::BW);
  ampVS_m.SetPropagatorS(PropConf::LASS);
  ampSS.SetPropagator(PropConf::LASS);

  const auto amps =
      std::make_tuple(ampVV_S, ampVV_P, ampVV_D, ampVS_p, ampVS_m, ampSS);

  SigPDF pdf_sig(amps, mother.getTauIdx(), mother.getDGIdx(),
                 mother.getDmIdx());
  pdf_sig.NormTime(par);
  pdf_sig.ResizeEvents(norm);
  pdf_sig.Norm(norm, par);

  Results results(pdf_sig, norm, mn_param);

  std::cout << "Loading covariance..." << std::endl;
  FitCovariance cov(inCov);
  results.LoadCov(cov);

  results.FitFrac();
  results.Phases();
  results.FitFracErr(nexperiments);

  // Get vector fractions
  const double f_vector_SD_val = results.GetFF_VV_SD();
  const double f_vector_SD_err = results.GetFF_VV_SD_Err();
  const double f_vector_P_val = results.GetFF_VV_P();
  const double f_vector_P_err = results.GetFF_VV_P_Err();

  // Get scalar fractions
  const double f_scalar_VS_plus_val = results.GetFF_VS_Plus();
  const double f_scalar_VS_plus_err = results.GetFF_VS_Plus_Err();
  const double f_scalar_VS_minus_val = results.GetFF_VS_Minus();
  const double f_scalar_VS_minus_err = results.GetFF_VS_Minus_Err();
  const double f_scalar_SS_val = results.GetFF_SS();
  const double f_scalar_SS_err = results.GetFF_SS_Err();
  const double f_scalar_val = results.GetFF_Scalar();
  const double f_scalar_err = results.GetFF_Scalar_Err();

  // Get phases
  const double delta_VV_D_val = results.GetArg_VV_D();
  const double delta_VV_D_err = results.GetArg_VV_D_Err();
  const double delta_VS_plus_val = results.GetArg_VS_Plus();
  const double delta_VS_plus_err = results.GetArg_VS_Plus_Err();
  const double delta_VS_minus_val = results.GetArg_VS_Minus();
  const double delta_VS_minus_err = results.GetArg_VS_Minus_Err();
  const double delta_SS_val = results.GetArg_SS();
  const double delta_SS_err = results.GetArg_SS_Err();

  // Print fit fractions
  const auto ff = results.FitFrac(mn_param);
  std::cout << "Fit fractions:\n" << ff << std::endl;

  // Save parameters
  json out;

  // add fit fractions and phases
  json f_VV_SD, f_VV_P, f_VS_plus, f_VS_minus, f_SS, f_scalar;
  json delta_VV_D, delta_VS_plus, delta_VS_minus, delta_SS;

  /* f_vector["val"] = f_vector_val; */
  /* f_vector["err"] = f_vector_err; */
  /* f_vector["fixed"] = true; */
  /* f_vector["ll"] = nullptr; */
  /* f_vector["ul"] = nullptr; */
  /* f_vector["latex"] = "$f_{\rm VV}$"; */

  f_scalar["val"] = f_scalar_val;
  f_scalar["err"] = f_scalar_err;
  f_scalar["fixed"] = false;
  f_scalar["ll"] = nullptr;
  f_scalar["ul"] = nullptr;
  f_scalar["latex"] = "$f_{\\text{scalar}}^{}$";

  f_VV_SD["val"] = f_vector_SD_val;
  f_VV_SD["err"] = f_vector_SD_err;
  f_VV_SD["fixed"] = false;
  f_VV_SD["ll"] = nullptr;
  f_VV_SD["ul"] = nullptr;
  f_VV_SD["latex"] = "$f_{VV}^{S+D}$";

  f_VV_P["val"] = f_vector_P_val;
  f_VV_P["err"] = f_vector_P_err;
  f_VV_P["fixed"] = false;
  f_VV_P["ll"] = nullptr;
  f_VV_P["ul"] = nullptr;
  f_VV_P["latex"] = "$f_{VV}^{P}$";

  f_VS_plus["val"] = f_scalar_VS_plus_val;
  f_VS_plus["err"] = f_scalar_VS_plus_err;
  f_VS_plus["fixed"] = false;
  f_VS_plus["ll"] = nullptr;
  f_VS_plus["ul"] = nullptr;
  f_VS_plus["latex"] = "$f_{VS}^{+}$";

  f_VS_minus["val"] = f_scalar_VS_minus_val;
  f_VS_minus["err"] = f_scalar_VS_minus_err;
  f_VS_minus["fixed"] = false;
  f_VS_minus["ll"] = nullptr;
  f_VS_minus["ul"] = nullptr;
  f_VS_minus["latex"] = "$f_{VS}^{-}$";

  f_SS["val"] = f_scalar_SS_val;
  f_SS["err"] = f_scalar_SS_err;
  f_SS["fixed"] = false;
  f_SS["ll"] = nullptr;
  f_SS["ul"] = nullptr;
  f_SS["latex"] = "$f_{SS}^{}$";

  delta_VV_D["val"] = delta_VV_D_val;
  delta_VV_D["err"] = delta_VV_D_err;
  delta_VV_D["fixed"] = false;
  delta_VV_D["ll"] = nullptr;
  delta_VV_D["ul"] = nullptr;
  delta_VV_D["latex"] = "$\\delta_{VV}^{D}$";

  delta_VS_plus["val"] = delta_VS_plus_val;
  delta_VS_plus["err"] = delta_VS_plus_err;
  delta_VS_plus["fixed"] = false;
  delta_VS_plus["ll"] = nullptr;
  delta_VS_plus["ul"] = nullptr;
  delta_VS_plus["latex"] = "$\\delta_{VS}^{+}$";

  delta_VS_minus["val"] = delta_VS_minus_val;
  delta_VS_minus["err"] = delta_VS_minus_err;
  delta_VS_minus["fixed"] = false;
  delta_VS_minus["ll"] = nullptr;
  delta_VS_minus["ul"] = nullptr;
  delta_VS_minus["latex"] = "$\\delta_{VS}^{-}$";

  delta_SS["val"] = delta_SS_val;
  delta_SS["err"] = delta_SS_err;
  delta_SS["fixed"] = false;
  delta_SS["ll"] = nullptr;
  delta_SS["ul"] = nullptr;
  delta_SS["latex"] = "$\\delta_{SS}^{}$";

  out["f_vv_SD"] = f_VV_SD;
  out["f_vv_P"] = f_VV_P;
  out["f_vs_plus"] = f_VS_plus;
  out["f_vs_minus"] = f_VS_minus;
  out["f_ss"] = f_SS;
  out["f_scalar"] = f_scalar;
  out["delta_vv_d"] = delta_VV_D;
  out["delta_vs_plus"] = delta_VS_plus;
  out["delta_vs_minus"] = delta_VS_minus;
  out["delta_ss"] = delta_SS;

  // Save values of interest to json file
  std::ofstream file_results(outResults);
  file_results << std::setw(4) << out << std::endl;

  return 0;
}
