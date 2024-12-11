#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "FitCovariance.h"
#include "FitUtil.h"
#include "boost/program_options/parsers.hpp"
#include <boost/program_options.hpp>

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnPrint.h"
#include "TRandom3.h"

#include "AmpSS.h"
#include "AmpVS.h"
#include "AmpVV.h"
#include "Data.h"
#include "Event.h"
#include "FitParameters.h"
#include "Meson.h"
#include "MinuitFcn.h"
#include "SigPDF.h"

using namespace ROOT::Minuit2;
namespace po = boost::program_options;

int main(const int argc, const char *argv[]) {
  std::string inMeson, inSig, inNorm, inPars, outPars, outCov;

  po::options_description desc{"Options"};

  auto op = desc.add_options();
  op("help,h", "Display usage");
  op("meson,m", po::value<std::string>(&inMeson), "Input meson (Bs/Du)");
  op("sigfile,f", po::value<std::string>(&inSig), "Input signal file");
  op("normfile,F", po::value<std::string>(&inNorm), "Input normalisation file");
  op("input-pars,i", po::value<std::string>(&inPars), "Input parameter file");
  op("output-pars,o", po::value<std::string>(&outPars),
     "Output parameter file");
  op("output-cov,c", po::value<std::string>(&outCov),
     "Output covariance matrix file");

  po::variables_map args;
  po::store(po::parse_command_line(argc, argv, desc), args);
  po::notify(args);

  if (args.count("help")) {
    std::cout << desc << std::endl;
    std::exit(0);
  }

  Meson mother(inMeson);

  std::vector<Event> data, norm;

  std::cout << "Loading data..." << std::endl;
  LoadToy(inSig, "fitTree", data, true);

  std::cout << "Loading integration events..." << std::endl;
  LoadToy(inNorm, "fitTree", norm, true);

  // Load fit parameters
  FitParameters mn_param;
  mn_param.LoadParFromJSON(inPars);

  // Avoid vanishing gradients
  for (unsigned int idx = x_VV_S; idx <= y_SS; idx++) {
    if (!mn_param.Parameter(idx).IsFixed()) {
      mn_param.SetValue(idx, gRandom->Gaus(mn_param.Value(idx), 0.1));
    }
  }

  // Fix D wave
  // mn_param.SetAndFix(Par::x_VV_D, 0.0);
  // mn_param.SetAndFix(Par::y_VV_D, 0.0);

  const std::vector<double> &par = mn_param.Params();

  std::cout << "Initial parameter state: " << mn_param << std::endl;

  // Load amplitudes
  AmpVV ampVV_S("K*(892)0 K*(892)0b [S]", SpinVV::S, x_VV_S, y_VV_S,
                absLambda_VV_S, argLambda_VV_S);

  AmpVV ampVV_P("K*(892)0 K*(892)0b [P]", SpinVV::P, x_VV_P, y_VV_P,
                absLambda_VV_P, argLambda_VV_P);

  AmpVV ampVV_D("K*(892)0 K*(892)0b [D]", SpinVV::D, x_VV_D, y_VV_D,
                absLambda_VV_D, argLambda_VV_D);

  AmpVS ampVS_p("K*(892) (Kpi)0 CP-odd", VSConf::Plus, x_VS_plus, y_VS_plus,
                absLambda_VS_plus, argLambda_VS_plus);

  AmpVS ampVS_m("K*(892) (Kpi)0 CP-even", VSConf::Minus, x_VS_minus, y_VS_minus,
                absLambda_VS_minus, argLambda_VS_minus);

  AmpSS ampSS("(Kpi)0 (Kpi)0bar", x_SS, y_SS, absLambda_SS, argLambda_SS);

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

  // Load PDFs
  SigPDF pdf_sig(amps, mother.getTauIdx(), mother.getDGIdx(),
                 mother.getDmIdx());
  pdf_sig.NormTime(par);
  pdf_sig.ResizeEvents(data, norm);

  // Load Fcn
  MinuitFcn fcn(data, norm, pdf_sig);

  // Double minimisation with strategy 2
  MnMigrad migrad(fcn, mn_param, 2);
  FunctionMinimum min = migrad();
  min = migrad();
  print_minimum_info(min);

  if (min.IsValid())
    mn_param.SaveParToJSON(outPars, min);
  if (min.HasValidCovariance()) {
    const FitCovariance &cov = min.UserCovariance();
    cov.SaveCov(outCov);
  }

  std::cout << "Fit params: " << mn_param << std::endl;

  std::cout << "Best -2logL: " << min.Fval() << std::endl;
  std::cout << "       ndof: " << mn_param.VariableParameters() << std::endl;

  return 0;
}
