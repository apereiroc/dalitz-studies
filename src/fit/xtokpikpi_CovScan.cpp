#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "boost/program_options/parsers.hpp"
#include <boost/program_options.hpp>

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnPrint.h"

#include "AmpSS.h"
#include "AmpVS.h"
#include "AmpVV.h"
#include "Data.h"
#include "Event.h"
#include "FitParameters.h"
#include "FitUtil.h"
#include "Meson.h"
#include "MinuitFcn.h"
#include "SigPDF.h"

using namespace ROOT::Minuit2;
namespace po = boost::program_options;

int main(const int argc, const char *argv[]) {
  std::string inMeson, inSig, inNorm, inParam, inPars, outPlot;

  po::options_description desc{"Options"};

  desc.add_options()("help,h", "Display usage")(
      "meson,m", po::value<std::string>(&inMeson), "Input meson (Bs/Du)")(
      "sigfile,f", po::value<std::string>(&inSig),
      "Input signal file")("normfile,F", po::value<std::string>(&inNorm),
                           "Input normalisation file")(
      "input-pars,i", po::value<std::string>(&inPars), "Input parameter file")(
      "parameter,p", po::value<std::string>(&inParam),
      "Parameter to be scanned")("output-pars,o",
                                 po::value<std::string>(&outPlot),
                                 "Output parameter file");

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

  // Create plot file
  //
  auto *file = new TFile(outPlot.c_str(), "recreate");
  auto *tree = new TTree("plot", "");

  double x, pdf;

  tree->Branch("x", &x);
  tree->Branch("pdf", &pdf);

  // Configuration
  const unsigned int npoints = 30;
  const unsigned int nfits_per_point = 30;

  unsigned int paramIdx = mn_param.Index(inParam);

  // const double ll = mn_param.Parameter(paramIdx).LowerLimit();
  // const double ul = mn_param.Parameter(paramIdx).UpperLimit();
  const double ll = -50.0;
  const double ul = +20.0;

  std::cout << std::endl
            << "Scanning parameter " << inParam << " from " << ll << " to "
            << ul << std::endl
            << std::endl;

  mn_param.Fix(paramIdx);

  bool minimum_found = false;

  for (unsigned int ipoint = 0; ipoint < npoints; ipoint++) {
    minimum_found = false;
    x = ll + (ul - ll) * ipoint / npoints;
    pdf = std::numeric_limits<double>::max();
    mn_param.SetValue(paramIdx, x);

    std::cout << "This is point " << ipoint + 1 << "/" << npoints << '\n';
    std::cout << "  Parameter " << inParam << " set to " << x << '\n';

    for (unsigned int ifit = 0; ifit < nfits_per_point; ifit++) {
      std::cout << "    This is fit " << ifit + 1 << "/" << nfits_per_point
                << '\n';
      RandomiseCovariantCoeffs(mn_param, time(NULL));
      MnMigrad migrad(fcn, mn_param, 2);
      FunctionMinimum min = migrad();
      min = migrad();

      if (min.IsValid()) {
        std::cout << "      Found a minimum! ✅\n";
        double minVal = min.Fval();
        minimum_found = true;

        if (minVal < pdf)
          pdf = minVal;
      } else {
        std::cout << "      Minimum not valid ❗️\n";
      }
    }

    if (minimum_found)
      tree->Fill();
  }
  tree->Write();
  file->Close();

  return 0;
}
