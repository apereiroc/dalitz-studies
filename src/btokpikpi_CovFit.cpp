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
#include "Minuit2/MnScan.h"
#include "TRandom3.h"

#include "AmpVV.h"
#include "Data.h"
#include "Event.h"
#include "FitParameters.h"
#include "MinuitFcn.h"
#include "Par.h"
#include "SigPDF.h"

using namespace ROOT::Minuit2;
namespace po = boost::program_options;

int main(const int argc, const char *argv[]) {
  std::string inSig, inNorm, inPars, outPars;

  po::options_description desc{"Options"};

  desc.add_options()("help,h", "Display usage")(
      "sigfile,f", po::value<std::string>(&inSig),
      "Input signal file")("normfile,F", po::value<std::string>(&inNorm),
                           "Input normalisation file")(
      "input-pars,i", po::value<std::string>(&inPars),
      "Input parameter file")("output-pars,o", po::value<std::string>(&outPars),
                              "Output parameter file");

  po::variables_map args;
  po::store(po::parse_command_line(argc, argv, desc), args);
  po::notify(args);

  if (args.count("help")) {
    std::cout << desc << std::endl;
    std::exit(0);
  }

  std::vector<Event> data, norm;

  std::cout << "Loading data..." << std::endl;
  LoadToy(inSig, "fitTree", data, true);

  std::cout << "Loading integration events..." << std::endl;
  LoadToy(inNorm, "fitTree", norm, true);

  // Load fit parameters
  FitParameters mn_param;
  mn_param.LoadParFromJSON(inPars);

  mn_param.SetValue(Par::abs_VV_P,
                    gRandom->Gaus(mn_param.Value(Par::abs_VV_P), 0.5));
  mn_param.SetValue(Par::abs_VV_D,
                    gRandom->Gaus(mn_param.Value(Par::abs_VV_D), 0.5));
  mn_param.SetValue(Par::arg_VV_D,
                    gRandom->Gaus(mn_param.Value(Par::arg_VV_D), 0.5));
  const std::vector<double> &par = mn_param.Params();

  std::cout << "Initial parameter state: " << mn_param << std::endl;

  // Load amplitudes
  AmpVV ampVV_S("K*(892)0 K*(892)0b [S]", SpinVV::S, abs_VV_S, arg_VV_S,
                absLambda_VV_S, argLambda_VV_S);

  AmpVV ampVV_P("K*(892)0 K*(892)0b [P]", SpinVV::P, abs_VV_P, arg_VV_P,
                absLambda_VV_P, argLambda_VV_P);

  AmpVV ampVV_D("K*(892)0 K*(892)0b [D]", SpinVV::D, abs_VV_D, arg_VV_D,
                absLambda_VV_D, argLambda_VV_D);

  ampVV_S.SetPropagator(PropConf::BW);
  ampVV_P.SetPropagator(PropConf::BW);
  ampVV_D.SetPropagator(PropConf::BW);

  const auto amps = std::make_tuple(ampVV_S, ampVV_P, ampVV_D);

  // Load PDFs
  SigPDF pdf_sig(amps, tau_Bs, DG_Bs, Dm_Bs);
  pdf_sig.NormTime(par);
  pdf_sig.ResizeEvents(data, norm);

  // Load Fcn
  MinuitFcn fcn(data, norm, pdf_sig);

  MnMigrad migrad(fcn, mn_param, 2);
  FunctionMinimum min = migrad();
  min = migrad();
  mn_param.SaveParToJSON(outPars, min);

  std::cout << "Fit params: " << mn_param << std::endl;

  return 0;
}
