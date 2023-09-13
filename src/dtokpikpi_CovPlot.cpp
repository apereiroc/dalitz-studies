#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "boost/program_options/parsers.hpp"
#include <boost/program_options.hpp>

#include "AmpVV.h"
#include "AmpVS.h"
#include "AmpSS.h"
#include "Data.h"
#include "Event.h"
#include "FitParameters.h"
#include "SigPDF.h"

using namespace ROOT::Minuit2;
namespace po = boost::program_options;
using json = nlohmann::ordered_json;

int main(const int argc, const char *argv[]) {
  std::string inSig, inNorm, inPars, outFile;

  po::options_description desc{"Options"};

  desc.add_options()("help,h", "Display usage")(
          "sigfile,f", po::value<std::string>(&inSig),
          "Input signal file")("normfile,F", po::value<std::string>(&inNorm),
                               "Input normalisation file")(
          "input-pars,i", po::value<std::string>(&inPars),
          "Input parameter file")(
          "output-file,o", po::value<std::string>(&outFile), "Output file");

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

  const std::vector<double> &par = mn_param.Params();

  // Load amplitudes
  AmpVV ampVV_S("K*(892)0 K*(892)0b [S]", SpinVV::S,
                abs_VV_S, arg_VV_S,
                absLambda_VV_S, argLambda_VV_S);

  AmpVV ampVV_P("K*(892)0 K*(892)0b [P]", SpinVV::P,
                abs_VV_P, arg_VV_P,
                absLambda_VV_P, argLambda_VV_P);

  AmpVV ampVV_D("K*(892)0 K*(892)0b [D]", SpinVV::D,
                abs_VV_D, arg_VV_D,
                absLambda_VV_D, argLambda_VV_D);

  AmpVS ampVS_p("K*(892) (Kpi)0 CP-odd", VSConf::Plus,
                abs_VS_plus, arg_VS_plus,
                absLambda_VS_plus, argLambda_VS_plus);

  AmpVS ampVS_m("K*(892) (Kpi)0 CP-even", VSConf::Minus,
                abs_VS_minus, arg_VS_minus,
                absLambda_VS_minus, argLambda_VS_minus);

  AmpSS ampSS("(Kpi)0 (Kpi)0bar",
              abs_SS, arg_SS,
              absLambda_SS, argLambda_SS);

  ampVV_S.SetPropagator(PropConf::BW);
  ampVV_P.SetPropagator(PropConf::BW);
  ampVV_D.SetPropagator(PropConf::BW);
  ampVS_p.SetPropagatorV(PropConf::BW);
  ampVS_p.SetPropagatorS(PropConf::LASS);
  ampVS_m.SetPropagatorV(PropConf::BW);
  ampVS_m.SetPropagatorS(PropConf::LASS);
  ampSS.SetPropagator(PropConf::LASS);

  const auto amps =
          std::make_tuple(ampVV_S, ampVV_P, ampVV_D,
                          ampVS_p, ampVS_m, ampSS);

  // Load PDFs
  SigPDF pdf_sig(amps, tau_Du, DG_Du, Dm_Du);
  pdf_sig.NormTime(par);
  pdf_sig.ResizeEvents(data, norm);
  pdf_sig.Norm(norm, par);
  std::cout << "Creating plot..." << std::endl;

  auto *const file = new TFile(outFile.c_str(), "recreate");
  auto *const treeData = new TTree("data", "");

  double cos1, cos2, phi, m1, m2, w, pdf, time;

  treeData->Branch("cos1", &cos1);
  treeData->Branch("cos2", &cos2);
  treeData->Branch("phi", &phi);
  treeData->Branch("m1", &m1);
  treeData->Branch("m2", &m2);
  treeData->Branch("time", &time);
  treeData->Branch("weight", &w);

  for (const auto &event: data) {
    w = 1.0;

    cos1 = event.GetCosTheta1();
    cos2 = event.GetCosTheta2();
    phi = event.GetChi();
    m1 = event.GetMassKpPim();
    m2 = event.GetMassKmPip();
    time = event.GetTime();

    treeData->Fill();
  }

  treeData->Write();

  auto *const treeFit = new TTree("fit", "");

  treeFit->Branch("cos1", &cos1);
  treeFit->Branch("cos2", &cos2);
  treeFit->Branch("phi", &phi);
  treeFit->Branch("m1", &m1);
  treeFit->Branch("m2", &m2);
  treeFit->Branch("time", &time);
  treeFit->Branch("pdf", &pdf);

  const double norm_pdf = pdf_sig.GetTimeIntegratedNorm(par);

  for (const auto &event: norm) {
    pdf = pdf_sig.GetTimeIntegratedPDF(event, par) / norm_pdf;

    cos1 = event.GetCosTheta1();
    cos2 = event.GetCosTheta2();
    phi = event.GetChi();
    m1 = event.GetMassKpPim();
    m2 = event.GetMassKmPip();
    time = event.GetTime();

    treeFit->Fill();
  }

  treeFit->Write();
  file->Close();

  return 0;
}
