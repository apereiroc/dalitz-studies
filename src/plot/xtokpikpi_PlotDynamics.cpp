#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "boost/program_options/parsers.hpp"
#include <boost/program_options.hpp>

#include "AmpSS.h"
#include "AmpVS.h"
#include "AmpVV.h"
#include "Data.h"
#include "Event.h"
#include "FitParameters.h"
#include "Meson.h"
#include "SigPDF.h"

using namespace ROOT::Minuit2;
namespace po = boost::program_options;
using json = nlohmann::ordered_json;

int main(const int argc, const char *argv[]) {
  std::string inMeson, inNorm, inPars, outFile;

  po::options_description desc{"Options"};

  auto op = desc.add_options();
  op("help,h", "Display usage");
  op("meson,m", po::value<std::string>(&inMeson), "Input meson (Bs/Du)");
  op("normfile,F", po::value<std::string>(&inNorm), "Input normalisation file");
  op("input-pars,i", po::value<std::string>(&inPars), "Input parameter file");
  op("output-file,o", po::value<std::string>(&outFile), "Output file");

  po::variables_map args;
  po::store(po::parse_command_line(argc, argv, desc), args);
  po::notify(args);

  if (args.count("help")) {
    std::cout << desc << std::endl;
    std::exit(0);
  }

  Meson mother(inMeson);

  std::vector<Event> norm;

  std::cout << "Loading integration events..." << std::endl;
  LoadToy(inNorm, "fitTree", norm, true);

  // Load fit parameters
  FitParameters mn_param;
  mn_param.LoadParFromJSON(inPars);
  const std::vector<double> &par = mn_param.Params();

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
  pdf_sig.ResizeEvents(norm);
  pdf_sig.Norm(norm, par);

  std::cout << "Inserting amps in the events...\n";
#pragma omp parallel for
  for (auto &event : norm) {
    pdf_sig.updateAmpsInEvent(event, par);
  }

  std::cout << "Creating plot...\n";

  TFile *const file = new TFile(outFile.c_str(), "recreate");

  auto *const treeVariables = new TTree("variables", "");

  double cos1, cos2, phi, mass1, mass2;

  treeVariables->Branch("cos1", &cos1);
  treeVariables->Branch("cos2", &cos2);
  treeVariables->Branch("phi", &phi);
  treeVariables->Branch("mass1", &mass1);
  treeVariables->Branch("mass2", &mass2);

  auto *const treeVV = new TTree("VV", "");

  double spin_S, spin_P, spin_D;
  double bf_S, bf_P, bf_D;

  treeVV->Branch("spin_S", &spin_S);
  treeVV->Branch("spin_P", &spin_P);
  treeVV->Branch("spin_D", &spin_D);

  treeVV->Branch("bf_S", &bf_S);
  treeVV->Branch("bf_P", &bf_P);
  treeVV->Branch("bf_D", &bf_D);

  for (const auto &event : norm) {
    // Fill phase spoace points
    cos1 = event.GetCosTheta1();
    cos2 = event.GetCosTheta2();
    phi = event.GetChi();
    mass1 = event.GetMassKpPim();
    mass2 = event.GetMassKmPip();

    treeVariables->Fill();

    // Fill dynamical information
    // const double detJ = get_detJ(mass1, mass2, mass_Bs);

    spin_S = std::abs(event.GetVV().GetAmpSpin(SpinVV::S, CPConf::A));
    spin_P = std::abs(event.GetVV().GetAmpSpin(SpinVV::P, CPConf::A));
    spin_D = std::abs(event.GetVV().GetAmpSpin(SpinVV::D, CPConf::A));

    bf_S = std::abs(event.GetVV().GetAmpBarrier(SpinVV::S, CPConf::A));
    bf_P = std::abs(event.GetVV().GetAmpBarrier(SpinVV::P, CPConf::A));
    bf_D = std::abs(event.GetVV().GetAmpBarrier(SpinVV::D, CPConf::A));

    treeVV->Fill();
  }

  treeVariables->Write();
  treeVV->Write();

  file->Close();

  return 0;
}
