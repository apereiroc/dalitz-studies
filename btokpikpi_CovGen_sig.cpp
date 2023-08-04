#include <boost/program_options.hpp>
#include <filesystem>
#include <iostream>

#include "AmpVV.h"
#include "Data.h"
#include "FitParameters.h"
#include "GenUtil.h"
#include "Par.h"
#include "Selection.h"
#include "SigPDF.h"

namespace po = boost::program_options;

int main(const int argc, const char *argv[]) {
  std::string inPars, outFile;
  unsigned int seed, nevents;
  double hmax;

  po::options_description desc{"Options"};

  desc.add_options()("help,h", "Display usage")(
      "input-pars,i", po::value<std::string>(&inPars), "Input parameter file")(
      "output-file,o", po::value<std::string>(&outFile), "Output file")(
      "seed,s", po::value<unsigned int>(&seed)->default_value(123),
      "Seed for RNG")("nevents,n",
                      po::value<unsigned int>(&nevents)->default_value(300),
                      "Number of events to be generated")(
      "hmax,h", po::value<double>(&hmax), "Maximum height of the distribution");

  po::variables_map args;
  po::store(po::parse_command_line(argc, argv, desc), args);
  po::notify(args);

  if (args.count("help")) {
    std::cout << desc << std::endl;
    std::exit(0);
  }

  // Load fit parameters
  FitParameters mn_param;
  mn_param.LoadParFromJSON(inPars);
  const std::vector<double> &par = mn_param.Params();

  // Load amplitudes
  AmpVV ampVV_S("K*(892)0 K*(892)0b [S]", SpinVV::S, abs_VV_S, abs_VV_S,
                absLambda_VV_S, argLambda_VV_S);

  AmpVV ampVV_P("K*(892)0 K*(892)0b [P]", SpinVV::P, abs_VV_P, abs_VV_P,
                absLambda_VV_P, argLambda_VV_P);

  AmpVV ampVV_D("K*(892)0 K*(892)0b [D]", SpinVV::D, abs_VV_D, abs_VV_D,
                absLambda_VV_D, argLambda_VV_D);

  ampVV_S.SetPropagator(PropConf::BW);
  ampVV_P.SetPropagator(PropConf::BW);
  ampVV_D.SetPropagator(PropConf::BW);

  const auto amps = std::make_tuple(ampVV_S, ampVV_P, ampVV_D);

  // Load PDFs
  SigPDF pdf_sig(amps, tau_Bs, DG_Bs, Dm_Bs);
  pdf_sig.NormTime(par);

  // Prepare phase space generator
  TLorentzVector P(0.0, 0.0, 0.0, mass_Bs);

  const std::vector<double> masses = {mass_Kp, mass_pip, mass_Kp, mass_pip};

  TGenPhaseSpace generator;
  gRandom->SetSeed(seed);
  generator.SetDecay(P, masses.size(), masses.data());

  std::vector<Event> gen;
  gen.reserve(nevents);

  // Start generation
  for (unsigned int i = 0; i < nevents; ++i) {
    std::cout << "Event: " << i << std::endl;
    while (true) {
      generate_flat_event(generator);

      // Selection criteria
      const double mass1_kpi =
          ((*generator.GetDecay(0)) + (*generator.GetDecay(1))).M();
      const double mass2_kpi =
          ((*generator.GetDecay(2)) + (*generator.GetDecay(3))).M();

      if (Selection_masses(mass1_kpi, mass2_kpi) == false)
        continue;

      const double tauLong = 1.0 / ((1.0 / mn_param.Value(Par::tau_Bs)) -
                                    fabs(mn_param.Value(Par::DG_Bs) / 2.0));
      const double time_gen = -tauLong * std::log(gRandom->Uniform());

      const int qtag_gen = (gRandom->Uniform() < 0.5) ? -1 : +1;

      Event event(*generator.GetDecay(0), *generator.GetDecay(1),
                  *generator.GetDecay(2), *generator.GetDecay(3), time_gen,
                  qtag_gen);

      const double height_gen =
          gRandom->Uniform(hmax) * std::exp(-time_gen / tauLong);

      // Resize amplitude container inside event
      pdf_sig.ResizeEvents(event);

      // Insert amplitudes into event
      pdf_sig.Amps(event, par);

      // Get PDF
      const double pdf_gen = pdf_sig.GetPDF(event, time_gen, qtag_gen, par);

      if (pdf_gen > hmax) {
        std::cout << "ERROR: Generated PDF value " << pdf_gen
                  << ", larger than maximum " << hmax << std::endl;
        std::exit(1);
      }

      if (height_gen < pdf_gen)
        event.SetGenPDF(pdf_gen);
      else
        continue;

      gen.push_back(event);
      break;
    }
  }

  const std::string fname =
      outFile.size() > 0
          ? outFile
          : "dat/bdstokppimkmpip-norm-" + std::to_string(seed) + ".root";

  SaveData(fname, gen);

  return 0;
}
