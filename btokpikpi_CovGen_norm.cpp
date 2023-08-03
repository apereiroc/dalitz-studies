#include <boost/program_options.hpp>
#include <filesystem>
#include <iostream>

// #include "Data.h"
#include "FitParameters.h"
#include "GenUtil.h"
// #include "SigPDF.h"
#include "Utils.h"
#incldue "Par.h"

namespace po = boost::program_options;

int main(const int argc, const char *argv[]) {
  std::string type, outFile, inPars, propagatorV, propagatorS, motherParticle;
  unsigned int seed, nevents;
  double hmax;

  po::options_description desc{"Options"};

  desc.add_options()("help,h", "Display usage")
    ("type,t", po::value<std::string>(&type)->default_value("sig"), "Type: signal (sig) or normalisation (norm)")
    ("input-pars,i", po::value<std::string>(&inPars), "Input parameter file")
    ("output-file,o", po::value<std::string>(&outFile), "Output file")
    ("mother-particle,p", po::value<std::string>(&motherParticle), "Mother particle (Bs/Bd)")
    ("seed,s", po::value<unsigned int>(&seed)->default_value(123), "Seed for RNG")
    ("nevents,n", po::value<unsigned int>(&nevents)->default_value(300), "Number of events to be generated")
    ("hmax", po::value<double>(&hmax), "Maximum height of the PDF")
    ("propagatorV", po::value<std::string>(&propagatorV)->default_value("BW"), "Vector mass propagator")
    ("propagatorS", po::value<std::string>(&propagatorS)->default_value("LASS"), "Scalar mass propagator");

  po::variables_map args;
  po::store(po::parse_command_line(argc, argv, desc), args);
  po::notify(args);

  if (args.count("help")) {
    std::cout << desc << std::endl;
    std::exit(0);
  }

  // Type must be sig or norm
  if (args.count("type") && (type != "sig" and type != "norm")) {
    std::cout << "ERROR: type must be sig or norm. Try again." << std::endl;
    std::exit(1);
  }


  // Load fit parameters
  FitParameters mn_param;
  mn_param.LoadParFromJSON(inPars);
  const std::vector<double> &par = mn_param.Params();

  std::cout << "Initial parameter state: " << mn_param << std::endl;

  // if (mn_param.Parameter(Par::NSig).IsFixed() == false)
  //   nevents = par[Par::NSig];

  // Load amplitudes
  AmpVV ampVV_S("K*(892)0 K*(892)0b [S]", SpinVV::S, 
                abs_VV_S, abs_VV_S,
                absLambda_VV_S, argLambda_VV_S);

  AmpVV ampVV_P("K*(892)0 K*(892)0b [P]", SpinVV::P, 
                abs_VV_P, abs_VV_P,
                absLambda_VV_S, argLambda_VV_S);

  AmpVV ampVV_D("K*(892)0 K*(892)0b [D]", SpinVV::D, 
                abs_VV_D, abs_VV_D,
                absLambda_VV_S, argLambda_VV_S);

  const auto helper = std::make_unique<Helper>();

  const PropConf propV = helper->parsePropagator(propagatorV);

  ampVV_S.SetPropagator(propV);
  ampVV_P.SetPropagator(propV);
  ampVV_D.SetPropagator(propV);

  const auto amps = std::make_tuple(ampVV_S, ampVV_P, ampVV_D);

  // Load PDFs
  SigPDF pdf_sig(amps);

  // For signal
  double amp2_max = hmax;

  if (type == "norm")
    amp2_max = 1.0;

  TLorentzVector P;

  if (motherParticle == "Bs")
    P.SetPxPyPzE(0.0, 0.0, 0.0, mass_Bs);
  else if (motherParticle == "Bd")
    P.SetPxPyPzE(0.0, 0.0, 0.0, mass_Bd);
  else {
    std::cout << "ERROR: mother particle must be Bs or Bd" << std::endl;
    std::exit(1);
  }

  const std::vector<double> masses = {mass_Kp, mass_pip, mass_Kp, mass_pip};

  TGenPhaseSpace generator;
  gRandom->SetSeed(seed);
  generator.SetDecay(P, masses.size(), masses.data());

  std::vector<Event> gen;
  for (unsigned int i = 0; i < nevents; ++i) {
    while (true) {
      generate_flat_event(generator);

      // Selection criteria
      const double mass1_kpi =
          ((*generator.GetDecay(0)) + (*generator.GetDecay(1))).M();
      const double mass2_kpi =
          ((*generator.GetDecay(2)) + (*generator.GetDecay(3))).M();

      if (Selection_masses(mass1_kpi, mass2_kpi) == false)
        continue;

      const double time_gen = gRandom->Uniform(time_min, time_max);

      const int qtag_gen = (gRandom->Uniform() < 0.5) ? -1 : +1;

      Event event(*generator.GetDecay(0), *generator.GetDecay(1),
                  *generator.GetDecay(2), *generator.GetDecay(3), time_gen,
                  qtag_gen);

      double pdf_gen = 1.0;

      event.SetGenPDF(pdf_gen);

      if (type == "sig") {
        const double height_gen = gRandom->Uniform(amp2_max);

        pdf_sig.ResizeEvents(event);
        pdf_sig.Amps(event, par);
        pdf_gen = pdf_sig.GetPDF(event, time_gen, qtag_gen, par);

        if (pdf_gen > amp2_max) {
          std::cout << "ERROR: Generated PDF value " << pdf_gen
                    << ", larger than maximum " << amp2_max << std::endl;
          amp2_max = 2.0 * pdf_gen;
          i = 0;
          gen.clear();
          // std::exit(1);
        }

        if (height_gen < pdf_gen)
          event.SetGenPDF(pdf_gen);
        else
          continue;
      }

      gen.push_back(event);

      break;
    }
  }

  const std::string fname = outFile.size() > 0 ? outFile : "dat/bdstokppimkmpip-" + type + "-" + std::to_string(seed) + ".root";

  SaveData(fname, gen);

  return 0;
}
