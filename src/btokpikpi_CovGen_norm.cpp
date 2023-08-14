#include <boost/program_options.hpp>
#include <filesystem>
#include <iostream>

#include "Data.h"
#include "FitParameters.h"
#include "GenUtil.h"
#include "Par.h"
#include "Selection.h"

namespace po = boost::program_options;

int main(const int argc, const char *argv[]) {
  std::string outFile;
  unsigned int seed, nevents;

  po::options_description desc{"Options"};

  desc.add_options()("help,h", "Display usage")(
      "output-file,o", po::value<std::string>(&outFile), "Output file")(
      "seed,s", po::value<unsigned int>(&seed)->default_value(123),
      "Seed for RNG")("nevents,n",
                      po::value<unsigned int>(&nevents)->default_value(300),
                      "Number of events to be generated");

  po::variables_map args;
  po::store(po::parse_command_line(argc, argv, desc), args);
  po::notify(args);

  if (args.count("help")) {
    std::cout << desc << std::endl;
    std::exit(0);
  }

  TLorentzVector P(0.0, 0.0, 0.0, mass_Bs);

  const std::vector<double> masses = {mass_Kp, mass_pip, mass_Kp, mass_pip};

  TGenPhaseSpace generator;
  gRandom->SetSeed(seed);
  generator.SetDecay(P, masses.size(), masses.data());

  std::vector<Event> gen;
  gen.reserve(nevents);
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
