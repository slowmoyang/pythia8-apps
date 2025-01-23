#include "Pythia8/Pythia.h"
#include "Pythia8/LHEF3.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8Plugins/HepMC3.h"
#include "Pythia8Plugins/JetMatching.h"
#include "Pythia8Plugins/aMCatNLOHooks.h"

#include "TFile.h"
#include "TTree.h"

#include "argparse/argparse.hpp"

#include <algorithm>
#include <memory>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <filesystem>
#include <cmath>

namespace fs = std::filesystem;


// Returns the number of events in the input LHE file
// Read each event and write them out again, also in reclustered form.
long countEvents(const std::string& input_file_path) {
  Pythia8::Reader reader{input_file_path};
  long num_events = 0;
  while (reader.readEvent()) {
    ++num_events;
  }
  return num_events;
}


// Get the inclusive x-section by summing over all process x-sections.
double getMG5AMCNLOCrossSection(Pythia8::Pythia& pythia) {
  double xs = 0.;
  for (int idx = 0; idx < pythia.info.nProcessesLHEF(); ++idx) {
    xs += pythia.info.sigmaLHEF(idx);
  }
  return xs;
}


// Returns a weight normalisation factor for Pythia8::Pythia8.info::accumulateXsec
// taken from https://github.com/HEPcodes/MG5aMC_PY8_interface/blob/c1b01c2721826ccd6eb0f583d9b6ea77df1a320e/MG5aMC_PY8_interface.cc
double getEventWeightNormalizer(
    Pythia8::Pythia& pythia,
    const double mg5amcnlo_xsecec_pb,
    const int num_events) {

  // https://pythia.org/latest-manual/LHA.html
  const bool is_weight_in_pb = std::abs(pythia.info.lhaStrategy()) == 4;
  const double norm_den = is_weight_in_pb ? 1.0 : mg5amcnlo_xsecec_pb;
  const double norm = (norm_den * Pythia8::PB2MB) / static_cast<double>(num_events);

  return norm;
}


// FIXME:
int getInternalMergingScheme(Pythia8::Pythia& pythia) {

  if (pythia.settings.flag("Merging:doUMEPSTree") or pythia.settings.flag("Merging:doUMEPSSubt")) {
    return 1;

  } else {
    bool do_unlops = false;
    do_unlops |= pythia.settings.flag("Merging:doUNLOPSTree");
    do_unlops |= pythia.settings.flag("Merging:doUNLOPSSubt");
    do_unlops |= pythia.settings.flag("Merging:doUNLOPSLoop");
    do_unlops |= pythia.settings.flag("Merging:doUNLOPSSubtNLO");

    if (do_unlops) {
      return 2;

    } else {
      return 0;

    }
  }
}


// Additional PDF/alphaS weight for internal merging.
// Additional weight due to random choice of reclustered/non-reclustered
// treatment. Also contains additional sign for subtractive samples.
double getMergingWeight(
    Pythia8::Pythia& pythia,
    std::shared_ptr<Pythia8::amcnlo_unitarised_interface> merging_hook) {
  double merging_weight = pythia.info.mergingWeightNLO();
  if (merging_hook.get()) {
    merging_weight = merging_hook->getNormFactor();
  }
  return merging_weight;
}


void run(
    const fs::path input_file_path,
    const fs::path card_file_path,
    const fs::path output_file_path,
    const int max_events,
    const bool store_djr
) {
  if (not fs::exists(input_file_path)) {
    throw std::runtime_error(std::string{"file not found: "} + input_file_path.c_str());
  }

  if (not fs::exists(card_file_path)) {
    throw std::runtime_error(std::string{"file not found: "} + card_file_path.c_str());
  }

  if (fs::exists(output_file_path)) {
    throw std::runtime_error(std::string{"file already exists: "} + output_file_path.c_str());
  }

  /////////////////////////////////////////////////////////////////////////////
  // Pythia
  /////////////////////////////////////////////////////////////////////////////
  Pythia8::Pythia pythia{};

  pythia.readString("Beams:frameType = 4");
  pythia.settings.word("Beams:LHEF", input_file_path.c_str());
  pythia.readFile(card_file_path.c_str(), 0);

  /////////////////////////////////////////////////////////////////////////////
  // Pythia8ToHepMC
  /////////////////////////////////////////////////////////////////////////////

  Pythia8::Pythia8ToHepMC hepmc_writer{output_file_path.c_str()};
  // Switch off warnings for parton-level events.
  hepmc_writer.set_print_inconsistency(false);
  hepmc_writer.set_free_parton_warnings(false);

  /////////////////////////////////////////////////////////////////////////////
  // NOTE: Jet Matching and Merging
  //
  // adapted from main89.cc in pythia8 examples
  /////////////////////////////////////////////////////////////////////////////
  const bool do_matching = pythia.settings.flag("JetMatching:merge");
  const bool do_merging = pythia.settings.word("Merging:Process") != "void";

  if (do_matching and do_merging) {
    throw std::runtime_error("both matching and merging are used");
  }

  std::shared_ptr<Pythia8::JetMatchingMadgraph> matching_hook = nullptr;
  if (do_matching) {
    matching_hook = std::make_shared<Pythia8::JetMatchingMadgraph>();
    pythia.setUserHooksPtr(matching_hook);
  }

  // NOTE: djr
  std::shared_ptr<TFile> djr_file = nullptr;
  std::shared_ptr<TTree> djr_tree = nullptr;
  std::vector<double> djr_vec;
  int djr_num_jets;
  double djr_weight;

  if (do_matching and store_djr) {
    const std::string djr_file_name = output_file_path.stem().c_str() + std::string{"_djr.root"};
    const fs::path djr_file_path = output_file_path.parent_path() / djr_file_name;

    djr_file = std::make_shared<TFile>(djr_file_path.c_str(), "RECREATE");
    djr_tree = std::make_shared<TTree>("tree", "tree");
    djr_tree->SetDirectory(djr_file.get());

    djr_tree->Branch("djr", "vector<double>", &djr_vec);
    djr_tree->Branch("weight", &djr_weight);
    djr_tree->Branch("num_jets", &djr_num_jets);
  }

  std::shared_ptr<Pythia8::amcnlo_unitarised_interface> merging_hook = nullptr;
  if (do_merging) {
    // Allow to set the number of addtional partons dynamically.
    // Store merging scheme.
    const int scheme = getInternalMergingScheme(pythia);
    merging_hook = std::make_shared<Pythia8::amcnlo_unitarised_interface>(scheme);
    pythia.setUserHooksPtr(merging_hook);

  }

  /////////////////////////////////////////////////////////////////////////////
  //
  /////////////////////////////////////////////////////////////////////////////

  pythia.init();

  /////////////////////////////////////////////////////////////////////////////
  //
  /////////////////////////////////////////////////////////////////////////////
  const double mg5amcnlo_xsec = getMG5AMCNLOCrossSection(pythia);
  // FIXME:
  const long num_events = (max_events >= 0) ? max_events : countEvents(input_file_path);

  const auto fail_tolerance = static_cast<uint32_t>(pythia.mode("Main:timesAllowErrors"));
  uint32_t fail_count = 0;

  double xsec = 0.;
  double xsec_error_squared = 0.;

  std::cout << "Start generating events" << std::endl;

  while (pythia.info.nSelected() < num_events) {
    const bool generation_okay = pythia.next();

    if (not generation_okay) {
      if (pythia.info.atEndOfFile()) {
        std::cout << "reached the end of the file" << std::endl;
        break;

      } else {
        if (++fail_count >= fail_tolerance) {
          throw std::runtime_error("aborted because of too many failures");
        }

        continue;
      }
    }

    ///////////////////////////////////////////////////////////////////////////
    //
    ///////////////////////////////////////////////////////////////////////////

    // FIXME: rename
    const double central_weight = pythia.info.weight();
    const double merging_weight = getMergingWeight(pythia, merging_hook);
    const double weight = central_weight * merging_weight;
    const double weight_norm = getEventWeightNormalizer(pythia, mg5amcnlo_xsec, num_events);

    xsec += weight * weight_norm;
    xsec_error_squared += std::pow(weight * weight_norm, 2);

    // protect against 0-weight from internal merging
    if (std::abs(merging_weight) == 0) {
      continue;
    }

    ///////////////////////////////////////////////////////////////////////////
    // write
    ///////////////////////////////////////////////////////////////////////////

    // hepmc
    hepmc_writer.writeNextEvent(pythia);

    // djr
    if (do_matching and store_djr) {
      const std::vector<double> djr_vec_src = matching_hook->getDJR();
      djr_num_jets = matching_hook->nMEpartons().first;
      djr_vec.clear();
      djr_vec.resize(static_cast<int>(djr_vec_src.size()));
      std::copy(djr_vec_src.begin(), djr_vec_src.end(), djr_vec.begin());
      djr_weight = weight;
      djr_tree->Fill();
    }
  } // generation loop

  // NOTE: finalize
  if (do_matching and store_djr) {
    djr_file->Write();
    djr_file->Close();
  }

  // NOTE: report
  const double xsec_error = std::sqrt(xsec_error_squared);

  std::cout << "nSelected / nTotal = "
            << pythia.info.nSelected() << " / " << num_events
            << std::endl;

  std::cout
      << "MG5_aMC@NLO's cross section: "
      << std::scientific << std::setprecision(8)
      << mg5amcnlo_xsec << " pb"
      << std::endl;

  std::cout
      << "Inclusive cross section: "
      << std::scientific << std::setprecision(8)
      << xsec * Pythia8::MB2PB
      << "  +-  "
      << xsec_error * Pythia8::MB2PB << " pb"
      << std::endl;
}


// main
int main(int argc, char* argv[]) {
  argparse::ArgumentParser parser{"pythia8-hadronizer"};

  parser.add_argument("-i", "--input")
    .required()
    .help("input LHE file");

  parser.add_argument("-c", "--card")
    .required()
    .help("pythia8 card file");

  parser.add_argument("-o", "--output")
    .required()
    .help("output HepMC3 file");

  parser.add_argument("-n", "--max-events")
    .scan<'d', long>()
    .required()
    .default_value(-1l)
    .help("if the value is negative, process all events");

  parser.add_argument("--djr")
    .flag()
    .help("store djr into a ROOT file when JetMatching is used. The root file name will be derived from the output file name.");

  try {
    parser.parse_args(argc, argv);

  } catch (const std::exception& err) {
    std::cerr << "😱😱😱: " << err.what() << std::endl;
    std::cerr << parser;
    return 1;
  }


  try {
    const fs::path input_file = parser.get<std::string>("input");
    const fs::path card_file = parser.get<std::string>("card");
    const fs::path output_file = parser.get<std::string>("output");
    const long max_events = parser.get<long>("max-events");
    const bool store_djr = parser.get<bool>("--djr");

    run(
      /*input_file=*/input_file,
      /*card_file=*/card_file,
      /*output_file=*/output_file,
      /*max_events=*/max_events,
      /*store_djr=*/store_djr
    );

  } catch (const std::exception& err) {
    std::cerr << "😱😱😱: " << err.what() << std::endl;
    return 1;
  }

  return 0;
}
