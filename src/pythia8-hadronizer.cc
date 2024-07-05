#include "Pythia8/Pythia.h"
#include "Pythia8/LHEF3.h"
#include "Pythia8Plugins/HepMC3.h"
#include "Pythia8Plugins/JetMatching.h"
#include "TFile.h"
#include "TTree.h"
#include <algorithm>
#include <memory>


const double pb_to_mb = 1e-9;
const double mb_to_pb = 1e9;

// Returns the number of events in the input LHE file
long countEvents(const std::string& input_file) {
  Pythia8::Reader reader{input_file};
  long num_events = 0;
  // Read each event and write them out again, also in reclustered form.
  while (reader.readEvent()) {
    ++num_events;
  }
  return num_events;
}


// Returns an inclusive cross section
double getMG5CrossSection(Pythia8::Pythia& pythia) {
  double xsec = 0.;
  for (int idx = 0; idx < pythia.info.nProcessesLHEF(); ++idx) {
    xsec += pythia.info.sigmaLHEF(idx);
  }
  return xsec;
}


// Returns a weight normalisation factor for Pythia8::Pythia8.info::accumulateXsec
double getWeightNorm(Pythia8::Pythia& pythia,
                     const double mg5_xsec_pb,
                     const int num_events) {
  const int abs_strategy = std::abs(pythia.info.lhaStrategy());
  if ((abs_strategy < 1) or (abs_strategy > 4)) {
    std::cerr << "" << std::endl;
    return 0.0;
  }

  // https://pythia.org/latest-manual/LHA.html
  const bool is_weight_in_pb = abs_strategy == 4;
  const double norm_den = is_weight_in_pb ? 1.0 : mg5_xsec_pb;
  const double norm = (norm_den * pb_to_mb) / static_cast<double>(num_events);

  return norm;
}

// main
int main(int argc, char* argv[]) {
  /////////////////////////////////////////////////////////////////////////////
  // Check that correct number of command-line arguments
  /////////////////////////////////////////////////////////////////////////////
  if (argc != 5) {
    std::cerr << "got wrong argc=" << argc
              << std::endl
              << "usage: " << argv[0] << " CARD_FILE INPUT_FILE OUTPUT_HEPMC_FILE OUTPUT_ROOT_FILE"
              << std::endl;
    return 1;
  }

  const std::string card_file{argv[1]};
  const std::string input_file{argv[2]};
  const std::string output_hepmc_file{argv[3]};
  const std::string output_root_file{argv[4]};

  std::cout << "- input_file: " << input_file << std::endl
            << "- card_file: " << card_file << std::endl
            << "- output_hepmc_file: " << output_hepmc_file << std::endl
            << "- output_root_file: " << output_root_file << std::endl;

  /////////////////////////////////////////////////////////////////////////////
  // Pythia
  /////////////////////////////////////////////////////////////////////////////

  Pythia8::Pythia pythia{};
  pythia.readString("Beams:frameType = 4");
  pythia.settings.word("Beams:LHEF", input_file);
  pythia.readFile(card_file, 0);
  pythia.readString("Weights:suppressAUX = on");

  /////////////////////////////////////////////////////////////////////////////
  // Pythia8ToHepMC
  /////////////////////////////////////////////////////////////////////////////

  Pythia8::Pythia8ToHepMC writer{output_hepmc_file};
  // Switch off warnings for parton-level events.
  writer.set_print_inconsistency(false);
  writer.set_free_parton_warnings(false);
  // Do not store the following information.
  writer.set_store_pdf(false);
  writer.set_store_proc(false);

  /////////////////////////////////////////////////////////////////////////////
  //
  /////////////////////////////////////////////////////////////////////////////

  auto jet_matching = std::make_shared<Pythia8::JetMatchingMadgraph>();
  if (pythia.settings.flag("JetMatching:merge")) {
    pythia.setUserHooksPtr(jet_matching);
  }

  /////////////////////////////////////////////////////////////////////////////
  //
  /////////////////////////////////////////////////////////////////////////////

  TFile djr_file{output_root_file.c_str(), "RECREATE"};
  TTree djr_tree{"tree", "tree"};
  djr_tree.SetDirectory(&djr_file);

  std::vector<double> djr_vec;
  int djr_num_jets;
  double djr_weight;
  djr_tree.Branch("djr", "vector<double>", &djr_vec);
  djr_tree.Branch("weight", &djr_weight);
  djr_tree.Branch("num_jets", &djr_num_jets);

  /////////////////////////////////////////////////////////////////////////////
  //
  /////////////////////////////////////////////////////////////////////////////

  // the number of events in the input LHE file
  const long num_events = countEvents(input_file);

  /////////////////////////////////////////////////////////////////////////////
  //
  /////////////////////////////////////////////////////////////////////////////
  std::cout << "Start generating events" << std::endl;

  pythia.init();

  // double counted cross section computed by MG5
  const double mg5_xsec_pb = getMG5CrossSection(pythia);

  const int max_abort = pythia.mode("Main:timesAllowErrors");
  int num_abort = 0;
  bool aborted = false;

  while ((pythia.info.nSelected() < num_events) and (not aborted)) {
    const bool generation_okay = pythia.next();

    if (generation_okay) {
      if (pythia.info.weightValueByIndex() == 0) {
        std::cerr << "got a zero-weight event" << std::endl;
        continue;
      }

      if (pythia.event.size() < 3) {
        std::cerr << "got a broken event: event.size()=" << pythia.event.size() << std::endl;
        continue;
      }

    } else {
      if (pythia.info.atEndOfFile()) {
        std::cout << "reached the end of file" << std::endl;
        break;

      } else {
        ++num_abort;

        if (num_abort >= max_abort) {
          std::cerr << "exceeded the maximum number of failed generations" << std::endl;
          aborted = true;
          break;

        } else {
          std::cerr << "failed to generate an event. skip this event." << std::endl;
          continue;

        }
      }
    }

    ///////////////////////////////////////////////////////////////////////////
    //
    ///////////////////////////////////////////////////////////////////////////

    // norm = (mb/pb) / N = 1e-9 / N
    // where
    //   - N: the number of events in LHA
    const double norm = getWeightNorm(pythia, mg5_xsec_pb, num_events);

    // normally 1
    const double central_weight = pythia.info.weight();

    // central_weight * norm,
    // where
    //   - central_weight: normally 1
    //   - norm: (mb/pb) / N
    const double weight = central_weight * norm;

    pythia.info.weightContainerPtr->accumulateXsec(norm);
    writer.setWeightNames(pythia.info.weightNameVector());
    writer.writeNextEvent(pythia);

    ///////////////////////////////////////////////////////////////////////////
    // djr
    ///////////////////////////////////////////////////////////////////////////

    const std::vector<double> djr_vec_src = jet_matching->getDJR();
    djr_num_jets = jet_matching->nMEpartons().first;
    djr_vec.clear();
    djr_vec.resize(static_cast<int>(djr_vec_src.size()));
    std::copy(djr_vec_src.begin(), djr_vec_src.end(), djr_vec.begin());
    djr_weight = weight;
    djr_tree.Fill();
  }

  /////////////////////////////////////////////////////////////////////////////
  //
  /////////////////////////////////////////////////////////////////////////////

  if (aborted) {
    std::cerr << " Run was not completed owing to too many aborted events" << std::endl;
    return 1;
  }

  pythia.stat();

  const double sigm_total = pythia.info.weightContainerPtr->getTotalXsec().at(0) * mb_to_pb;
  const double error_total = pythia.info.weightContainerPtr->getSampleXsecErr().at(0) * mb_to_pb;

  std::cout << "MG5's double counted cross section: " << mg5_xsec_pb << " pb" << std::endl;
  std::cout << "Inclusive cross section: " << sigm_total << "  +-  " << error_total << " pb " << std::endl;

  djr_file.Write();
  djr_file.Close();

  return 0;
}
