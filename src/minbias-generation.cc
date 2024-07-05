// adapted from main41.cc
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC3.h"
#include <iostream>
#include <string>


int main(int argc, char* argv[]) {
  if (argc != 4) {
    std::cout << "expected 3 arguments but got " << argc << std::endl;
    std::cerr << "Usage: generateMinBias CARD_FILE SEED OUTPUT_FILE" << std::endl;
    return 1;
  }

  const std::string card_file{argv[1]};
  const std::string seed{argv[2]};
  const std::string output_file{argv[3]};

  std::cout << "card_file: " << card_file << std::endl
            << "seed: " << seed << std::endl
            << "output_file: " << output_file << std::endl;

  Pythia8::Pythia pythia{};
  pythia.readFile(card_file);
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = " + seed);
  pythia.init();

  Pythia8::Pythia8ToHepMC hepmc3_writer{output_file};

  for (int event = 0; event < pythia.mode("Main:numberOfEvents"); ++event) {
    if (not pythia.next()) {
      continue;
    }
    hepmc3_writer.writeNextEvent(pythia);
  }
  pythia.stat();

  return 0;
}
