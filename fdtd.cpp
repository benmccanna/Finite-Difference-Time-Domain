#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <vector>
#include "scenario.cpp"

using namespace std;

// How often (in frames) to record the state of the simulation
const int outputInterval = 1;
const double courant = 1.0;
const int defaultSource = 150;

void simulate(Scenario& scene,
              int duration,
              int outputInterval,
              string basesname);
void outputField(vector<double>& F, ofstream& output);

int main() {
  int size, dur;

  size = 501;
  dur = 1001;

  Scenario basic;
  basic.size = size;
  simulate(basic, dur, outputInterval, "output/PEC-PMC");

  class : public TotalScattered,
          public LossLayer,
          public AbsorbingBoundaries {
  } scene1;
  scene1.size = size;
  scene1.electricLoss = 0.02;
  simulate(scene1, dur, outputInterval,
          "output/TFSF-Lossy-Dielectric-ABC");
          
  class : public TotalScattered,
          public Dielectric,
          public AdvectionABC {
  } scene2;
  scene2.size = size;
  scene2.dielectricPermittivity = 9.0;
  simulate(scene2, dur, outputInterval,
          "output/TFSF-Lossy-Dielectric-AABC");

  Standing standingWaves;
  standingWaves.size = size;
  standingWaves.sourceNode = floor(size / 2);
  double harmonic = 1.0;
  standingWaves.frequency = 0.5 * harmonic / double(size - 1);
  simulate(standingWaves, dur, outputInterval,
            "output/Harmonic-" + to_string(harmonic) + "(Frequency-" + 
            to_string(standingWaves.frequency) + ")");

}

/* Runs simulations and saves the results to files. Takes a `Scenario` which
describes how to run the simulation; `duration`, the number of steps to
simulate; `outputInterval`, the frequency with which to write the state of
the simulation to our output files; and `basename`, the base file name to
write the results of the simulation to. */
void simulate(Scenario& scene,
              int duration,
              int outputInterval,
              string basename) {
  cout << "Running simulation " << basename << "... " << endl;

  /* Initialise the scene and make sure it doesn't have state left over from a
  previous simulation */
  scene.init();

  // Open output files for the field arrays
  string filename;
  ofstream outE, outH;
  filename = basename + "-(Ey-Field).dat";
  outE.open(filename.c_str());
  filename = basename + "-(Hx-Field).dat";
  outH.open(filename.c_str());

  // The actual simulation:
  for (int tIndex = 0; tIndex < duration; tIndex++) {
    // Step the simulation forward
    scene.step(tIndex);

    // Output rows of data to our two output files every `outputInterval` steps
    if (tIndex % outputInterval == 0) {
      outputField(scene.Ey, outE);
      outputField(scene.Hx, outH);
    }
  }

  // Close our two output files
  outE.close();
  outH.close();

  cout << "Done." << endl;
}

// Output a row `F` to file `output`
void outputField(vector<double>& F, ofstream& output) {
  int size = F.size();
  for (int zIndex = 0; zIndex < size; zIndex++) {
    output << F[zIndex] << "\t";
  }
  output << endl;
}
