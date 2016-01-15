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

  size = 201;

  // Scenario basic;
  // Scenario basic;
  // basic.size = size;
  // basic.sourceNode = size / 4;
  // simulate(basic, dur, outputInterval, "output/basic");

  size = 401;
  dur = 4000;

  Standing standingWaves;
  standingWaves.size = size;
  standingWaves.sourceNode = floor(size / 2);
  int harmonic = 1;
  double precis = 0.0004;
  double f0 = 0.5 * harmonic / double(size - 1);
  double freq = f0 - precis;

  for (int i=0; i < 3; i++) {
    standingWaves.frequency = freq;
    simulate(standingWaves, dur, outputInterval,
             "output/PEC-Boundaries-(Frequency-" + to_string(freq) + ")");
    freq += precis;
  }

  size = 201;
  dur = 500;

//   class : public TotalScattered, public Dielectric, public AdvectionABC {
//   } incident;
//   incident.size = size;
//   incident.sourceNode = size / 4;
//   incident.dLeft = 100;
//   incident.dRight = size - 1;
//   incident.dielectricPermittivity = 9.0;
//   simulate(incident, dur, outputInterval,
//           "output/Incident-on-perfect-dielectric");

//   class : public LossLayer,
//           public TotalScattered,
//           public Dielectric,
//           public AdvectionABC {
//   } lossyDielectric;
//   lossyDielectric.size = size;
//   lossyDielectric.sourceNode = size / 4;
//   lossyDielectric.dLeft = 100;
//   lossyDielectric.dRight = size - 1;
//   lossyDielectric.lLeft = 180;
//   lossyDielectric.lRight = size - 1;
//   lossyDielectric.dielectricPermittivity = 3.0;
//   lossyDielectric.dielectricPermeability = 3.0;
//   lossyDielectric.electricLoss = 0.02;
//   lossyDielectric.magneticLoss = 0.08;
//   simulate(lossyDielectric, dur, outputInterval,
//           "output/Dielectric-with-lossy-region");

//   lossyDielectric.magneticLoss = 0.02;
//   simulate(lossyDielectric, dur, outputInterval,
//           "output/Dielectric-with-lossy-region-matched");
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
