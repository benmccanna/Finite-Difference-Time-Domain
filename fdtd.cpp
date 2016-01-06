#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <vector>
#include "scenario.cpp"

using namespace std;

const int dur = 500;         /* how many frames to run the simulation for */
const int outputInterval = 1; /* how often (in frames) to record the state of the simulation */
const double courant = 1.0;
const int defaultSource = 150;

void simulate(Scenario &scene, int duration, int outputInterval, string basename);
void outputField(vector<double> &F, ofstream &output);

int main() {
  
  //Basic example - free space resonator with one source node
  Scenario basic;
  basic.sourceNode = defaultSource;
  simulate(basic, dur, outputInterval, "output/basic");
  
  class : public TotalScattered
          , public DielectricInterface
          , public AdvectionABC
          {} incident;
  simulate(incident, dur, outputInterval, "output/tfsf-aabc-dielectric");
  
  class : public AbsorbingBoundaries
          , public DielectricInterface
          , public LossyInterface
          {} scene1;
  scene1.sourceNode = defaultSource;
  scene1.dielectricPermittivity = 9.0;
  scene1.electricLoss = 0.0;
  scene1.cour = courant;
  simulate(scene1, dur, outputInterval, "output/abc-dielectric");
  
  class : public AdvectionABC
          , public DielectricInterface
          , public LossyInterface
          {} scene2;
  scene2.sourceNode = defaultSource;
  scene2.dielectricPermittivity = 9.0;
  scene2.electricLoss = 0.0;
  scene2.cour = courant;
  simulate(scene2, dur, outputInterval, "output/aabc-dielectric");

  return 0;

}

/* Runs simulations and saves the results to files. Takes a `Scenario` which 
describes how to run the simulation; `duration`, the number of steps to
simulate; `outputInterval`, the frequency with which to write the state of
the simulation to our output files; and `basename`, the base file name to
write the results of the simulation to. */
void simulate(Scenario &scene, int duration, int outputInterval, string basename) {
  
  cout << "Running simulation " << basename  << "... " << endl;
  
  /* Initialise the scene and make sure it doesn't have state left over from a
  previous simulation */
  scene.init();
  
  // Open output files for the field arrays
  string filename;
  ofstream outE, outH;
  filename = basename + "-Ex.dat";
  outE.open(filename.c_str());
  filename = basename + "-Hy.dat";
  outH.open(filename.c_str());
  
  // The actual simulation:
  for(int tIndex = 0; tIndex < duration; tIndex++) {
    // Step the simulation forward
    scene.step(tIndex);
    
    // Output rows of data to our two output files every `outputInterval` steps
    if(tIndex % outputInterval == 0) {
      outputField(scene.Ex, outE);
      outputField(scene.Hy, outH);
    }
  }
  
  // Close our two output files
  outE.close();
  outH.close();
  
  cout << "Done." << endl;
}

// Output a row `F` to file `output`
void outputField(vector<double> &F, ofstream &output) {
  int size = F.size();
  for(int zIndex = 0; zIndex < size; zIndex++) {
    output << F[zIndex] << "\t";
  }
  output << endl;
  
}
