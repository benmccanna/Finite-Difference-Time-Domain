#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <string>
#include <sstream>
#include <stdlib.h>

using namespace std;

const int size = 200;         /* spatial size of the simulation */
const int dur = 1000;        /* length of time to run the simulation */
const int outputInterval = 4; /* how often (in frames) to record the state of the simulation */
const double imp0 = 377.0;    /* impedance of free space */

class Scenario {
  public:
  
    double *Ex, *Hy;

    Scenario(int s = 200, int sn = 50) {
      size = s;
      sourceNode = sn;
    }
    
    bool init() {
      //Initialise fields to zero
      Ex = (double *)calloc(size, sizeof(double));
      Hy = (double *)calloc(size, sizeof(double));
  
      //Catch errors in memory allocation
      if(Ex == NULL || Hy == NULL) {
        cout << "Memory allocation failed" << endl;
        return false;
      }
      
      return true;
    }
    
    void step(double tIndex) {
      updateFields();
      Ex[sourceNode] += source(tIndex);
    }
    
    
  protected:
  
    int size, sourceNode;
  
    double source(double tIndex) {
      double srcT = tIndex - 30.;
      return exp(-srcT*srcT / 100.);
    }
    
    //Update the fields according to the Yee algorithm, ignoring the boundaries for now
    void updateFields() {
    
      for (int zIndex = 0; zIndex < size - 1; zIndex++) {
        Hy[zIndex] += (Ex[zIndex + 1] - Ex[zIndex])/imp0;
      }
      
      for (int zIndex = 1; zIndex < size; zIndex++) {
        Ex[zIndex] += (Hy[zIndex] - Hy[zIndex - 1])*imp0;
      }
    
    }
    
};

void simulate(Scenario &scene, int duration, int outputInterval, string basename);
void outputField(double *F, ofstream &output);


int main() {
  
  int sourceNode = 50;
  
  Scenario scene(size, sourceNode);
  simulate(scene, dur, outputInterval, "basic");
  
  return 0;

}

void simulate(Scenario &scene, int duration, int outputInterval, string basename) {
  
  scene.init();
  string filename;
  
  ofstream outE, outH;
  
  filename = basename + "-Ex.dat";
  outE.open(filename.c_str());
  
  filename = basename + "-Hy.dat";
  outH.open(filename.c_str());
  
  for(int tIndex = 0; tIndex < duration; tIndex++) {
    scene.step(tIndex);
    
    if(tIndex % outputInterval == 0) {
      outputField(scene.Ex, outE);
      outputField(scene.Hy, outH);
    }
  }
  
  outE.close();
  outH.close();
}

void outputField(double *F, ofstream &output) {

  for (int zIndex = 0; zIndex < size; zIndex++) {
    output << F[zIndex] << "\t";
  }
  output << endl;
  
}
