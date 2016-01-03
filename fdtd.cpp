#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <string>
#include <sstream>
#include <stdlib.h>

using namespace std;

const int size = 200;
const int maxT = 1000;
const int nameSize = 50;
const int outputInterval = 4;
const double imp0 = 377.0;

void updateFields(double *Ex, double *Hy);
void outputField(double *F, ofstream& output);
double sourceFunction(double t);

struct Grid {
  double *ex;
  double *hz;
  int sizeX;
  int timeIndex, maxT;
  double cour;
};

int main() {
  
  //Pointers to field arrays
  double *Ex, *Hy;
  string basename = "sim";
  stringstream filename;
  
  ofstream output;
  output.open("sim.dat");
  
  //Allocates cleared memory for the field arrays
  Ex = (double *)calloc(size, sizeof(double));
  Hy = (double *)calloc(size, sizeof(double));
  
  //Catch errors in memory allocation
  if(Ex == NULL || Hy == NULL) {
    cout << "Memory allocation failed" << endl;
    return 1;
  }
  
  for (int tIndex = 0; tIndex < maxT; tIndex++) {

    updateFields(Ex, Hy);
    Ex[50] += sourceFunction(tIndex);
    
    if(tIndex % outputInterval == 0) {
      filename << basename << tIndex << ".dat";
      outputField(Hy, output);
      filename.str(string());
    }

  }
  
  output.close();
  
  return 0;

}

void updateFields(double *Ex, double *Hy) {

  for (int zIndex = 0; zIndex < size - 1; zIndex++) {
    Hy[zIndex] += (Ex[zIndex + 1] - Ex[zIndex])/imp0;
  }
  
  for (int zIndex = 1; zIndex < size; zIndex++) {
    Ex[zIndex] += (Hy[zIndex] - Hy[zIndex - 1])*imp0;
  }

}

void outputField(double *F, ofstream& output) {
  
//  ofstream output;
//  output.open(filename.c_str());
  
  for (int zIndex = 0; zIndex < size; zIndex++) {

    output << F[zIndex] << "\t";

  }
  
  output << endl;
  
  // output.close();
  // cout << "done";
}

//Input wave for the E field
double sourceFunction(double t) {
  double srcT = t - 30.;
  return exp(-srcT*srcT / 100.);
}