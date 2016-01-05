#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <vector>

using namespace std;

const int defaultSize = 200;  // Default size of the simulation
const double imp0 = 377.0;    // Impedance of free space

/* The Scenario class defines a minimal simulation, including the Yee algorithm
for updating fields in time. We inherit from this class to define more complex
simulations. */
class Scenario {
  public:
    int size = defaultSize; // Size of the space being simulated
    int sourceNode = 50;     // Location of the source node
    double cour = 1.0;      // Courant number, determines how fast waves travel
    vector<double> Ex, Hy;  // Electric and magnetic fields
    
    /* Init makes sure the scenario is ready to be run - it should also make
    sure there's no state left over from previous simulations. By default, the
    only state is in the two fields. */
    virtual void init() {
      // Initialise fields to zero
      Ex = vector<double> (size, 0.0);
      Hy = vector<double> (size, 0.0);
      
      if(sourceNode < 2 || sourceNode > size - 2) {
        cout << "Field source out of range of the simulation." << endl
          << "Please enter a valid source node (integer between 2 and "
          << size - 2 << ")" << endl;
        cin >> sourceNode;
      }
    }
    
    /* The step function moves the simulation forward one step. */
    void step(double tIndex) {
      updateFields();
      source(tIndex);
    }

  protected:
    /* Update the interior fields according to the Yee algorithm with boundary
    nodes updated according to boundary conditions specified in
    boundaryConditions() */
    virtual void updateFields() {
      //double coeff, loss;
      
      Hy.at(size - 1) = Hy.at(size - 2);
      
      for(int zIndex = 0; zIndex < size - 1; zIndex++) {
        // loss = hLoss(zIndex);
        // Hy[zIndex] *= (1.0 - loss) / (1.0 + loss);
        // coeff = cour / imp0 / permeability(zIndex) / (1.0 + loss);
        // Hy[zIndex] += coeff * (Ex[zIndex + 1] - Ex[zIndex]);
        
        Hy[zIndex] += (Ex[zIndex + 1] - Ex[zIndex]) / imp0;
      }
      
      Ex.at(0) = Ex.at(1);
      
      for(int zIndex = 1; zIndex < size; zIndex++) {
        // loss = eLoss(zIndex);
        // Ex[zIndex] *= (1.0 - loss) / (1.0 + loss);
        // coeff = cour * imp0 / permittivity(zIndex) / (1.0 + loss);
        // Ex[zIndex] += coeff * (Hy[zIndex] - Hy[zIndex - 1]);
        
        Ex[zIndex] += (Hy[zIndex] - Hy[zIndex - 1]) * imp0;
      }
      
      
      //boundaryConditions();
    }

    // Field excited at specified node
    virtual void source(double tIndex) {
      double srcT = tIndex - 30.0;
      Ex[sourceNode] += exp(-srcT*srcT / 100.0);
    }
    
    /* The left and right boundaries model perfect electric and magnetic
    conductors respectively, with the electric field fixed at zero at the
    leftmost node and the magnetic field fixed at zero at the rightmost once.
    Together this makes the simulated region act like a resonant cavity */
    virtual void boundaryConditions() {
      Ex[0] = 0.0;
      Hy[size - 1] = 0.0;
    }
    
    // Relative permittivity and permeability of material
    virtual double permittivity(double zIndex) {
      return 1.0;
    }
    virtual double permeability(double zIndex) {
      return 1.0;
    }
    
    /* Loss factor in the electric and magnetic fields due to electric and
    magnetic conductivity respectively */
    virtual double eLoss(double zIndex) {
      return 0.0;
    }
    virtual double hLoss(double zIndex) {
      return 0.0;
    }

};

/* Adds Absorbing boundary conditions (ABCs), allowing energy to leave the
simulation */
class NaiveAbsorbingBoundaries: virtual public Scenario {
  public:
  void init() {
    Scenario::init();
    HyOldRight = Hy[size - 2];
    ExOldLeft = Ex[1];
  }

  protected:
    void boundaryConditions() {
      Hy[size - 1] = HyOldRight;
      HyOldRight = Hy[size - 2];
      Ex[0] = ExOldLeft;
      ExOldLeft = Ex[1];
    }
    
  private:
    double HyOldRight, ExOldLeft;

};

/* Adds an interface between free space and a dielectric medium of a different
relative permittivity and permeability */
class DielectricInterface: virtual public Scenario{
  
  public:
    int interfaceIndex = size/2;
    double dielectricPermittivity = 1.0;
    double dielectricPermeability = 1.0;
  
  protected:
    double permittivity(double zIndex) {
      if(zIndex < interfaceIndex) {
        return 1.0;
      }
      if(zIndex == interfaceIndex) {
        return 0.5 * (1.0 + dielectricPermittivity);
      }
      return dielectricPermittivity;
    }
    
    double permeability(double zIndex) {
      if(zIndex < interfaceIndex) {
        return 1.0;
      }
      return dielectricPermeability;
    }

};

/* Adds ABCs using a first order discretised version of the advection equation,
letting less of the field reflect than in the naive case */
class AdvectionAbsorbingBoundaries1: virtual public Scenario {
  public:
    void init() {
      Scenario::init();
      HyOldRight = Hy[size - 2];
      ExOldLeft = Ex[1];
      double temp = sqrt(permeability(0) * permittivity(0)) / cour;
      cout << "temp: " << temp << endl;
      coeffLeft = (1.0 - temp) / (1.0 + temp);
      temp = sqrt(permeability(size - 1) * permittivity(size - 1)) / cour;
      cout << "temp: " << temp << endl;
      coeffRight = (1.0 - temp) / (1.0 + temp);
      
      cout << "coeffLeft: " << coeffLeft << endl;
      cout << "coeffRight: " << coeffRight << endl;
    }

  protected:
    void boundaryConditions() {
      Hy[size - 1] = HyOldRight + coeffRight * (Hy[size - 2] - Hy[size - 1]);
      HyOldRight = Hy[size - 2];
      Ex[0] = ExOldLeft + coeffLeft * (Ex[1] - Ex[0]);
      ExOldLeft = Ex[1];
    }
    
  private:
    double HyOldRight, ExOldLeft, coeffLeft, coeffRight;

};