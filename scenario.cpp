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
    
    // Field excited at specified node
    virtual void source(double tIndex) {
      double srcT = tIndex - 30.0;
      Ex[sourceNode] += exp(-srcT*srcT / 100.0);
    }
    
    /* The step function moves the simulation forward one step. */
    void step(double tIndex) {
      updateFields(tIndex);
      source(tIndex);
    }

  protected:
  
    /* Update the interior fields according to the Yee algorithm */
    
    void updateHy() {
      double coeff, loss;

      /* Magnetic field constant at rightmost node */
      for(int zIndex = 0; zIndex < size - 1; zIndex++) {
        loss = hLoss(zIndex);
        Hy[zIndex] *= (1.0 - loss) / (1.0 + loss);
        coeff = cour / imp0 / permeability(zIndex) / (1.0 + loss);
        Hy[zIndex] += coeff * (Ex[zIndex + 1] - Ex[zIndex]);
      }
    }
    
    void updateEx() {
      double coeff, loss;

      /* Electric field constant at leftmost node */
      for(int zIndex = 1; zIndex < size; zIndex++) {
        loss = eLoss(zIndex);
        Ex[zIndex] *= (1.0 - loss) / (1.0 + loss);
        coeff = cour * imp0 / permittivity(zIndex) / (1.0 + loss);
        Ex[zIndex] += coeff * (Hy[zIndex] - Hy[zIndex - 1]);
      }
    }
    
    virtual void abcValues() {}
    virtual void tfsfHy(double tIndex) {}
    virtual void abcHy() {}
    virtual void tfsfEx(double tIndex) {}
    virtual void abcEx() {}

  
    /* This function will be used to apply different boundary conditions. In the
    basic case one field is fixed at zero at each end, and the simulated region
    behaves like a resonant cavity */
    virtual void updateFields(double tIndex) {
      abcValues();
      updateHy();
      tfsfHy(tIndex);
      abcHy();
      updateEx();
      tfsfEx(tIndex);
      abcEx();
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

/* Adds an interface between free space and a dielectric medium of a different
relative permittivity and permeability */
class DielectricInterface: virtual public Scenario{
  
  public:
    int dielectricIndex = size/2;
    double dielectricPermittivity = 2.0;
    double dielectricPermeability = 1.0;
  
  protected:
    double permittivity(double zIndex) {
      if(zIndex < dielectricIndex) {
        return 1.0;
      }
      return dielectricPermittivity;
    }
    
    double permeability(double zIndex) {
      if(zIndex < dielectricIndex) {
        return 1.0;
      }
      if(zIndex == dielectricIndex) {
        return 0.5 * (1.0 + dielectricPermeability);
      }
      return dielectricPermeability;
    }

};

/* Adds a lossy region with non-zero electric and/or magnetic conductivity */
class LossyInterface: virtual public Scenario{
  
  public:
    int lossyIndex = size/2;
    double electricLoss = 0.01;
    double magneticLoss = 0.0;
    
  protected:
    double eLoss(double zIndex) {
      if(zIndex < lossyIndex) {
        return 0.0;
      }
      return electricLoss;
    }
    
    double hLoss(double zIndex) {
      if(zIndex < lossyIndex) {
        return 0.0;
      }
      if(zIndex == lossyIndex) {
        return 0.5*electricLoss;
      }
      return electricLoss;
    }
  
};

/* Adds the simplest form of Absorbing boundary conditions (ABCs), allowing
energy to leave the simulation */
class AbsorbingBoundaries: virtual public Scenario {
  /* The previously fixed field values are now assigned to the past value of
  their nearest neighbour, eliminating reflection from the boundaries in
  free space scenarios */
  protected:
    virtual void abcValues() {
      HyOldRight = Hy[size - 2];
      ExOldLeft = Ex[1];
    }
    virtual void abcHy() {
      Hy[size - 1] = HyOldRight;
    }
    virtual void abcEx() {
      Ex[0] = ExOldLeft;
    }

    double HyOldRight, ExOldLeft;

};

/* Adds ABCs using a first order discretised version of the advection equation,
letting less of the field reflect than in the naive case */
class AdvectionABC: virtual public AbsorbingBoundaries {
  public:
    void init() {
      double temp;
      Scenario::init();
      
      temp = sqrt(permeability(0) * permittivity(0)) / cour;
      coeffLeft = (1.0 - temp) / (1.0 + temp);
      
      temp = sqrt(permeability(size - 1) * permittivity(size - 1)) / cour;
      coeffRight = (1.0 - temp) / (1.0 + temp);
    }

  /* Fields on the boundary are now updated according to their present values,
  the local relative permittivity and permeability, and the past and present
  values of their nearest neighbours */
  protected:
    void abcHy() {
      Hy[size - 1] = HyOldRight + coeffRight * (Hy[size - 2] - Hy[size - 1]);
    }
    void abcEx() {
      Ex[0] = ExOldLeft + coeffLeft * (Ex[1] - Ex[0]);
    }

    double coeffLeft, coeffRight;

};

/* Replaces the source in the basic setup with an incident wave travelling
rightwards from a point using a Total-Field/Scattered-Field boundary at that
point */
class TotalScattered: virtual public Scenario {
  public:
    int tfsfIndex = int(size/4);
    
    /* Overwriting the basic additive source node */
    void source(double zIndex) {}
    
    double tfsfSource(double tIndex) {
      double srcT = tIndex - 30.0;
      return exp(-srcT*srcT / 100.0);
    }
    
    
  /* Corrections to the update equation to inject a unidirectional wave pulse */
  protected:
    void tfsfHy(double tIndex) {
      Hy[tfsfIndex - 1] -= tfsfSource(tIndex) / imp0;
    }
    
    void tfsfEx(double tIndex) {
      Ex[tfsfIndex] += tfsfSource(tIndex + 1.0);
    }
  
};