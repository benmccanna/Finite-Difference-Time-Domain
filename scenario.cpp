#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <vector>

using namespace std;

const double imp0 = 377.0;  // Impedance of free space
const double pi = 4.0 * atan(1);

/* The Scenario class defines a minimal simulation, including the Yee algorithm
for updating fields in time. We inherit from this class to define more complex
simulations. */
class Scenario {
 public:
  int size;                  // Size of the space being simulated
  int sourceNode;            // Location of the source node
  double sourcePeak = 30.0;  // Time of the source peak
  double cour = 1.0;         // Courant number, determines how fast waves travel
  vector<double> Ey, Hx;     // Electric and magnetic field arrays

  /* Init makes sure the scenario is ready to be run - it should also make
  sure there's no state left over from previous simulations. By default, the
  only state is in the two fields. */
  virtual void init() {
    /* Initialise fields to zero */
    Ey = vector<double>(size, 0.0);
    Hx = vector<double>(size, 0.0);
    checkVars();
  }

  /* The step function moves the simulation forward one step. */
  void step(int tIndex) { updateFields(tIndex); }

 protected:
  virtual void checkVars() {
    /* Some boundary conditions fail if an additive source is on or next to the
    boundary */
    while (sourceNode < 2 || sourceNode > size - 2) {
      cout << "Set field source position - any integer from 2 to " << size - 2
           << endl;
      cin >> sourceNode;
    }
    checkDielectric();
    checkLossLayer();
  }
  virtual void checkDielectric() {}
  virtual void checkLossLayer() {}

  /* Update the interior fields according to the Yee algorithm */
  virtual void updateHx() {
    double coeff, loss;
    /* Magnetic field constant at rightmost node */
    for (int zIndex = 0; zIndex < size - 1; zIndex++) {
      loss = hLoss(zIndex);
      Hx[zIndex] *= (1.0 - loss) / (1.0 + loss);
      coeff = cour / imp0 / permeability(zIndex) / (1.0 + loss);
      Hx[zIndex] += coeff * (Ey[zIndex + 1] - Ey[zIndex]);
    }
  }
  virtual void updateEy() {
    double coeff, loss;
    /* Electric field constant at leftmost node */
    for (int zIndex = 1; zIndex < size; zIndex++) {
      loss = eLoss(zIndex);
      Ey[zIndex] *= (1.0 - loss) / (1.0 + loss);
      coeff = cour * imp0 / permittivity(zIndex) / (1.0 + loss);
      Ey[zIndex] += coeff * (Hx[zIndex] - Hx[zIndex - 1]);
    }
  }

  /* Field excited at specified node */
  virtual double source(int tIndex) {
    double srcT = double(tIndex) - sourcePeak;
    return exp(-srcT * srcT / 100.0);
  }
  virtual void updateSource(int tIndex) { Ey[sourceNode] += source(tIndex); }

  virtual void abcValues() {}
  virtual void tfsfCorrection(int tIndex) {}
  virtual void abcHx() {}
  virtual void abcEy() {}

  /* This function is a wrapper for all the possible field update functions
  required by various scenarios, allowing derived classes to avoid conflicting
  inheritance */
  virtual void updateFields(int tIndex) {
    abcValues();
    updateHx();
    tfsfCorrection(tIndex);
    abcHx();
    updateEy();
    abcEy();
    updateSource(tIndex);
  }

  /* Relative permittivity and permeability of material */
  virtual double permittivity(int zIndex) { return 1.0; }
  virtual double permeability(int zIndex) { return 1.0; }

  /* Loss factor in the electric and magnetic fields due to electric and
  magnetic conductivity respectively */
  virtual double eLoss(int zIndex) { return 0.0; }
  virtual double hLoss(int zIndex) { return 0.0; }
};

/* Scenario optimised for generating standing waves. This class is self
contained, it should only be used on its own */
class Standing : virtual public Scenario {
 public:
  void init() {
    Ey = vector<double>(size, 0.0);
    Hx = vector<double>(size - 1, 0.0);
    checkVars();
  }
  double frequency;
  double source(int tIndex) {
    double t = double(tIndex);
    return sin(2.0 * pi * frequency * t);
  }

 protected:
  void updateHx() {
    double coeff, loss;
    for (int zIndex = 0; zIndex < size - 1; zIndex++) {
      loss = hLoss(zIndex);
      Hx[zIndex] *= (1.0 - loss) / (1.0 + loss);
      coeff = cour / imp0 / permeability(zIndex) / (1.0 + loss);
      Hx[zIndex] += coeff * (Ey[zIndex + 1] - Ey[zIndex]);
    }
  }
  void updateEy() {
    double coeff, loss;
    /* Electric field constant at left and rightmost node */
    for (int zIndex = 1; zIndex < size - 1; zIndex++) {
      loss = eLoss(zIndex);
      Ey[zIndex] *= (1.0 - loss) / (1.0 + loss);
      coeff = cour * imp0 / permittivity(zIndex) / (1.0 + loss);
      Ey[zIndex] += coeff * (Hx[zIndex] - Hx[zIndex - 1]);
    }
  }
};

/* Adds a dielectric medium between dLeft and dRight */
class Dielectric : virtual public Scenario {
 public:
  int dLeft, dRight;
  double dielectricPermittivity = 2.0;
  double dielectricPermeability = 1.0;

 protected:
  void checkDielectric() {
    while (dLeft < 0 || dLeft > size - 1) {
      cout << "Set left edge of dielectric - any integer from 0 to " << size - 1
           << endl;
      cin >> dLeft;
    }
    while (dRight < dLeft || dRight > size - 1) {
      cout << "Set right edge of dielectric - any integer from " << dLeft
           << " to " << size - 1 << endl;
      cin >> dRight;
    }
  }

  double permittivity(int zIndex) {
    if (zIndex < dLeft || zIndex > dRight) {
      return 1.0;
    }
    return dielectricPermittivity;
  }
  double permeability(int zIndex) {
    if (zIndex < dLeft || zIndex > dRight) {
      return 1.0;
    }
    /* Averages the two permeability values on the boundaries so that the
    interfaces are aligned with a particular (magnetic field) node */
    if (zIndex == dLeft || zIndex == dRight) {
      return 0.5 * (1.0 + dielectricPermeability);
    }
    return dielectricPermeability;
  }
};

/* Adds a lossy region with non-zero electric and/or magnetic conductivity */
class LossLayer : virtual public Scenario {
 public:
  int lLeft, lRight;
  double electricLoss = 0.01;
  double magneticLoss = 0.0;

 protected:
  void checkLossLayer() {
    while (lLeft < 0 || lLeft > size - 1) {
      cout << "Set left edge of lossy region - any integer from 0 to "
           << size - 1 << endl;
      cin >> lLeft;
    }
    while (lRight < lLeft || lRight > size - 1) {
      cout << "Set right edge of lossy region - any integer from " << lLeft
           << " to " << size - 1 << endl;
      cin >> lRight;
    }
  }

  double eLoss(int zIndex) {
    if (zIndex < lLeft || zIndex > lRight) {
      return 0.0;
    }
    return electricLoss;
  }
  double hLoss(int zIndex) {
    if (zIndex < lLeft || zIndex > lRight) {
      return 0.0;
    }
    /* Again, to ensure the location of the interfaces is unambiguous */
    if (zIndex == lLeft || zIndex == lRight) {
      return 0.5 * magneticLoss;
    }
    return magneticLoss;
  }
};

/* Adds the simplest form of Absorbing boundary conditions (ABCs), allowing
energy to leave the simulation. Not compatible with AdvectionABC */
class AbsorbingBoundaries : virtual public Scenario {
  /* The previously fixed field values are now assigned to the past value of
  their nearest neighbour, eliminating reflection from the boundaries in
  free space scenarios */
 protected:
  virtual void abcValues() {
    HxOldRight = Hx[size - 2];
    EyOldLeft = Ey[1];
  }
  virtual void abcHx() { Hx[size - 1] = HxOldRight; }
  virtual void abcEy() { Ey[0] = EyOldLeft; }
  double HxOldRight, EyOldLeft;
};

/* Adds ABCs using a first order discretised version of the advection equation,
letting less of the field reflect than in the naive case. Not compatible with
AbsorbingBoundaries */
class AdvectionABC : virtual public AbsorbingBoundaries {
 public:
  void init() {
    double temp;
    Scenario::init();

    temp = sqrt(permeability(0) * permittivity(0)) / cour;
    coeffLeft = (1.0 - temp) / (1.0 + temp);

    temp = sqrt(permeability(size - 1) * permittivity(size - 1)) / cour;
    coeffRight = (1.0 - temp) / (1.0 + temp);
  }

 protected:
  /* Fields on the boundary are now updated according to their present values,
  the local relative permittivity and permeability, and the past and present
  values of their nearest neighbours */
  void abcHx() {
    Hx[size - 1] = HxOldRight + coeffRight * (Hx[size - 2] - Hx[size - 1]);
  }
  void abcEy() { Ey[0] = EyOldLeft + coeffLeft * (Ey[1] - Ey[0]); }
  double coeffLeft, coeffRight;
};

/* Replaces the source in the basic setup with an incident wave travelling
rightwards from a point using a Total-Field/Scattered-Field boundary at that
point */
class TotalScattered : virtual public Scenario {
 protected:
  /* Corrects the update equation to make the source wave unidirecitonal */
  void tfsfCorrection(int tIndex) {
    Hx[sourceNode - 1] -= source(tIndex - 1) / imp0;
  }
};