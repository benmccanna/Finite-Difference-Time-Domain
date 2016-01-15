# Simulates Electromagnetic fields in 1 dimension using the Finite-Difference
Time-Domain Method.

This program requires a C++11 compiler and a recent version of GNUPlot. 

## Compiling
Run `make`

## Running the simulations
Run `./fdtd.out`. This will run several different scenarios (depending on what
is written in 'fdtd.cpp' and write output `.dat` files to `output/`.

## Plotting the results of the simulations
Run `make images`. This uses GNUPlot to draw images from the `.dat` files in
`output/`.

## Structure of the implementation
We define a Scenario base class - this defines the basic algorithm, the simplest
boundary conditions, the default source function, along with many stub functions
that can be overridden for the more complex scenarios. 

We then have several other classes that are derived from the Scenario class,
adding additional properties to the simulation. For example,
`Dielectric` adds a layer of a dielectric medium (leaving the rest as free space)
and `AbsorbingBoundaries` adds the simplest absorbing boundary conditions. When
we come to run simulations, we use anonymous classes and multiple inheritance to
compose the properties we desire, then set fields to configure the simulation
e.g:
```C++
class : public AbsorbingBoundaries,
        public Dielectric {} simulation;
simulation.size = 200; // Size of the space being simulated
simulation.sourceNode = 50; // excite the field at point 50
simulation.dLeft = 100; // position of the left edge of the dielectric medium
simulation.dRight = 200; // position of the right edge of the dielectric medium
simulation.dielectricPermittivity = 9.0; // Permitivity of the dielectric medium
```

The 'building blocks' for scenarios and configuration properties available for
them are documented along with thier implementations in `scenario.cpp`.