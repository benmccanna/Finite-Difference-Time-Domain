# Simulates Magic.

This program requires a C++11 compiler and a recent version of GNUPlot for its
foul majicks. 

## Compiling
Run `make`

## Running the simulations
Run `./fdtd.out`. This will run several different scenarios and write output
`.dat` files to `output/`.

## Plotting the results of the simulations
Run `make images`. This uses GNUPlot to draw images from the `.dat` files in
`output/`.

## How it do
We define a Scenario base class - this implements a very simple simulation,
calling into stub functions that can be overridden for more complex scenarios. 

We then have several other classes that are derived from the Scenario class,
adding additional properties to the simulation. For example,
`DielectricInterface` adds an interface from free space to a dielectric medium,
and `NaiveAbsorbingBoundaries` adds absorbing boundary conditions. When we come
to run simulations, we use anonymous classes and multiple inheritance to
compose the properties we desire, then set fields to configure the simulation
e.g:
```C++
class : public NaiveAbsorbingBoundaries,
        public DielectricInterface {} simulation;
simulation.size = 200; // Size of the space being simulated
simulation.sourceNode = 50; // excite the field at point 50
simulation.interfaceIndex = 100; // position of the dielectric medium
simulation.dielectricPermittivity = 9.0; // Permitivity of the dielectric medium
```

The 'building blocks' for scenarios and configuration properties available for
them are documented along with thier implementations in `scenario.cpp`.