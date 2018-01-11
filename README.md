# RObustly Complete control Synthesis (ROCS)

**ROCS** is a robust symbolic control synthesis tool for nonlinear systems. It provides two control synthesis methods from linear temporal logic (LTL) formulas:
- Abstraction-based control synthesis (to be added): state and control space are uniformly discretized.
- Specification-guided control synthesis: state space are partitioned adaptively with respect to given specifications and local dynamics of the systems.

## Installation
### Prerequisite packages and softwares
 - boost c++ libraries (http://www.boost.org): libraries for graph searching, ode solvers;
 - armadillo (http://arma.sourceforge.net): a c++ linear algebra library;
 - matlab: control simulation and graphics.
### Configuration
 - Update the environment variable DYLD_LIBRARY_PATH in terminal:
```
export DYLD_LIBRARY_PATH=“/Applications/MATLAB_R2016a.app/bin/maci64:/Applications/MATLAB_R2016a.app/sys/os/maci64:$DYLD_LIBRARY_PATH”
```


## Usage
- Clone this repositary. Source code will be downloaded and ready for use.
- Make sure the prerequisites are installed.
- Go to the folder *examples*, modify the *Makefile* according to your actual installed paths of prerequisites, e.g.
```
BOOSTROOT = /usr/local/Cellar/boost/1.62.0
MATLABROOT = /Applications/MATLAB_R2016a.app
LINALGROOT = /usr/local/Cellar/armadillo/7.500.2
```
- Type `make` in command window to compile examples. For example, for dcdc example, type `make dcdc` in command window.
- Execute compiled examples by typing the name of the executable, e.g. `./dcdc`.
- Run *simulation.m* in the same example folder using matlab.

## Examples
Simulations of the examples can be found under the folder *examples*:
- *examples/dcdc*: DCDC converter invariance control.
- *examples/car*: Robot car reach avoid control: reach a target region while avoiding obstacles.
- *examples/ipdl*: Inverted pendulum reach stay control: reach and stay in the upright position.
- *examples/poly*: Polynomial systems region of attraction computation.
- *examples/temperature*: Room temperature control (4-mode system): reach and stay in desired room and heater temperatures.

## Create your own examples
To create your own examples, you have to provide a *main file* specifying all necessary problem settings and an additional file of system dynamics. For the DCDC converter example,
- main file: *DemoDcdc.cpp*
- vectorfield file: *dcdc.hpp*

## Features included
- Basic invariance and reachability control for discrete-time systems.
- Buchi, coBuchi and reach-stay (a simpler version of coBuchi) control synthesis algorithms.

## Features under developement
- Control synthesis from more complex LTL formulas.
- Integrate continuous-time validated ode solver.


## Note
- Only tested on Mac OSX 10.11/12.