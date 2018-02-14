# RObustly Complete control Synthesis (ROCS) #

**ROCS** is a robust symbolic control synthesis tool for nonlinear systems. It provides two control synthesis methods from linear temporal logic (LTL) formulas:
- Abstraction-based control synthesis (to be added): state and control space are uniformly discretized.
- Specification-guided control synthesis: state space are partitioned adaptively with respect to given specifications and local dynamics of the systems.

**ROCS** is currently a package of c++ APIs, which only need to be included and no recompiliation required.



## Installation ##
### Prerequisite packages and softwares ###
 - boost c++ libraries (http://www.boost.org): libraries for graph searching, ode solvers;
 - armadillo (http://arma.sourceforge.net): a c++ linear algebra library;
 - matlab: control simulation and graphics.

### Configuration of Matlab ###
 - Matlab 2017a and earlier versions are supported.
 - Follow the instruction Matlab documentation on **Set run-time library path on Mac and Linux systems**, e.g.,
   - On Mac OSX: update the environment variable DYLD_LIBRARY_PATH by placing the following commands in a startup script (~/.profile):
```
export DYLD_LIBRARY_PATH=“matlabroot/bin/maci64:matlabroot/sys/os/maci64:$DYLD_LIBRARY_PATH”
```
   - On Linux:
```
export LD_LIBRARY_PATH=“matlabroot/bin/glnxa64:matlabroot/sys/os/glnxa64:$LD_LIBRARY_PATH”
```

### Installation of **Armadillo** ###
The installation differs on Mac OSX and Linux.
 - For Mac users, recommend to use *Homebrew* with simple command in terminal:
 ```
 $brew install armadillo
 ```
 - For Linux users, packages such as *cmake*, *openBLAS*, *LAPACK* and *ARPACK*. Follow the instruction on http://arma.sourceforge.net/download.html.



## Usage ##
- Clone this repositary.
- Make sure the prerequisites are installed.
- Prepare the configuration file for build:
  + The default configuration file *configMake* in the root path is for Mac OSX, the user only need to modify the paths for the external packages, e.g.,
  ```
  BSTDIR := /usr/local/opt/boost
  MATDIR := /Applications/MATLAB_R2016a.app
  LADIR := /usr/local/opt/armadillo
  ```
  + for Linux users, rename the file *Make-linux* to *configMake*. If all the prerequisites are installed in the default paths, only modify the matlab path, e.g., 
  ```
  MATDIR = /usr/local/MATLAB/R2016b
  ```
- To run examples, go to the folder *examples*:
  + To build all the examples: type `$ make examples` in terminal.
  + To build one of the examples: go to the corresponding subfolder and type `$ make` in terminal.
- Execute compiled examples by typing the name of the generated executable, e.g. `./dcdc`. **Data will be generated in each subfolder corresponding to the specific example**.
- Run *simulation.m* in the same example folder using matlab.



## Examples ##
Simulations of the examples can be found under the folder *examples*:
- *examples/dcdc*: DCDC converter invariance control.
- *examples/car*: Robot car reach avoid control: reach a target region while avoiding obstacles.
- *examples/ipdl*: Inverted pendulum reach stay control: reach and stay in the upright position.
- *examples/temp*: Room temperature control (4-mode system): reach and stay in desired room and heater temperatures.



## Create your own examples ##
To create your own examples, you have to provide a *main file* specifying all necessary problem settings and an additional file of system dynamics. For the DCDC converter example,
- main file: *DemoDcdc.cpp*
- vectorfield file: *dcdc.hpp*



## Features included ##
- Basic invariance and reachability control for discrete-time systems.
- Buchi, coBuchi and reach-stay (a simpler version of coBuchi) control synthesis algorithms.



## Features under developement ##
- Control synthesis from more complex LTL formulas.
- Integrate continuous-time validated ode solver.



## Note ##
Tested on
- Mac OSX 10.11/12.
- Ubuntu 16.04.3 LTS.